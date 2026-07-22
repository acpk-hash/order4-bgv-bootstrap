#!/usr/bin/env python3
"""SAT/SMT-style search for auxiliary-radix cleaning circuits.

The target function is the bounded overflow cleaning map

    x = hi * aux + lo  ->  lo mod p,    hi, lo in [-B, B].

The first backend is deliberately conservative: generate a library of candidate
signals, then ask Z3 to find a sparse linear combination that is exact on the
real support.  This is a useful first SAT experiment because additions and
scalar multiplications are cheap in FHE, while the generated signals can be
ranked by an estimated multiplication/depth cost.

If Z3 is unavailable, the script still computes interpolation and linear-span
certificates with modular Gaussian elimination.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class Feature:
    label: str
    values: tuple[int, ...]
    kind: str
    exponent: int | None
    depth_est: int
    mult_est: int


def centered_mod(a: int, p: int) -> int:
    a %= p
    return a - p if a > p // 2 else a


def inv_mod(a: int, p: int) -> int:
    return pow(a % p, p - 2, p)


def support_table(p: int, B: int, aux: int) -> tuple[list[int], list[int], list[dict], dict[int, list[int]]]:
    table: dict[int, int] = {}
    witness: dict[int, list[int]] = {}
    collisions: dict[int, set[int]] = {}

    for hi in range(-B, B + 1):
        for lo in range(-B, B + 1):
            x = (hi * aux + lo) % p
            y = lo % p
            witness.setdefault(x, [hi, lo])
            if x in table and table[x] != y:
                collisions.setdefault(x, {table[x]}).add(y)
            table[x] = y

    if collisions:
        return [], [], [], {x: sorted(vals) for x, vals in collisions.items()}

    points = sorted(table)
    values = [table[x] for x in points]
    support = [{"x": x, "hi": witness[x][0], "lo": witness[x][1], "y": table[x]} for x in points]
    return points, values, support, {}


def newton_interpolation_coeffs(points: list[int], values: list[int], p: int) -> list[int]:
    coeffs = [v % p for v in values]
    n = len(points)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            denom = (points[i] - points[i - j]) % p
            coeffs[i] = ((coeffs[i] - coeffs[i - 1]) * inv_mod(denom, p)) % p
    return coeffs


def newton_to_monomial(points: list[int], newton: list[int], p: int) -> list[int]:
    poly = [0]
    basis = [1]
    for i, coeff in enumerate(newton):
        if coeff:
            if len(poly) < len(basis):
                poly.extend([0] * (len(basis) - len(poly)))
            for j, b in enumerate(basis):
                poly[j] = (poly[j] + coeff * b) % p

        if i + 1 < len(newton):
            root = points[i] % p
            nxt = [0] * (len(basis) + 1)
            for j, b in enumerate(basis):
                nxt[j] = (nxt[j] - root * b) % p
                nxt[j + 1] = (nxt[j + 1] + b) % p
            basis = nxt
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def degree(poly: list[int]) -> int:
    for i in range(len(poly) - 1, -1, -1):
        if poly[i] != 0:
            return i
    return -1


def polynomial_eval(poly: list[int], x: int, p: int) -> int:
    acc = 0
    for coeff in reversed(poly):
        acc = (acc * x + coeff) % p
    return acc


def ps_cost_for_degree(d: int) -> dict[str, int]:
    if d <= 0:
        return {"degree": d, "m": 0, "k": 0, "mults": 0}
    best: tuple[int, int, int] | None = None
    for m in range(math.ceil(math.log2(d)) + 1):
        k = math.ceil(d / (2**m))
        mults = (k - 1) if m == 0 else k + m - 2
        mults += math.ceil(d / max(k, 1)) - 1
        cand = (mults, m, k)
        if best is None or cand < best:
            best = cand
    assert best is not None
    mults, m, k = best
    return {"degree": d, "m": m, "k": k, "mults": mults}


def binary_power_cost(exp: int) -> tuple[int, int]:
    if exp <= 1:
        return 0, 0
    depth = exp.bit_length() - 1
    mults = depth + max(exp.bit_count() - 1, 0)
    return depth, mults


def unique_features(features: Iterable[Feature], max_features: int | None = None) -> list[Feature]:
    seen: dict[tuple[int, ...], Feature] = {}
    for feat in features:
        old = seen.get(feat.values)
        if old is None or (feat.depth_est, feat.mult_est, len(feat.label)) < (
            old.depth_est,
            old.mult_est,
            len(old.label),
        ):
            seen[feat.values] = feat
    out = list(seen.values())
    out.sort(key=lambda f: (f.depth_est, f.mult_est, f.kind, f.exponent if f.exponent is not None else -1, f.label))
    if max_features is not None:
        return out[:max_features]
    return out


def candidate_exponents(p: int, degree_limit: int, feature_set: str, extra: list[int]) -> list[int]:
    exps: set[int] = set()
    if feature_set in {"monomial", "mixed"}:
        exps.update(range(0, degree_limit + 1))
    if feature_set in {"cheap", "mixed"}:
        exps.update({0, 1})
        k = 1
        while k <= max(p - 1, degree_limit):
            exps.add(k)
            k *= 2
        for j in range(0, degree_limit + 1):
            if p - 1 - j >= 0:
                exps.add(p - 1 - j)
        for j in range(1, degree_limit + 1):
            exps.add(j * 2)
            exps.add(j * 3)
    exps.update(e for e in extra if e >= 0)
    return sorted(exps)


def generate_features(
    p: int,
    points: list[int],
    degree_limit: int,
    feature_set: str,
    extra_exponents: list[int],
    max_features: int | None,
) -> list[Feature]:
    feats: list[Feature] = []
    for exp in candidate_exponents(p, degree_limit, feature_set, extra_exponents):
        values = tuple(pow(x, exp, p) for x in points)
        depth, mults = binary_power_cost(exp)
        feats.append(Feature(label=f"x^{exp}", values=values, kind="power", exponent=exp, depth_est=depth, mult_est=mults))
    return unique_features(feats, max_features=max_features)


def solve_linear_span(features: list[Feature], target: list[int], p: int) -> list[int] | None:
    m = len(target)
    n = len(features)
    rows = [[features[j].values[i] % p for j in range(n)] + [target[i] % p] for i in range(m)]
    pivot_cols: list[int] = []
    r = 0
    for c in range(n):
        pivot = None
        for i in range(r, m):
            if rows[i][c] % p:
                pivot = i
                break
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        inv = inv_mod(rows[r][c], p)
        rows[r] = [(v * inv) % p for v in rows[r]]
        for i in range(m):
            if i != r and rows[i][c] % p:
                factor = rows[i][c] % p
                rows[i] = [(rows[i][j] - factor * rows[r][j]) % p for j in range(n + 1)]
        pivot_cols.append(c)
        r += 1
        if r == m:
            break

    for i in range(m):
        if all(rows[i][c] % p == 0 for c in range(n)) and rows[i][-1] % p:
            return None

    coeffs = [0] * n
    for row_idx, col_idx in enumerate(pivot_cols):
        coeffs[col_idx] = rows[row_idx][-1] % p
    return coeffs


def verify_solution(features: list[Feature], coeffs: list[int], target: list[int], p: int) -> bool:
    for i, y in enumerate(target):
        acc = 0
        for coeff, feat in zip(coeffs, features):
            acc = (acc + coeff * feat.values[i]) % p
        if acc != y % p:
            return False
    return True


def z3_sparse_solve(
    features: list[Feature],
    target: list[int],
    p: int,
    max_terms: int,
    timeout_ms: int,
) -> tuple[list[int] | None, int | None, str]:
    try:
        import z3  # type: ignore
    except ImportError:
        return None, None, "z3_not_available"

    n = len(features)
    coeff = [z3.Int(f"c_{i}") for i in range(n)]
    used = [z3.Bool(f"u_{i}") for i in range(n)]

    for k in range(1, max_terms + 1):
        solver = z3.Solver()
        solver.set(timeout=timeout_ms)
        for i in range(n):
            solver.add(coeff[i] >= 0, coeff[i] < p)
            solver.add(used[i] == (coeff[i] != 0))
        solver.add(z3.Sum([z3.If(u, 1, 0) for u in used]) <= k)
        for row, y in enumerate(target):
            expr = z3.Sum([(features[i].values[row] % p) * coeff[i] for i in range(n)])
            solver.add(expr % p == y % p)
        status = solver.check()
        if status == z3.sat:
            model = solver.model()
            out = [model.eval(c, model_completion=True).as_long() % p for c in coeff]
            return out, k, "sat"
        if status == z3.unknown:
            return None, k, f"unknown:{solver.reason_unknown()}"
    return None, max_terms, "unsat"


def summarize_solution(features: list[Feature], coeffs: list[int], p: int, limit: int = 40) -> list[dict]:
    terms = []
    for coeff, feat in zip(coeffs, features):
        if coeff % p:
            terms.append(
                {
                    "coeff_mod_p": coeff % p,
                    "coeff_centered": centered_mod(coeff, p),
                    "feature": feat.label,
                    "kind": feat.kind,
                    "exponent": feat.exponent,
                    "depth_est": feat.depth_est,
                    "mult_est": feat.mult_est,
                }
            )
    terms.sort(key=lambda t: (t["depth_est"], t["mult_est"], t["feature"]))
    return terms[:limit]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--B", type=int, required=True)
    parser.add_argument("--aux", type=int, required=True)
    parser.add_argument("--degree-limit", type=int, default=64)
    parser.add_argument("--feature-set", choices=["monomial", "cheap", "mixed"], default="mixed")
    parser.add_argument("--extra-exponents", type=int, nargs="*", default=[])
    parser.add_argument("--max-features", type=int)
    parser.add_argument("--z3-sparse", action="store_true")
    parser.add_argument("--max-terms", type=int, default=12)
    parser.add_argument("--timeout-ms", type=int, default=30_000)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    points, target, support, collisions = support_table(args.p, args.B, args.aux)
    if collisions:
        result = {
            "p": args.p,
            "B": args.B,
            "aux": args.aux,
            "valid": False,
            "collisions": collisions,
        }
        text = json.dumps(result, indent=2, sort_keys=True)
        if args.output:
            args.output.write_text(text + "\n")
        print(text)
        return

    newton = newton_interpolation_coeffs(points, target, args.p)
    poly = newton_to_monomial(points, newton, args.p)
    interp_degree = max((i for i, c in enumerate(poly) if c % args.p), default=-1)
    interp_nonzero = sum(1 for c in poly if c % args.p)
    interp_ok = all(polynomial_eval(poly, x, args.p) == y % args.p for x, y in zip(points, target))

    features = generate_features(
        args.p,
        points,
        args.degree_limit,
        args.feature_set,
        args.extra_exponents,
        args.max_features,
    )
    linear_coeffs = solve_linear_span(features, target, args.p)
    linear_ok = linear_coeffs is not None and verify_solution(features, linear_coeffs, target, args.p)

    z3_coeffs = None
    z3_terms = None
    z3_status = "disabled"
    if args.z3_sparse:
        z3_coeffs, z3_terms, z3_status = z3_sparse_solve(
            features,
            target,
            args.p,
            max_terms=args.max_terms,
            timeout_ms=args.timeout_ms,
        )
        if z3_coeffs is not None and not verify_solution(features, z3_coeffs, target, args.p):
            raise RuntimeError("Z3 returned a model that failed Python verification")

    chosen = z3_coeffs if z3_coeffs is not None else linear_coeffs
    selected_terms = summarize_solution(features, chosen, args.p) if chosen is not None else []

    result = {
        "p": args.p,
        "B": args.B,
        "aux": args.aux,
        "valid": True,
        "support_size": len(points),
        "support_min_centered": centered_mod(min(points), args.p),
        "support_max_centered": centered_mod(max(points), args.p),
        "interpolation": {
            "degree": interp_degree,
            "nonzero_terms": interp_nonzero,
            "verified": interp_ok,
            "ps_cost": ps_cost_for_degree(interp_degree),
        },
        "feature_search": {
            "feature_set": args.feature_set,
            "degree_limit": args.degree_limit,
            "num_features": len(features),
            "linear_span_sat": linear_ok,
            "z3_status": z3_status,
            "z3_terms_bound": z3_terms,
            "selected_terms_count": len([c for c in (chosen or []) if c % args.p]),
            "selected_terms_preview": selected_terms,
        },
    }

    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n")
    if args.json:
        print(text)
        return

    print(f"p={args.p} B={args.B} aux={args.aux} support={len(points)} valid=True")
    print(
        "interpolation: "
        f"degree={interp_degree} nonzero={interp_nonzero} "
        f"ps_mults={result['interpolation']['ps_cost']['mults']} verified={interp_ok}"
    )
    print(
        "feature_search: "
        f"features={len(features)} span_sat={linear_ok} "
        f"z3={z3_status} selected_terms={result['feature_search']['selected_terms_count']}"
    )
    for term in selected_terms[:12]:
        print(
            f"  {term['coeff_centered']:>8} * {term['feature']:<16} "
            f"depth~{term['depth_est']} mult~{term['mult_est']}"
        )


if __name__ == "__main__":
    main()
