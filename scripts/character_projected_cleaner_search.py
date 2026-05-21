#!/usr/bin/env python3
"""Character-projected search for bounded-support cleaning polynomials.

For a fixed finite field F_p, support radius B, and auxiliary radix A, the
cleaner is the unique polynomial P_A with

    P_A(h A + lo) = lo,        |h|, |lo| <= B.

This script turns the previous order-four observation into a reusable search
framework.  For each candidate A it interpolates P_A, decomposes its nonzero
monomials by residue classes modulo d,

    P_A(X) = sum_r X^r Q_r(X^d),

and estimates whether a character-projected evaluator can beat a generic
Paterson--Stockmeyer evaluation.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent))

from sat_cleaning_circuit_search import (  # noqa: E402
    centered_mod,
    degree,
    newton_interpolation_coeffs,
    newton_to_monomial,
    polynomial_eval,
    ps_cost_for_degree,
    support_table,
)
from structured_aux_cleaner_experiment import (  # noqa: E402
    divisors,
    low_order_candidates,
    multiplicative_order,
)


def inv_mod(a: int, p: int) -> int:
    return pow(a % p, p - 2, p)


def parse_int_list(items: list[str]) -> list[int]:
    out: list[int] = []
    for item in items:
        for part in item.split(","):
            part = part.strip()
            if part:
                out.append(int(part, 0))
    return out


def small_positive_sqrt_minus_one(p: int) -> list[int]:
    if p % 4 != 1:
        return []
    # For p prime and p = 1 mod 4, sqrt(-1) is a^(p-1)/4 for a quadratic
    # non-residue.  Try small bases and return both roots.
    roots: set[int] = set()
    for a in range(2, min(p, 128)):
        r = pow(a, (p - 1) // 4, p)
        if (r * r + 1) % p == 0:
            roots.add(r)
            roots.add((-r) % p)
            break
    return sorted(roots, key=lambda x: abs(centered_mod(x, p)))


def addition_chain_cost(targets: set[int]) -> int:
    """Minimal multiplications to make all small powers in targets.

    The search is intentionally tiny: we only use it for residue/order powers
    such as {3,4} or {7,8}.  It keeps evaluator estimates stable without
    importing a full addition-chain library.
    """

    targets = {t for t in targets if t > 1}
    if not targets:
        return 0
    max_t = max(targets)
    start = frozenset({1})
    frontier = {start}
    seen = {start}
    for depth in range(1, max_t + 1):
        nxt: set[frozenset[int]] = set()
        for state in frontier:
            vals = sorted(state)
            for a in vals:
                for b in vals:
                    c = a + b
                    if c > max_t or c in state:
                        continue
                    ns = frozenset(set(state) | {c})
                    if targets.issubset(ns):
                        return depth
                    if ns not in seen:
                        seen.add(ns)
                        nxt.add(ns)
        frontier = nxt
    return max_t - 1


def interpolate_cleaner(p: int, B: int, aux: int) -> dict:
    points, target, _support, collisions = support_table(p, B, aux)
    if collisions:
        return {
            "aux": aux,
            "centered_aux": centered_mod(aux, p),
            "valid": False,
            "collisions": len(collisions),
        }
    newton = newton_interpolation_coeffs(points, target, p)
    poly = newton_to_monomial(points, newton, p)
    verified = all(polynomial_eval(poly, x, p) == y % p for x, y in zip(points, target))
    if not verified:
        raise RuntimeError(f"interpolation failed for p={p} B={B} aux={aux}")
    terms = [(i, c % p) for i, c in enumerate(poly) if c % p]
    return {
        "aux": aux,
        "centered_aux": centered_mod(aux, p),
        "valid": True,
        "support": len(points),
        "degree": degree(poly),
        "nonzero_terms": len(terms),
        "terms": terms,
        "generic_ps": ps_cost_for_degree(degree(poly)),
        "multiplicative_order": multiplicative_order(aux, p),
        "sqrt_minus_one": (aux * aux + 1) % p == 0,
    }


def residue_decomposition(terms: list[tuple[int, int]], d: int) -> dict:
    classes: dict[int, list[tuple[int, int, int]]] = {r: [] for r in range(d)}
    for exp, coeff in terms:
        r = exp % d
        q_exp = (exp - r) // d
        classes[r].append((q_exp, exp, coeff))
    nonempty = {r: vals for r, vals in classes.items() if vals}
    q_degrees = {
        r: max(q_exp for q_exp, _exp, _coeff in vals) for r, vals in nonempty.items()
    }
    return {
        "modulus": d,
        "num_residue_classes": len(nonempty),
        "residue_counts": {str(r): len(vals) for r, vals in nonempty.items()},
        "q_degrees": {str(r): q_degrees[r] for r in sorted(q_degrees)},
    }


def evaluator_cost_for_decomposition(terms: list[tuple[int, int]], d: int) -> dict:
    decomp = residue_decomposition(terms, d)
    q_degrees = {int(r): deg for r, deg in decomp["q_degrees"].items()}
    residues = sorted(q_degrees)

    # Need X^d and the residue powers used to multiply nonconstant Q_r(X^d).
    power_targets = {d}
    for r, qdeg in q_degrees.items():
        if r > 1 and qdeg > 0:
            power_targets.add(r)
    power_mults = addition_chain_cost(power_targets)

    q_eval_mults = sum(ps_cost_for_degree(qdeg)["mults"] for qdeg in q_degrees.values())
    component_mults = sum(1 for r, qdeg in q_degrees.items() if r != 0 and qdeg > 0)
    total = power_mults + q_eval_mults + component_mults
    return {
        **decomp,
        "power_targets": sorted(power_targets),
        "power_precompute_mults": power_mults,
        "q_eval_mults": q_eval_mults,
        "component_mults": component_mults,
        "estimated_mults": total,
    }


def candidate_decompositions(row: dict, max_d: int) -> list[dict]:
    if not row.get("valid"):
        return []
    order = row.get("multiplicative_order") or 1
    ds = {2, 4}
    ds.update(d for d in divisors(order) if 1 < d <= max_d)
    ds = {d for d in ds if d <= max_d}
    out = []
    for d in sorted(ds):
        out.append(evaluator_cost_for_decomposition(row["terms"], d))
    out.sort(key=lambda x: (x["estimated_mults"], x["num_residue_classes"], x["modulus"]))
    return out


def default_candidates(p: int, B: int, max_order: int) -> list[int]:
    out: set[int] = {2 * B + 1, p - (2 * B + 1)}
    out.update(small_positive_sqrt_minus_one(p))
    out.update(low_order_candidates(p, max_order))
    # A few powers often worth checking around the HElib-friendly integer range.
    for a in [16, 32, 64, 128, 256, 512, 1024, 4096]:
        if 1 < a < p:
            out.add(a)
            out.add((-a) % p)
    out.discard(0)
    out.discard(1)
    return sorted(out)


def strip_terms(row: dict) -> dict:
    row = dict(row)
    row.pop("terms", None)
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--B", type=int, required=True)
    parser.add_argument("--aux", nargs="*", default=[])
    parser.add_argument("--max-order", type=int, default=64)
    parser.add_argument("--max-decomposition", type=int, default=16)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    candidates = set(default_candidates(args.p, args.B, args.max_order))
    candidates.update(a % args.p for a in parse_int_list(args.aux))
    candidates = sorted(a for a in candidates if a not in {0, 1})

    rows = []
    for aux in candidates:
        row = interpolate_cleaner(args.p, args.B, aux)
        if row.get("valid"):
            decomps = candidate_decompositions(row, args.max_decomposition)
            row["best_decomposition"] = decomps[0] if decomps else None
            row["decompositions"] = decomps[:5]
            row["score"] = (
                row["best_decomposition"]["estimated_mults"] if row["best_decomposition"] else row["generic_ps"]["mults"],
                row["nonzero_terms"],
                row["generic_ps"]["mults"],
                abs(row["centered_aux"]),
            )
        rows.append(row)

    valid = [r for r in rows if r.get("valid")]
    valid.sort(key=lambda r: r["score"])
    invalid = [r for r in rows if not r.get("valid")]

    result = {
        "p": args.p,
        "B": args.B,
        "support_size": (2 * args.B + 1) ** 2,
        "relation_A_p": {
            "field_requirement": "p should be prime or the search must work over a field/ring with inverses for interpolation",
            "sqrt_minus_one_requirement": "A^2=-1 mod p exists for prime p iff p=1 mod 4",
            "support_injectivity": "A must avoid (h-h')A+(lo-lo')=0 mod p for |h-h'|,|lo-lo'|<=2B except the zero difference",
            "helib_practical_choice": "Use the small positive representative of A; p-A has the same field action but worse HElib scaling.",
        },
        "num_candidates": len(candidates),
        "num_valid": len(valid),
        "num_invalid": len(invalid),
        "best": [strip_terms(r) for r in valid[: args.top]],
        "invalid": invalid[: min(args.top, 20)],
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    if args.json:
        print(text)
        return

    print(
        f"p={args.p} B={args.B} candidates={len(candidates)} "
        f"valid={len(valid)} invalid={len(invalid)}"
    )
    print("rank aux centered order degree terms genericPS best_d est_mults residues")
    for i, r in enumerate(valid[: args.top], start=1):
        bd = r["best_decomposition"]
        print(
            f"{i:>3} {r['aux']:>8} {r['centered_aux']:>8} "
            f"{r['multiplicative_order']:>5} {r['degree']:>6} "
            f"{r['nonzero_terms']:>5} {r['generic_ps']['mults']:>9} "
            f"{bd['modulus']:>6} {bd['estimated_mults']:>9} "
            f"{bd['residue_counts']}"
        )


if __name__ == "__main__":
    main()
