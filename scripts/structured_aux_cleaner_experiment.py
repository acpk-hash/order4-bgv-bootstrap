#!/usr/bin/env python3
"""Structured auxiliary-radix experiment for bounded-support cleaners."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent))

from aux_cleaning_global_sweep import evaluate_aux, tags_for_aux  # noqa: E402
from sat_cleaning_circuit_search import (  # noqa: E402
    centered_mod,
    ps_cost_for_degree,
)


def divisors(n: int) -> list[int]:
    out = []
    for i in range(1, int(math.isqrt(n)) + 1):
        if n % i == 0:
            out.append(i)
            if i * i != n:
                out.append(n // i)
    return sorted(out)


def factor_distinct(n: int) -> list[int]:
    fac = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            fac.append(d)
            while n % d == 0:
                n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        fac.append(n)
    return fac


def multiplicative_order(a: int, p: int) -> int:
    order = p - 1
    for q in factor_distinct(p - 1):
        while order % q == 0 and pow(a, order // q, p) == 1:
            order //= q
    return order


def primitive_root(p: int) -> int:
    fac = factor_distinct(p - 1)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fac):
            return g
    raise ValueError(f"no primitive root found for p={p}")


def low_order_candidates(p: int, max_order: int) -> set[int]:
    g = primitive_root(p)
    out = set()
    for d in divisors(p - 1):
        if 1 < d <= max_order:
            step = (p - 1) // d
            for k in range(d):
                out.add(pow(g, step * k, p))
    out.discard(0)
    out.discard(1)
    return out


def small_window_candidates(p: int, lo: int, hi: int) -> set[int]:
    out = set()
    for a in range(max(2, lo), min(p, hi + 1)):
        out.add(a)
        out.add((-a) % p)
        out.add(pow(a, p - 2, p))
        out.add((-pow(a, p - 2, p)) % p)
    out.discard(0)
    out.discard(1)
    return out


def polynomial_root_candidates(p: int) -> set[int]:
    # Small set of algebraic shapes that often induce symmetry in S_A,B.
    # A brute-force scan over F_p is cheap for p=65537 and keeps the script
    # dependency-free.
    polys = [
        ("x^2+1", lambda x: x * x + 1),
        ("x^2+x+1", lambda x: x * x + x + 1),
        ("x^2-x+1", lambda x: x * x - x + 1),
        ("x^4+1", lambda x: pow(x, 4, p) + 1),
        ("x^4-1", lambda x: pow(x, 4, p) - 1),
    ]
    out = set()
    tags: dict[int, list[str]] = {}
    for name, poly in polys:
        for x in range(2, p):
            if poly(x) % p == 0:
                out.add(x)
                tags.setdefault(x, []).append(name)
    return out


def parse_aux_list(items: list[str]) -> set[int]:
    out = set()
    for item in items:
        for part in item.split(","):
            part = part.strip()
            if part:
                out.add(int(part, 0))
    return out


def centered_abs_aux(aux: int, p: int) -> int:
    return abs(centered_mod(aux, p))


def add_cost_fields(row: dict, p: int) -> dict:
    if not row.get("valid"):
        return row
    ps = ps_cost_for_degree(row["degree"])
    row = dict(row)
    row["ps_cost"] = ps
    row["centered_aux_abs"] = centered_abs_aux(row["aux"], p)
    row["score"] = (
        row["nonzero_terms"],
        ps["mults"],
        row["degree"],
        row["centered_aux_abs"],
    )
    row["tags"] = tags_for_aux(row["aux"], p)
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--B", type=int, required=True)
    parser.add_argument("--max-order", type=int, default=256)
    parser.add_argument("--small-lo", type=int, default=35)
    parser.add_argument("--small-hi", type=int, default=512)
    parser.add_argument(
        "--include-small-window",
        action="store_true",
        help="Also test a, -a, a^-1, -a^-1 for small-lo <= a <= small-hi.",
    )
    parser.add_argument("--aux", nargs="*", default=[])
    parser.add_argument(
        "--candidate-cap",
        type=int,
        help="Debug guard: keep only the first N sorted candidates after construction.",
    )
    parser.add_argument("--top", type=int, default=30)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    candidates = set()
    candidates.update(parse_aux_list(args.aux))
    candidates.update(low_order_candidates(args.p, args.max_order))
    candidates.update(polynomial_root_candidates(args.p))
    if args.include_small_window:
        candidates.update(small_window_candidates(args.p, args.small_lo, args.small_hi))
    candidates = {a % args.p for a in candidates if a % args.p not in {0, 1}}
    candidates = sorted(candidates)
    if args.candidate_cap is not None:
        candidates = candidates[: args.candidate_cap]

    rows = [add_cost_fields(evaluate_aux(args.p, args.B, a), args.p) for a in candidates]
    valid = [r for r in rows if r.get("valid")]
    valid.sort(key=lambda r: r["score"])
    invalid = [r for r in rows if not r.get("valid")]

    best = valid[: args.top]
    baseline_aux = [35, 256, (-256) % args.p]
    baselines = [add_cost_fields(evaluate_aux(args.p, args.B, a), args.p) for a in baseline_aux]

    result = {
        "p": args.p,
        "B": args.B,
        "num_candidates": len(candidates),
        "num_valid": len(valid),
        "num_invalid": len(invalid),
        "candidate_sources": {
            "low_order_max": args.max_order,
            "small_window": [args.small_lo, args.small_hi] if args.include_small_window else None,
            "explicit_aux": sorted(parse_aux_list(args.aux)),
        },
        "best": best,
        "baselines": baselines,
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
    print("rank aux centered|aux| degree nonzero ps_mults tags")
    for i, r in enumerate(best, start=1):
        print(
            f"{i:>3} {r['aux']:>8} {r['centered_aux_abs']:>8} "
            f"{r['degree']:>6} {r['nonzero_terms']:>7} "
            f"{r['ps_cost']['mults']:>8} {','.join(r['tags'])}"
        )


if __name__ == "__main__":
    main()
