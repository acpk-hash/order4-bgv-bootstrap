#!/usr/bin/env python3
"""Sweep auxiliary radices for bounded-support cleaning polynomials.

For fixed p and B, the cleaning problem is

    x = lo + A * hi  ->  lo mod p,      hi, lo in [-B, B].

When the support map is injective modulo p, there is a unique polynomial of
degree < (2B+1)^2 that agrees with this function on the support.  This script
enumerates candidate A values and reports the exact interpolation degree and
monomial sparsity.  Use --mode full for a true finite search over F_p^*, and
--mode algebraic for large p where only low-degree algebraic candidates are
needed as a cheap prefilter.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Iterable

sys.path.append(str(Path(__file__).resolve().parent))

from sat_cleaning_circuit_search import (  # noqa: E402
    degree,
    newton_interpolation_coeffs,
    newton_to_monomial,
    polynomial_eval,
    support_table,
)


def inv_mod(a: int, p: int) -> int:
    return pow(a % p, p - 2, p)


def multiplicative_order(a: int, p: int, max_order: int = 32) -> int | None:
    cur = 1
    for k in range(1, max_order + 1):
        cur = (cur * a) % p
        if cur == 1:
            return k
    return None


def tags_for_aux(a: int, p: int) -> list[str]:
    tags: list[str] = []
    if (a * a + 1) % p == 0:
        tags.append("sqrt_minus_one")
    if (a * a + a + 1) % p == 0:
        tags.append("order_3_root")
    order = multiplicative_order(a, p, max_order=16)
    if order is not None and order > 1:
        tags.append(f"order_{order}")
    if not tags:
        tags.append("generic")
    return tags


def algebraic_candidates(p: int, max_order: int) -> list[int]:
    out: set[int] = set()
    for a in range(2, p):
        order = multiplicative_order(a, p, max_order=max_order)
        if order is not None and order > 1:
            out.add(a)
        if (a * a + 1) % p == 0:
            out.add(a)
        if (a * a + a + 1) % p == 0:
            out.add(a)
    return sorted(out)


def parse_aux_list(items: Iterable[str]) -> list[int]:
    out: list[int] = []
    for item in items:
        for part in item.split(","):
            part = part.strip()
            if part:
                out.append(int(part, 0))
    return out


def evaluate_aux(p: int, B: int, aux: int) -> dict:
    points, target, _support, collisions = support_table(p, B, aux)
    if collisions:
        return {
            "aux": aux,
            "valid": False,
            "collisions": len(collisions),
            "tags": tags_for_aux(aux, p),
        }

    newton = newton_interpolation_coeffs(points, target, p)
    poly = newton_to_monomial(points, newton, p)
    deg = degree(poly)
    nonzero = sum(1 for c in poly if c % p)
    verified = all(polynomial_eval(poly, x, p) == y % p for x, y in zip(points, target))
    if not verified:
        raise RuntimeError(f"interpolation failed for p={p} B={B} aux={aux}")

    return {
        "aux": aux,
        "valid": True,
        "support": len(points),
        "degree": deg,
        "nonzero_terms": nonzero,
        "tags": tags_for_aux(aux, p),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--B", type=int, required=True)
    parser.add_argument("--mode", choices=["full", "algebraic", "list"], default="full")
    parser.add_argument("--max-order", type=int, default=16)
    parser.add_argument("--aux", nargs="*", default=[])
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.mode == "full":
        candidates = list(range(2, args.p))
    elif args.mode == "algebraic":
        candidates = algebraic_candidates(args.p, args.max_order)
    else:
        candidates = parse_aux_list(args.aux)

    rows = [evaluate_aux(args.p, args.B, a) for a in candidates]
    valid = [r for r in rows if r["valid"]]
    valid.sort(
        key=lambda r: (
            r["degree"],
            r["nonzero_terms"],
            0 if "sqrt_minus_one" in r["tags"] else 1,
            r["aux"],
        )
    )

    generic = [r for r in valid if "generic" in r["tags"]]
    structured = [r for r in valid if "generic" not in r["tags"]]
    result = {
        "p": args.p,
        "B": args.B,
        "mode": args.mode,
        "num_candidates": len(candidates),
        "num_valid": len(valid),
        "best": valid[: args.top],
        "best_generic": generic[: min(args.top, 10)],
        "best_structured": structured[: min(args.top, 10)],
    }

    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n")
    if args.json:
        print(text)
        return

    print(
        f"p={args.p} B={args.B} mode={args.mode} "
        f"candidates={len(candidates)} valid={len(valid)}"
    )
    print("best candidates: aux degree nonzero tags")
    for r in valid[: args.top]:
        print(f"{r['aux']:>8} {r['degree']:>6} {r['nonzero_terms']:>7} {','.join(r['tags'])}")
    if structured:
        print("best structured: aux degree nonzero tags")
        for r in structured[: min(args.top, 10)]:
            print(f"{r['aux']:>8} {r['degree']:>6} {r['nonzero_terms']:>7} {','.join(r['tags'])}")


if __name__ == "__main__":
    main()
