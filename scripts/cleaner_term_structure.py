#!/usr/bin/env python3
"""Inspect monomial support structure of bounded-support cleaner polynomials."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent))

from sat_cleaning_circuit_search import (  # noqa: E402
    newton_interpolation_coeffs,
    newton_to_monomial,
    support_table,
)


def analyze_aux(p: int, B: int, aux: int, residue_mod: int) -> dict:
    points, target, _support, collisions = support_table(p, B, aux)
    if collisions:
        return {"aux": aux, "valid": False, "collisions": len(collisions)}

    poly = newton_to_monomial(points, newton_interpolation_coeffs(points, target, p), p)
    nonzero_exponents = [i for i, c in enumerate(poly) if c % p]
    residue_counts = {
        str(r): sum(1 for i in nonzero_exponents if i % residue_mod == r)
        for r in range(residue_mod)
    }
    return {
        "aux": aux,
        "valid": True,
        "degree": max(nonzero_exponents, default=-1),
        "nonzero_terms": len(nonzero_exponents),
        "residue_mod": residue_mod,
        "residue_counts": residue_counts,
        "first_nonzero_exponents": nonzero_exponents[:32],
        "last_nonzero_exponents": nonzero_exponents[-32:],
    }


def parse_aux_list(items: list[str]) -> list[int]:
    out = []
    for item in items:
        for part in item.split(","):
            part = part.strip()
            if part:
                out.append(int(part, 0))
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--B", type=int, required=True)
    parser.add_argument("--aux", nargs="+", required=True)
    parser.add_argument("--residue-mod", type=int, default=4)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    result = {
        "p": args.p,
        "B": args.B,
        "analyses": [
            analyze_aux(args.p, args.B, aux % args.p, args.residue_mod)
            for aux in parse_aux_list(args.aux)
        ],
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    if args.json:
        print(text)
        return

    for row in result["analyses"]:
        if not row["valid"]:
            print(f"aux={row['aux']} invalid collisions={row['collisions']}")
            continue
        print(
            f"aux={row['aux']} degree={row['degree']} "
            f"terms={row['nonzero_terms']} residues={row['residue_counts']}"
        )


if __name__ == "__main__":
    main()
