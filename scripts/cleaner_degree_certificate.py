#!/usr/bin/env python3
"""Emit a minimal-degree certificate for a bounded-support cleaner.

For an injective support S with n points over F_p, the interpolating polynomial
of degree < n is unique.  Therefore the degree of this canonical interpolant is
also the minimum possible degree of any univariate exact cleaner on S.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent))

from sat_cleaning_circuit_search import (  # noqa: E402
    centered_mod,
    degree,
    newton_interpolation_coeffs,
    newton_to_monomial,
    polynomial_eval,
    support_table,
)


def make_certificate(p: int, B: int, aux: int) -> dict:
    points, target, support, collisions = support_table(p, B, aux)
    if collisions:
        return {
            "p": p,
            "B": B,
            "aux": aux,
            "valid": False,
            "reason": "support map is not injective",
            "collisions": collisions,
        }

    newton = newton_interpolation_coeffs(points, target, p)
    poly = newton_to_monomial(points, newton, p)
    deg = degree(poly)
    nonzero = sum(1 for c in poly if c % p)
    verified = all(polynomial_eval(poly, x, p) == y % p for x, y in zip(points, target))
    lead = poly[deg] % p if deg >= 0 else 0
    centered_points = [centered_mod(x, p) for x in points]

    return {
        "p": p,
        "B": B,
        "aux": aux,
        "valid": True,
        "support_size": len(points),
        "support_formula": "(2B+1)^2",
        "support_points_are_distinct": len(points) == len(support),
        "support_min_centered": min(centered_points),
        "support_max_centered": max(centered_points),
        "canonical_interpolant": {
            "degree": deg,
            "nonzero_terms": nonzero,
            "leading_coeff_mod_p": lead,
            "leading_coeff_centered": centered_mod(lead, p),
            "verified_on_support": verified,
        },
        "minimal_degree_certificate": {
            "minimum_degree": deg,
            "field": f"F_{p}",
            "argument": (
                "There is a unique interpolating polynomial of degree < |S| "
                "on |S| distinct support points. If a lower-degree exact "
                "cleaner Q existed, then P-Q would have |S| roots and degree "
                "< |S| over a field, hence P=Q, contradicting deg(P)>deg(Q)."
            ),
            "lower_degree_exact_cleaner_exists": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--B", type=int, required=True)
    parser.add_argument("--aux", type=int, required=True)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    cert = make_certificate(args.p, args.B, args.aux)
    text = json.dumps(cert, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    if args.json:
        print(text)
        return

    if not cert["valid"]:
        print(f"invalid support for p={args.p} B={args.B} aux={args.aux}")
        return
    interp = cert["canonical_interpolant"]
    print(f"p={args.p} B={args.B} aux={args.aux}")
    print(f"support={cert['support_size']} distinct={cert['support_points_are_distinct']}")
    print(
        "canonical interpolant: "
        f"degree={interp['degree']} nonzero={interp['nonzero_terms']} "
        f"verified={interp['verified_on_support']}"
    )
    print(f"minimum exact univariate degree = {interp['degree']}")


if __name__ == "__main__":
    main()
