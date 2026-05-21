#!/usr/bin/env python3
"""Idea A: tower/normal-basis analysis for the HElib set-E EvalMap.

This script checks whether the current HElib thin bootstrapping linear maps
have exploitable subfield/tower structure that can reduce the number of
homomorphic automorphisms.  It is deliberately a structure-level experiment:
if the count-level result is not positive, a ciphertext prototype would not be
well motivated.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from linear_transform_factor_planner import factor, helib_bsgs_auts


@dataclass(frozen=True)
class DimensionAnalysis:
    dim: int
    cyclotomic_factor: int
    generator_order: int
    p_order_mod_factor: int
    step2_matrix_size: int
    field_degree_of_constants: int
    lives_in_proper_subfield: bool
    diagonal_offsets: int
    helib_bsgs_automorphisms: int
    tower_trace_automorphisms_if_trace_only: int
    conclusion: str


def multiplicative_order(a: int, modulus: int) -> int:
    if math.gcd(a, modulus) != 1:
        raise ValueError(f"{a} is not invertible modulo {modulus}")
    cur = 1
    for k in range(1, modulus + 1):
        cur = (cur * a) % modulus
        if cur == 1:
            return k
    raise RuntimeError(f"no order found for {a} mod {modulus}")


def euler_phi_prime(q: int) -> int:
    return q - 1


def subfield_degree_for_root(p: int, ell: int, exponent: int = 1) -> int:
    """Degree of zeta_ell^exponent over F_p for prime ell."""

    reduced_order = ell // math.gcd(ell, exponent)
    if reduced_order == 1:
        return 1
    return multiplicative_order(p, reduced_order)


def thin_step2_size(p: int, mvec: list[int], dim: int) -> int:
    """Reproduce HElib ThinEvalMap's reduced size for set E."""

    # For set E, dvec[0] = 1 and dvec[1] = ord_p(mvec[1]).
    if dim == 0:
        return euler_phi_prime(mvec[0])
    if dim == 1:
        return euler_phi_prime(mvec[1]) // multiplicative_order(p, mvec[1])
    raise ValueError(dim)


def analyze_set_e() -> dict:
    p = 65537
    mvec = [97, 523]
    gens = [48117, 5239]
    generator_orders = [96, 29]
    m = math.prod(mvec)
    ord_p_m = multiplicative_order(p, m)
    phi_m = math.prod(q - 1 for q in mvec)
    nslots = phi_m // ord_p_m

    dims: list[DimensionAnalysis] = []
    for dim, ell in enumerate(mvec):
        D = thin_step2_size(p, mvec, dim)
        p_ord = multiplicative_order(p, ell)
        cofactor = m // ell
        inflate = dim == 1
        exponent_multiplier = cofactor * (ord_p_m if inflate else 1)

        # Since ell is prime and cofactor is coprime to ell, only the optional
        # inflate-by-d factor can change the root order.  For set E it does not
        # kill the primitive 523-root because gcd(18, 523) = 1.
        const_degree = subfield_degree_for_root(p, ell, exponent_multiplier)
        proper_subfield = const_degree < ord_p_m
        bsgs = helib_bsgs_auts(D)["automorphisms"]
        trace_tower = sum(r - 1 for r in factor(D))

        if proper_subfield and D == generator_orders[dim]:
            conclusion = (
                "constants lie in a proper subfield, but the Vandermonde block "
                "has all diagonal offsets nonzero; this can reduce plaintext "
                "constant arithmetic only after a basis change, not HElib "
                "automorphism/key-switch count"
            )
        elif D < generator_orders[dim]:
            conclusion = (
                "dimension is already quotient-reduced by the slot field degree; "
                "constants still require the full slot field, so no tower saving "
                "is visible"
            )
        else:
            conclusion = "no useful proper-subfield structure detected"

        dims.append(
            DimensionAnalysis(
                dim=dim,
                cyclotomic_factor=ell,
                generator_order=generator_orders[dim],
                p_order_mod_factor=p_ord,
                step2_matrix_size=D,
                field_degree_of_constants=const_degree,
                lives_in_proper_subfield=proper_subfield,
                diagonal_offsets=D,
                helib_bsgs_automorphisms=bsgs,
                tower_trace_automorphisms_if_trace_only=trace_tower,
                conclusion=conclusion,
            )
        )

    total_bsgs = sum(d.helib_bsgs_automorphisms for d in dims)
    total_offsets = sum(d.diagonal_offsets for d in dims)
    return {
        "idea": "tower/normal-basis EvalMap sparsification",
        "parameter_set": {
            "p": p,
            "m": m,
            "mvec": mvec,
            "gens": gens,
            "generator_orders": generator_orders,
            "ord_p_m": ord_p_m,
            "phi_m": phi_m,
            "nslots": nslots,
        },
        "dimension_analysis": [asdict(d) for d in dims],
        "aggregate": {
            "total_dense_diagonal_offsets": total_offsets,
            "total_helib_bsgs_automorphisms": total_bsgs,
            "automorphism_saving_from_ideaA_alone": 0,
            "why": (
                "Idea A finds proper-subfield constants only in the D=96 block. "
                "The full Vandermonde/EvalMap block still has every diagonal "
                "offset nonzero, so a normal-basis change alone does not reduce "
                "HElib rotations/key-switches.  The optimistic tower-trace count "
                "applies to trace-like maps, not to the full invertible EvalMap."
            ),
        },
        "recommended_next_step": (
            "Do not implement Idea A as a standalone ciphertext optimization. "
            "Use its subfield observation to support a fused Good-Thomas/Rader "
            "executor for the D=96 block, where the operation count already "
            "drops at plaintext level."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    report = analyze_set_e()
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    if args.json:
        print(text)
        return

    ps = report["parameter_set"]
    print(
        f"set E: p={ps['p']} m={ps['m']} ord_p_m={ps['ord_p_m']} "
        f"nslots={ps['nslots']}"
    )
    for row in report["dimension_analysis"]:
        print(
            f"dim={row['dim']} ell={row['cyclotomic_factor']} "
            f"D={row['step2_matrix_size']} const_deg={row['field_degree_of_constants']} "
            f"proper_subfield={row['lives_in_proper_subfield']} "
            f"diag_offsets={row['diagonal_offsets']} "
            f"bsgs_auts={row['helib_bsgs_automorphisms']} "
            f"trace_only_tower_auts={row['tower_trace_automorphisms_if_trace_only']}"
        )
        print(f"  {row['conclusion']}")
    ag = report["aggregate"]
    print(
        f"aggregate: dense_offsets={ag['total_dense_diagonal_offsets']} "
        f"helib_bsgs_auts={ag['total_helib_bsgs_automorphisms']} "
        f"ideaA_aut_saving={ag['automorphism_saving_from_ideaA_alone']}"
    )
    print(report["recommended_next_step"])


if __name__ == "__main__":
    main()

