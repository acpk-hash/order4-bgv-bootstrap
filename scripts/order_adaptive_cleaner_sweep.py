#!/usr/bin/env python3
"""Order-adaptive cleaner sweep for the conservative Galois-structured paper route.

This script does not claim that higher-order cleaners are integrated into the
ciphertext pipeline.  Its purpose is narrower: for a fixed prime p and a sweep
of support radii B, quantify when higher-order subgroup structure becomes
feasible and how much projected evaluator cost it could remove relative to the
current generic and order-four paths.

For each B, the script reports:
- the current implemented order-four candidate A with A^2 = -1 mod p;
- feasibility thresholds for dyadic orders 4, 8, 16 using the signed-circulant
  model discussed in the paper;
- the best existing character-projected decomposition returned by the current
  interpolation/search code;
- a conservative planner recommendation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.append(str(HERE))
ARTIFACT_SCRIPTS = HERE.parent / "github_order4_cleaner" / "scripts"
if ARTIFACT_SCRIPTS.exists():
    sys.path.append(str(ARTIFACT_SCRIPTS))

from character_projected_cleaner_search import (  # noqa: E402
    candidate_decompositions,
    interpolate_cleaner,
    small_positive_sqrt_minus_one,
)
from sat_cleaning_circuit_search import centered_mod, ps_cost_for_degree  # noqa: E402


def dyadic_signed_circulant_feasible(p: int, B: int, order: int) -> dict:
    assert order >= 4 and order & (order - 1) == 0
    dim = order // 2
    if order == 4:
        aux = 256 if p == 65537 else None
    elif order == 8:
        aux = 16 if p == 65537 else None
    elif order == 16:
        aux = 4 if p == 65537 else None
    else:
        aux = None

    if aux is None:
        return {
            "order": order,
            "dimension": dim,
            "aux": None,
            "sum_abs_powers": None,
            "max_B_sufficient": None,
            "feasible_by_bound": False,
        }

    total = sum(abs(aux) ** j for j in range(dim))
    max_b = ((p - 1) // total - 1) // 2
    feasible = (2 * B + 1) * total < p
    return {
        "order": order,
        "dimension": dim,
        "aux": aux,
        "sum_abs_powers": total,
        "max_B_sufficient": int(max_b),
        "feasible_by_bound": feasible,
    }


def choose_recommendation(B: int, feasible: list[dict]) -> str:
    order16 = next((r for r in feasible if r["order"] == 16), None)
    order8 = next((r for r in feasible if r["order"] == 8), None)
    order4 = next((r for r in feasible if r["order"] == 4), None)
    if order16 and order16["feasible_by_bound"]:
        return "order-16 feasible in the sufficient-bound model; use only as a tiny-local-support theory candidate"
    if order8 and order8["feasible_by_bound"]:
        return "order-8 becomes feasible; this is the first non-implemented higher-order target worth planner attention"
    if order4 and order4["feasible_by_bound"]:
        return "order-4 is the stable global choice; higher orders are not support-safe under the sufficient bound"
    return "no dyadic structured cleaner is support-safe under the current sufficient bound"


def sweep_row(p: int, B: int) -> dict:
    roots = small_positive_sqrt_minus_one(p)
    if not roots:
        raise RuntimeError(f"p={p} has no sqrt(-1) candidate")
    chosen_aux = min(roots, key=lambda x: abs(centered_mod(x, p)))
    row = interpolate_cleaner(p, B, chosen_aux)
    if not row.get("valid"):
        raise RuntimeError(f"chosen aux={chosen_aux} invalid for B={B}")

    decomps = candidate_decompositions(row, 16)
    best = decomps[0] if decomps else None
    feasible = [dyadic_signed_circulant_feasible(p, B, order) for order in (4, 8, 16)]
    return {
        "B": B,
        "support_size": (2 * B + 1) ** 2,
        "chosen_order4_aux": chosen_aux,
        "chosen_order4_aux_centered": centered_mod(chosen_aux, p),
        "degree": row["degree"],
        "nonzero_terms": row["nonzero_terms"],
        "generic_ps_mults": ps_cost_for_degree(row["degree"])["mults"],
        "best_character_decomposition": best,
        "dyadic_feasibility": feasible,
        "planner_recommendation": choose_recommendation(B, feasible),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, default=65537)
    parser.add_argument("--B-list", nargs="*", type=int, default=[1, 3, 5, 7, 9, 17])
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    rows = [sweep_row(args.p, B) for B in args.B_list]
    result = {
        "p": args.p,
        "rows": rows,
        "notes": {
            "scope": "offline planner evidence only; not a ciphertext implementation claim",
            "sufficient_bound": "(2B+1) * sum_{j=0}^{d-1} |A|^j < p for the d=order/2 signed-circulant model",
            "current_positive_ciphertext_claim": "only the order-4 A=256 cleaner/evaluator is implemented and benchmarked in HElib",
        },
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    if args.json:
        print(text)
        return

    print(f"p={args.p}")
    print("B support degree terms genericPS best_mod best_mults recommendation")
    for row in rows:
        best = row["best_character_decomposition"] or {}
        print(
            f"{row['B']:>2} {row['support_size']:>7} {row['degree']:>6} {row['nonzero_terms']:>5} "
            f"{row['generic_ps_mults']:>9} {best.get('modulus', '-'):>7} {best.get('estimated_mults', '-'):>10} "
            f"{row['planner_recommendation']}"
        )


if __name__ == "__main__":
    main()
