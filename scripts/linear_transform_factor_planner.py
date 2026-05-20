#!/usr/bin/env python3
"""Cost sketch for factorized EvalMap linear transforms.

This is not a final runtime benchmark.  It is a planning tool for deciding
which HElib bootstrapping dimensions are worth replacing by a structured
Vandermonde/CRT executor instead of the current generic diagonal MatMul1DExec.

The proposed algorithmic change is:

  current:    dense diagonal linear transform with BSGS/hoisting
  proposed:   exact factorized Vandermonde transform using radix/CRT stages

The factorized estimate counts automorphism classes for a Cooley-Tukey style
stage decomposition.  For prime lengths it also reports a Rader route that
turns a length-D transform into a length-(D-1) cyclic convolution.  The model is
deliberately simple; its role is to identify promising dimensions for a real
microbenchmark and C++ prototype.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass


@dataclass(frozen=True)
class ParamSet:
    label: str
    p: int
    mvec: tuple[int, ...]
    ords: tuple[int, ...]


PARAMS = {
    "A": ParamSet("A", 17, (29, 1321), (28, 55)),
    "B": ParamSet("B", 127, (37, 1531), (36, 34)),
    "C": ParamSet("C", 257, (43, 1289), (42, 46)),
    "D": ParamSet("D", 8191, (43, 1051), (42, 75)),
    "E": ParamSet("E", 65537, (97, 523), (96, 29)),
}


def factor(n: int) -> list[int]:
    out: list[int] = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            out.append(d)
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out.append(n)
    return out


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def helib_bsgs_auts(D: int, minimal: bool = True) -> dict:
    # HElib uses BSGS for MatMul1D when D > 50, or in minimal mode when D > 8.
    use_bsgs = D > 50 or (minimal and D > 8)
    if not use_bsgs:
        return {"mode": "linear", "g": 0, "h": D, "automorphisms": max(D - 1, 0)}
    g = math.isqrt(D)
    if g * g < D:
        g += 1
    h = (D + g - 1) // g
    return {"mode": "bsgs", "g": g, "h": h, "automorphisms": max(g - 1, 0) + max(h - 1, 0)}


def radix_stage_cost(radices: list[int]) -> dict:
    # One radix-r stage can be implemented as r-1 masked rotations plus
    # plaintext twiddles when the stage is fused over all blocks.
    return {
        "radices": radices,
        "stages": len(radices),
        "automorphisms": sum(r - 1 for r in radices),
    }


def cooley_tukey_plan(D: int) -> dict:
    return radix_stage_cost(factor(D))


def primitive_prime_rader_plan(ell: int) -> dict | None:
    if not is_prime(ell) or ell <= 3:
        return None
    conv_len = ell - 1
    conv = cooley_tukey_plan(conv_len)
    # Two extra permutations/masked rotations are a conservative placeholder:
    # isolate the DC term and move between natural and generator order.
    return {
        "radices": conv["radices"],
        "stages": conv["stages"],
        "automorphisms": conv["automorphisms"] + 2,
        "convolution_length": conv_len,
    }


def format_plan(D: int, ell: int) -> str:
    base = helib_bsgs_auts(D)
    ct = cooley_tukey_plan(D)
    rd = primitive_prime_rader_plan(ell) if D == ell - 1 else None
    pieces = [
        f"D={D}",
        f"HElib {base['mode']} auts~{base['automorphisms']} g={base['g']} h={base['h']}",
        f"radix {ct['radices']} auts~{ct['automorphisms']}",
    ]
    if rd is not None:
        pieces.append(
            f"primitive-ell Rader conv={rd['convolution_length']} radix={rd['radices']} auts~{rd['automorphisms']}"
        )
    return " | ".join(pieces)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--set", choices=sorted(PARAMS), default="E")
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()

    labels = sorted(PARAMS) if args.all else [args.set]
    for label in labels:
        ps = PARAMS[label]
        print(f"set {label}: p={ps.p}, mvec={ps.mvec}, ords={ps.ords}")
        for ell, D in zip(ps.mvec, ps.ords):
            print("  " + format_plan(abs(D), ell))


if __name__ == "__main__":
    main()
