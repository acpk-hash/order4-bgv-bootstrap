#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Security estimation under the rotated primal hybrid of ePrint 2026/279.

  sage -python rot_estimate.py --sets sets.json --out results/rot.jsonl
  sage -python rot_estimate.py --sweep-q 32768 1150 1265 5 --h 120
  sage -python rot_estimate.py --sweep-h 32768 1334 250 440 10

Requires the attack author's repository:

    git clone https://github.com/TabOg/mlwe-hybrids
    cd mlwe-hybrids/PrimalHybrid
    git -C lattice_estimator checkout 3e48ef4     # same commit as the stock tool

and expects it at --primal-dir, default ~/mlwe-hybrids/PrimalHybrid.

Two points of method, both of which are easy to get wrong.

First, the baseline and the attack must come from the same optimiser. In that
code `rot_search_space = search_space * poly_degree` and
`rot_hit_probability = 1 - (1 - hit_probability) ** poly_degree`, so
poly_degree = 1 is exactly the plain hybrid and poly_degree = n is the rotated
one. Comparing the rotated cost against a differently-optimised baseline mixes
two effects.

Second, the post-attack security of an instance is the minimum over every
attack considered, not the stock figure minus a delta. The stock estimator's
minimum frequently comes from an attack cheaper than this tool's own plain
hybrid, and subtracting a delta from it double counts. On two of our instances
the rotated cost is higher than the plain one, so a delta would even be
negative.
"""
from __future__ import print_function
import argparse, json, math, os, sys, time


def load_estimator(primal_dir):
    primal_dir = os.path.expanduser(primal_dir)
    est_dir = os.path.join(primal_dir, "lattice_estimator")
    sys.path.insert(0, est_dir)
    sys.path.insert(0, primal_dir)
    os.chdir(primal_dir)
    from estimator import LWE, ND, RC          # noqa
    from lwe_rot_primal import rot_primal_hybrid  # noqa
    from sage.all import oo                    # noqa
    return LWE, ND, RC, rot_primal_hybrid, oo


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primal-dir", default="~/mlwe-hybrids/PrimalHybrid")
    ap.add_argument("--sets", help="JSON list of {tag,n,log2q,h}")
    ap.add_argument("--out", default="-")
    ap.add_argument("--sigma", type=float, default=3.2)
    ap.add_argument("--sweep-q", nargs=4, type=int, metavar=("N", "LO", "HI", "STEP"))
    ap.add_argument("--sweep-h", nargs=4, type=int, metavar=("N", "LOG2Q", "LO", "HI"))
    ap.add_argument("--h", type=int, default=120, help="weight for --sweep-q")
    ap.add_argument("--step", type=int, default=10, help="step for --sweep-h")
    ap.add_argument("--heuristics", default="square root,estimator")
    a = ap.parse_args()

    LWE, ND, RC, rot_primal_hybrid, oo = load_estimator(a.primal_dir)
    COST = RC.MATZOV

    def lg(res):
        try:
            v = res["rop"]
            return None if v == oo else float(v.log2())
        except Exception:
            return None

    tasks = []
    if a.sets:
        for s in json.load(open(a.sets)):
            tasks.append((s["tag"], s["n"], s["log2q"], s["h"]))
    if a.sweep_q:
        n, lo, hi, st = a.sweep_q
        for q in range(lo, hi + 1, st):
            tasks.append(("q%d-n%d" % (q, n), n, q, a.h))
    if a.sweep_h:
        n, q, lo, hi = a.sweep_h
        for h in range(lo, hi + 1, a.step):
            tasks.append(("h%d-q%d" % (h, q), n, q, h))
    if not tasks:
        ap.error("give --sets, --sweep-q or --sweep-h")

    out = sys.stdout if a.out == "-" else open(a.out, "a")
    heurs = [x.strip() for x in a.heuristics.split(",") if x.strip()]

    for tag, n, log2q, h in tasks:
        hp = (h + 1) // 2                    # ceil(h/2) positive entries
        params = LWE.Parameters(
            n=n, q=2 ** log2q,
            Xs=ND.SparseTernary(p=hp, m=h - hp, n=n),
            Xe=ND.DiscreteGaussian(stddev=a.sigma))
        rec = dict(tag=tag, n=n, log2q=log2q, h=h, sigma=a.sigma,
                   cost_model="RC.MATZOV")
        for heur in heurs:
            t0 = time.time()
            try:
                # poly_degree = n is the rotated attack, poly_degree = 1 the
                # plain hybrid from the very same optimiser
                r = rot_primal_hybrid(params, babai=True, mitm=True,
                                      poly_degree=n, mitm_heuristic=heur,
                                      red_cost_model=COST)
                key = "rot_" + heur.replace(" ", "_")
                rec[key] = lg(r)
                for extra in ("zeta", "beta"):
                    if extra in r:
                        rec[key + "_" + extra] = int(r[extra])
            except Exception as e:
                rec["err_" + heur.replace(" ", "_")] = "%s: %s" % (type(e).__name__, e)
            rec["sec_" + heur.replace(" ", "_")] = round(time.time() - t0, 1)
        vals = [v for k, v in rec.items()
                if k.startswith("rot_") and not k.endswith(("_zeta", "_beta"))
                and isinstance(v, float)]
        rec["rot_min"] = min(vals) if vals else None
        rec["ceiling_log2n"] = round(math.log(n, 2), 2)
        out.write(json.dumps(rec, sort_keys=True) + "\n")
        out.flush()
        print("%-16s n=%-6d log2q=%-5d h=%-4d rot_min=%s"
              % (tag, n, log2q, h,
                 ("%.4f" % rec["rot_min"]) if rec["rot_min"] else "FAIL"),
              file=sys.stderr)


if __name__ == "__main__":
    main()
