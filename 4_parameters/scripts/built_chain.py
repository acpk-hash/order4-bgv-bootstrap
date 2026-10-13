#!/usr/bin/env sage -python
"""Per-attack security at the modulus chain HElib actually built (camera-ready
Table 3, column lambda*; full version Table 10, App. D).

Main-key instance: n = phi(m), q = 2^(integer log2 Q of the appendix table), sparse ternary
secret of weight h' = 120, sigma = 3.2, m = n samples, cost model MATZOV
(the estimator default at the pinned commit). Hybrid columns are the primal
hybrid at three switch settings, (mitm, babai) = FF, TF, TT.

    cd 4_parameters/scripts       # uses ../lattice-estimator-3e48ef4 (pinned copy)
    sage -python built_chain.py 128A > ../transcripts/built_chain/128A.json

Each set takes 10-40 min single-threaded.
"""
import json
import math
import os
import sys
import time

sys.path.insert(0, os.environ.get("LATTICE_ESTIMATOR",
                                  os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                                               "lattice-estimator-3e48ef4")))
from estimator import LWE, ND  # noqa: E402
from estimator.io import Logging  # noqa: E402

Logging.set_level(Logging.LEVEL0)

# set: (m, phi(m), built log2 Q, integer log2 q used) -- the chains printed by
# the general-ring runs (camera-ready Table 5), and the integer modulus of the appendix table (128D: the
# displayed 1267.5 was rounded up to 1268, the conservative direction).
SETS = {
    "80A": (37637, 36960, 1524.2, 1524), "80B": (43033, 41140, 1642.05, 1642),
    "80C": (45193, 44100, 1786.76, 1787), "80D": (50731, 50112, 2026.15, 2026),
    "128A": (65047, 62776, 1406.93, 1407), "128B": (53983, 52488, 1153.03, 1153),
    "128C": (55459, 53640, 1193.06, 1193), "128D": (57521, 56832, 1267.48, 1268),
    "128E": (54359, 51480, 1132.56, 1133),
}
H_MAIN, SIGMA = 120, 3.2


def bits(cost):
    v = float(cost["rop"])
    return math.log2(v) if 0 < v < float("inf") else None


def main(name):
    m, n, logq, logq_int = SETS[name]
    params = LWE.Parameters(n=n, q=2 ** logq_int,
                            Xs=ND.SparseTernary(H_MAIN // 2, H_MAIN - H_MAIN // 2, n),
                            Xe=ND.DiscreteGaussian(SIGMA), m=n, tag=name)
    attacks = {
        "usvp": lambda: LWE.primal_usvp(params),
        "bdd": lambda: LWE.primal_bdd(params),
        "dual": lambda: LWE.dual(params),
        "hyb_FF": lambda: LWE.primal_hybrid(params, mitm=False, babai=False),
        "hyb_TF": lambda: LWE.primal_hybrid(params, mitm=True, babai=False),
        "hyb_TT": lambda: LWE.primal_hybrid(params, mitm=True, babai=True),
    }
    out = dict(set=name, m=m, n=n, built_log2Q=logq, log2q_used=logq_int,
               h_main=H_MAIN, sigma=SIGMA, estimator_commit="3e48ef4")
    t0 = time.time()
    for k, f in attacks.items():
        try:
            out[k] = bits(f())
        except Exception as e:  # record, do not hide
            out[k] = None
            out[k + "_error"] = repr(e)
    vals = {k: out[k] for k in attacks if out[k] is not None}
    out["lambda_star"] = min(vals.values())
    out["binding"] = min(vals, key=vals.get)
    out["seconds"] = round(time.time() - t0, 1)
    print(json.dumps(out))


if __name__ == "__main__":
    main(sys.argv[1])
