#!/usr/bin/env python
# Security estimation for BGV bootstrapping parameter sets with sparse-secret
# key encapsulation.  Two LWE instances per set (MHWW24 Sec 4.3 model):
#   (i)  main key      : (n = phi(m), Q = full modulus chain, h' = 120)
#   (ii) encapsulated  : (n = phi(m), q_boot = q0*R,          h  = 7/12/26)
# Scheme security = min(i, ii).
#
# We record WHICH attack attains the minimum, because both eprint 2026/279
# (Remark 4) and eprint 2026/366 (Remarks 1.1, 3.7) state that their ring
# speedup applies only to the guessing step of a HYBRID attack: when the
# best attack is plain uSVP/BDD/dual (guessing dimension k = 0) the ring
# structure "cannot be leveraged to improve the probability".
import sys, os, time, math, json, argparse

EST_PATH = os.path.expanduser("~/sagesetup/lattice-estimator")
sys.path.insert(0, EST_PATH)

from estimator import *  # noqa
from estimator.io import Logging
Logging.set_level(Logging.LEVEL0)

SIGMA = 3.2  # HElib default stdev

HYBRID_ATTACKS = {"primal_hybrid", "primal_hybrid_nb", "dual_hybrid"}


def lam_of(rop):
    try:
        v = float(rop)
    except Exception:
        return None
    if v <= 0 or math.isinf(v):
        return None
    return math.log2(v)


def build(n, log2q, h, m_samples, label):
    q = 2 ** int(round(log2q))
    p1 = h // 2
    m1 = h - p1
    Xs = ND.SparseTernary(p1, m1, n)
    Xe = ND.DiscreteGaussian(SIGMA)
    kw = dict(n=n, q=q, Xs=Xs, Xe=Xe, tag=label)
    if m_samples is not None:
        kw["m"] = m_samples
    return LWE.Parameters(**kw)


def estimate(n, log2q, h, m_samples, label, quick=False):
    params = build(n, log2q, h, m_samples, label)
    out = {}
    t0 = time.time()
    # NOTE: LWE.primal_hybrid does not optimise over the `babai` switch, and on
    # very sparse instances babai=False can be several bits cheaper, so we run
    # both and keep the minimum.
    fns = [("primal_usvp", LWE.primal_usvp), ("dual", LWE.dual)]
    if not quick:
        fns += [
            ("primal_bdd", LWE.primal_bdd),
            ("primal_hybrid", LWE.primal_hybrid),
            ("primal_hybrid_nb", lambda pp: LWE.primal_hybrid(pp, babai=False)),
            ("dual_hybrid", LWE.dual_hybrid),
        ]
    for name, fn in fns:
        try:
            r = fn(params)
            out[name] = lam_of(r["rop"])
            # record guessing dimension zeta / eta where available
            for extra in ("zeta", "eta", "beta"):
                if extra in r:
                    try:
                        out["%s_%s" % (name, extra)] = int(r[extra])
                    except Exception:
                        pass
        except Exception as e:
            out["ERR_" + name] = str(e)[:100]
    cand = {k: v for k, v in out.items()
            if v is not None and not k.startswith("ERR_") and "_" not in k.replace("primal_", "P").replace("dual_", "D")}
    # simpler: explicit attack keys
    attack_keys = ["primal_usvp", "primal_bdd", "primal_hybrid", "primal_hybrid_nb",
                   "dual", "dual_hybrid"]
    cand = {k: out[k] for k in attack_keys if out.get(k) is not None}
    if cand:
        best = min(cand, key=cand.get)
        out["min"] = cand[best]
        out["argmin"] = best
        out["hybrid_wins"] = best in HYBRID_ATTACKS
        # how much cheaper is the hybrid than the best non-hybrid?
        nonhyb = {k: v for k, v in cand.items() if k not in HYBRID_ATTACKS}
        hyb = {k: v for k, v in cand.items() if k in HYBRID_ATTACKS}
        if nonhyb and hyb:
            out["hybrid_margin"] = round(min(nonhyb.values()) - min(hyb.values()), 2)
    else:
        out["min"] = None
        out["argmin"] = None
    out["secs"] = round(time.time() - t0, 1)
    return out


def show(tag, n, log2q, h, res):
    parts = []
    for k in ["primal_usvp", "primal_bdd", "primal_hybrid", "dual", "dual_hybrid"]:
        if res.get(k) is not None:
            parts.append("%s=%.1f" % (k, res[k]))
    print(
        "%-30s n=%-6d log2q=%-8.2f h=%-4d -> MIN=%s via %-13s | %s | [%.0fs]"
        % (tag, n, log2q, h,
           ("%.1f" % res["min"]) if res["min"] else "FAIL",
           res.get("argmin") or "-",
           ", ".join(parts), res["secs"]),
        flush=True,
    )


# MHWW24 (eprint 2024/115) Table 4, Set V / V-A: p^r=65537, M=50731, d=18,
#   log2(Q)=2036, lambda'=82.3 (main, h'=120)
#   encapsulated: h=24 -> 136.8 ; h=12 -> 82.3 ; log2(q0 R)=81.2
VALIDATION = [
    ("MHWW24 V-A main  exp 82.3", 50112, 2036.0, 120),
    ("MHWW24 V-A enc h24 exp 136.8", 50112, 81.2, 24),
    ("MHWW24 V-A enc h12 exp 82.3", 50112, 81.2, 12),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--sets", default=None)
    ap.add_argument("--msamples", default="n", choices=["n", "inf", "2n"])
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--only", default=None, help="comma-separated family filter")
    ap.add_argument("--slice", default=None, help="a:b index slice applied after --only")
    ap.add_argument("--encap-only", action="store_true",
                    help="skip the main-key instance (already computed elsewhere)")
    ap.add_argument("--main-only", action="store_true",
                    help="skip the encapsulated instance")
    a = ap.parse_args()

    def resolve(n):
        if a.msamples == "inf":
            return None
        if a.msamples == "n":
            return n
        return 2 * n

    if a.validate:
        print("=" * 118)
        print("STAGE 1 validation, m_samples=%s  (reproduce MHWW24 Table 4 Set V-A)" % a.msamples)
        print("=" * 118, flush=True)
        for tag, n, lq, h in VALIDATION:
            res = estimate(n, lq, h, resolve(n), tag, a.quick)
            show(tag, n, lq, h, res)
        return

    sets = json.load(open(a.sets))
    if a.only:
        fams = set(a.only.split(","))
        sets = [s for s in sets if s["family"] in fams]
    if a.slice:
        lo, hi = a.slice.split(":")
        sets = sets[int(lo) if lo else None:int(hi) if hi else None]
    rows = []
    print("=" * 118, flush=True)
    for s in sets:
        n = s["phim"]
        if a.encap_only:
            rm = {"min": None, "argmin": None, "secs": 0.0}
        else:
            rm = estimate(n, s["log2Q"], s["hmain"], resolve(n), s["name"] + "/main", a.quick)
            show(s["name"] + " MAIN", n, s["log2Q"], s["hmain"], rm)
        if a.main_only:
            re_ = {"min": None, "argmin": None, "secs": 0.0}
        else:
            re_ = estimate(n, s["log2qboot"], s["h"], resolve(n), s["name"] + "/encap", a.quick)
            show(s["name"] + " ENCAP", n, s["log2qboot"], s["h"], re_)
        row = dict(s)
        row["main"] = rm
        row["encap"] = re_
        vals = [x for x in [rm.get("min"), re_.get("min")] if x is not None]
        row["lam_min"] = min(vals) if vals else None
        row["binding"] = "main" if (rm.get("min") is not None and row["lam_min"] == rm.get("min")) else "encap"
        rows.append(row)
        if a.out:
            json.dump(rows, open(a.out, "w"), indent=1)
    print("=" * 118)
    print("%-24s %-7s %8s %8s %8s  %-9s %-14s %-14s" %
          ("set", "phi(m)", "main", "encap", "MIN", "binding", "main_argmin", "encap_argmin"))
    for r in rows:
        print("%-24s %-7d %8s %8s %8s  %-9s %-14s %-14s" % (
            r["name"], r["phim"],
            "%.1f" % r["main"]["min"] if r["main"]["min"] else "-",
            "%.1f" % r["encap"]["min"] if r["encap"]["min"] else "-",
            "%.1f" % r["lam_min"] if r["lam_min"] else "-",
            r["binding"], r["main"].get("argmin") or "-", r["encap"].get("argmin") or "-"), flush=True)


if __name__ == "__main__":
    main()
