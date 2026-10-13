#!/usr/bin/env python3
"""Security at the built modulus chain (camera-ready Table 3, column lambda*;
full version Table 10, App. D).

    python3 3_experiments/tests/test_security.py                  # check the stored transcripts
    python3 3_experiments/tests/test_security.py --rerun 80B      # recompute one set first (SageMath,
                                                       # Lattice Estimator 3e48ef4, ~20 min)

Every per-attack cost in 4_parameters/transcripts/built_chain/<set>.json must equal the
paper to within tolerance.security_abs_bits (0.02 bit), and lambda* must be
the minimum over the attacks.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = json.loads((ROOT / "3_experiments" / "paper_values.json").read_text())
REF = PAPER["security_built_chain"]
TOL = PAPER["tolerance"]["security_abs_bits"]
DIR = ROOT / "4_parameters" / "transcripts" / "built_chain"
COLS = ["usvp", "bdd", "dual", "hyb_FF", "hyb_TF", "hyb_TT", "lambda_star"]


def rerun(name, outdir):
    est = os.environ.get("LATTICE_ESTIMATOR", str(ROOT / "4_parameters" / "lattice-estimator-3e48ef4"))
    if not Path(est, "estimator").is_dir():
        sys.exit(f"Lattice Estimator not found at {est}")
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"  recomputing {name} (about 20 min) ...", flush=True)
    res = subprocess.run(["sage", "-python", str(ROOT / "4_parameters" / "scripts" / "built_chain.py"), name],
                         capture_output=True, text=True, env={**os.environ, "LATTICE_ESTIMATOR": est})
    if res.returncode:
        sys.exit(res.stderr)
    (outdir / f"{name}.json").write_text(res.stdout)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rerun", nargs="*", metavar="SET", help="recompute these sets with SageMath")
    ap.add_argument("--dir", type=Path, default=None,
                    help="where to read results (default: stored transcripts, or out/security after --rerun)")
    a = ap.parse_args()
    d = a.dir or (ROOT / "out" / "security" if a.rerun else DIR)
    for s in a.rerun or []:
        rerun(s, d)

    fails, checked = [], 0
    print(f"Security at the built chain, main key h' = 120 ({d.relative_to(ROOT)})")
    print(f"  {'set':>5} " + " ".join(f"{c:>11}" for c in COLS))
    for s, ref in REF.items():
        if s.startswith("_"):
            continue
        f = d / f"{s}.json"
        if not f.exists() or not f.read_text().strip():
            continue
        got = json.loads(f.read_text())
        cells = []
        for c in COLS:
            v = got.get(c)
            if c in ref:
                checked += 1
                good = v is not None and abs(v - ref[c]) <= TOL + 0.005
                if not good:
                    fails.append(f"{s}.{c}: {v} vs paper {ref[c]}")
                cells.append(f"{v:8.2f}{'  ' if good else ' X'} ")
            else:
                cells.append(f"{v:8.2f}   " if v is not None else f"{'-':>11}")
        attacks = [got[c] for c in COLS[:-1] if got.get(c) is not None]
        if abs(min(attacks) - got["lambda_star"]) > 1e-9:
            fails.append(f"{s}: lambda* is not the minimum over the attacks")
        print(f"  {s:>5} " + " ".join(cells))
    print(f"\n  {checked} values compared with the paper (tolerance {TOL} bit; X = outside)")
    if fails or not checked:
        print("FAILED" if fails else "NOTHING CHECKED")
        for f in fails:
            print("  " + f)
        return 1
    print("PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
