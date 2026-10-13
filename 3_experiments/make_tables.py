#!/usr/bin/env python3
"""Regenerate the ciphertext tables of the paper from bootstrapping logs.

    python3 3_experiments/make_tables.py                        # logs of the reported runs
    python3 3_experiments/make_tables.py --logs out          # your runs from 3_experiments/run.sh
    python3 3_experiments/make_tables.py --logs out po2      # one experiment only
    python3 3_experiments/make_tables.py --check shipped        # exit 1 unless every row matches the paper
    python3 3_experiments/make_tables.py --logs out --check rerun

Experiments (table numbers as in the camera-ready version and the full version):
    po2      Table 5   m = 2^16, ten primes, OFF vs ON
    general  Table 6   nine general cyclotomic sets, OFF vs ON
    m50731   Table 8   p = 65537, m = 50731: Ma et al. vs ours

A value is the mean over the per-trial lines

    time for linear1 = .., linear2 = .., extract = .., total = ..

printed by the drivers. A run counts only if it printed the decryption check
"### bts finished, everything ok ###". Every row shows the speedup next to the
paper's value and the relative deviation; the speedup (a ratio of two runs on
one machine) is the quantity the paper claims. Absolute seconds depend on the
machine and are never checked.

--check shipped  every present row within tolerance.shipped_logs_rel (0.5%)
--check rerun    every present row within tolerance.rerun_speedup_rel (15%) and
                 above the per-experiment floor; rows not run are skipped
"""
import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAPER = json.loads((ROOT / "3_experiments" / "paper_values.json").read_text())
TOL = PAPER["tolerance"]
TRIAL = re.compile(r"time for linear1 = ([\d.]+), linear2 = ([\d.]+), "
                   r"extract = ([\d.]+), total = ([\d.]+)")
BITS = re.compile(r"INFO: total bits = ([\d.]+)")
OK_MARK = "everything ok"


def parse(path):
    """Mean timings of one log, or None if absent, or {'error': ...}."""
    if not path.exists():
        return None
    text = path.read_text(errors="replace")
    if OK_MARK not in text:
        return {"error": "no decryption check ('everything ok' missing)"}
    trials = [tuple(map(float, m.groups())) for m in TRIAL.finditer(text)]
    if len(trials) > 1:      # the last line is the driver's own average
        trials = trials[:-1]
    if not trials:
        return {"error": "no timing lines"}
    n = len(trials)
    mean = [sum(t[i] for t in trials) / n for i in range(4)]
    b = BITS.search(text)
    return dict(linear1=mean[0], linear2=mean[1], extract=mean[2], total=mean[3],
                logQ=float(b.group(1)) if b else float("nan"), n=n)


def ok(r):
    return r is not None and "error" not in r


def where(key):
    w = PAPER[key]["_where"]
    return w['cr']


def rows_po2(logs):
    d = logs / "po2"
    for p, ref in PAPER["po2"].items():
        if p.startswith("_"):
            continue
        off, on = parse(d / f"R_p{p}_OFF.out"), parse(d / f"R_p{p}_ON.out")
        if off is None:
            off, on = parse(d / f"V_p{p}_OFF.out"), parse(d / f"V_p{p}_ON.out")
        yield p, ref, off, on, (off["extract"] / on["extract"] if ok(off) and ok(on) else None)


def rows_general(logs):
    d = logs / "general"
    for s, ref in PAPER["general"].items():
        if s.startswith("_"):
            continue
        off, on = parse(d / f"{s}_off.log"), parse(d / f"{s}_on.log")
        yield s, ref, off, on, (off["extract"] / on["extract"] if ok(off) and ok(on) else None)


def rows_m50731(logs):
    d = logs / "m50731"
    runs = {k: parse(d / f"{f}.log") for k, f in
            (("ma", "ma_largep"), ("order4", "order4"), ("combined", "combined"))}
    base = runs["ma"]
    for k, r in runs.items():
        sp = base["total"] / r["total"] if ok(base) and ok(r) else None
        yield k, PAPER["m50731"][k], base, r, sp


def dev(x, ref):
    return (x - ref) / ref


def status(r):
    if r is None:
        return "not run"
    return "FAIL: " + r["error"]


def show(name, logs):
    print()
    if name == "po2":
        print(f"po2: {where('po2')}\n  {PAPER['po2']['_caption']}")
        print(f"  {'p':>6} {'logQ':>6} {'OFF(s)':>8} {'ON(s)':>8} {'speedup':>8} {'paper':>7} {'dev':>7}  n")
        for p, ref, off, on, sp in rows_po2(logs):
            if sp is None:
                print(f"  {p:>6}  {status(off if not ok(off) else on)}")
                continue
            print(f"  {p:>6} {off['logQ']:6.0f} {off['extract']:8.2f} {on['extract']:8.2f} "
                  f"{sp:7.2f}x {ref['speedup']:6.2f}x {100*dev(sp, ref['speedup']):+6.1f}%  "
                  f"{min(off['n'], on['n'])}")
    elif name == "general":
        print(f"general: {where('general')}\n  {PAPER['general']['_caption']}")
        print(f"  {'set':>5} {'logQ':>7} {'OFF(s)':>8} {'ON(s)':>8} {'speedup':>8} "
              f"{'totOFF':>8} {'totON':>8} {'paper':>7} {'dev':>7}  n")
        for s, ref, off, on, sp in rows_general(logs):
            if sp is None:
                print(f"  {s:>5}  {status(off if not ok(off) else on)}")
                continue
            print(f"  {s:>5} {off['logQ']:7.1f} {off['extract']:8.2f} {on['extract']:8.2f} "
                  f"{sp:7.2f}x {off['total']:8.2f} {on['total']:8.2f} {ref['speedup']:6.2f}x "
                  f"{100*dev(sp, ref['speedup']):+6.1f}%  {min(off['n'], on['n'])}")
    elif name == "m50731":
        print(f"m50731: {where('m50731')}\n  {PAPER['m50731']['_caption']}")
        print(f"  {'row':>9} {'linear':>8} {'extract':>8} {'total':>8} {'speedup':>8} "
              f"{'paper':>7} {'dev':>7}  n")
        for k, ref, base, r, sp in rows_m50731(logs):
            if sp is None:
                print(f"  {k:>9}  {status(r if not ok(r) else base)}")
                continue
            print(f"  {k:>9} {r['linear1'] + r['linear2']:8.2f} {r['extract']:8.2f} "
                  f"{r['total']:8.2f} {sp:7.2f}x {ref['speedup']:6.2f}x "
                  f"{100*dev(sp, ref['speedup']):+6.1f}%  {r['n']}")


def check(name, logs, mode):
    """Return (passed, failed, skipped) lists of row labels."""
    gen = {"po2": rows_po2, "general": rows_general, "m50731": rows_m50731}[name]
    floor = {"po2": TOL["rerun_speedup_min_po2"],
             "general": TOL["rerun_speedup_min_general"],
             "m50731": TOL["rerun_speedup_min_m50731_order4"]}[name]
    passed, failed, skipped = [], [], []
    for label, ref, a, b, sp in gen(logs):
        tag = f"{name}/{label}"
        if name == "m50731" and label == "ma":
            if ok(b):
                passed.append(tag)      # baseline: only needs to decrypt correctly
            elif b is None:
                skipped.append(tag)
            else:
                failed.append(f"{tag} ({b['error']})")
            continue
        if sp is None:
            bad = [r for r in (a, b) if r is not None and "error" in r]
            (failed if bad else skipped).append(tag if not bad else f"{tag} ({bad[0]['error']})")
            continue
        d = abs(dev(sp, ref["speedup"]))
        if mode == "shipped":
            good = d <= TOL["shipped_logs_rel"]
        else:
            good = d <= TOL["rerun_speedup_rel"] and sp >= floor
        (passed if good else failed).append(
            tag if good else f"{tag} (speedup {sp:.2f}x vs paper {ref['speedup']:.2f}x)")
    return passed, failed, skipped


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--logs", type=Path, default=ROOT / "5_results" / "paper_runs",
                    help="directory holding po2/, general/, m50731/ (default: shipped logs)")
    ap.add_argument("--check", choices=["shipped", "rerun"],
                    help="exit with status 1 if a present row is outside tolerance")
    ap.add_argument("experiments", nargs="*", metavar="{po2,general,m50731}")
    a = ap.parse_args()
    a.experiments = a.experiments or ["po2", "general", "m50731"]
    for e in a.experiments:
        if e not in ("po2", "general", "m50731"):
            ap.error(f"unknown experiment {e!r}")
    logs = a.logs.resolve()
    for e in a.experiments:
        show(e, logs)
    if not a.check:
        return 0
    allp, allf, alls = [], [], []
    for e in a.experiments:
        p, f, s = check(e, logs, a.check)
        allp += p; allf += f; alls += s
    print(f"\ncheck ({a.check}): {len(allp)} passed, {len(allf)} failed, {len(alls)} not run")
    for f in allf:
        print(f"  FAIL {f}")
    return 1 if allf or not allp else 0


if __name__ == "__main__":
    sys.exit(main())
