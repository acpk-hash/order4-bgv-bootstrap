#!/usr/bin/env python3
"""Tabulate Tables 5, 6 and 7 of the paper from bootstrapping logs.

    python3 ae/make_tables.py                    # shipped logs (results/paper, results/)
    python3 ae/make_tables.py --logs ae/out      # your own runs from ae/run.sh

Every number is the mean of the per-trial lines
    time for linear1 = .., linear2 = .., extract = .., total = ..
that the drivers print, and a run counts only if it printed the
decryption check "### bts finished, everything ok ###". Each row is printed
next to the value in the paper; "dev" is the relative deviation of the
OFF/ON speedup, which is the quantity the paper claims. Absolute seconds
depend on the machine and should only be compared within one table.
"""
import argparse
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAPER = json.loads((ROOT / "ae" / "paper_values.json").read_text())
TRIAL = re.compile(r"time for linear1 = ([\d.]+), linear2 = ([\d.]+), "
                   r"extract = ([\d.]+), total = ([\d.]+)")
BITS = re.compile(r"INFO: total bits = ([\d.]+)")


def parse(path):
    """Return dict(linear1, linear2, extract, total, logQ, n) or None."""
    if not path.exists():
        return None
    text = path.read_text(errors="replace")
    if "everything ok" not in text:
        return {"error": "no decryption check"}
    trials = [tuple(map(float, m.groups())) for m in TRIAL.finditer(text)]
    # the drivers print one line per trial and then the average line
    if len(trials) > 1:
        trials = trials[:-1]
    n = len(trials)
    mean = [sum(t[i] for t in trials) / n for i in range(4)]
    b = BITS.search(text)
    return dict(linear1=mean[0], linear2=mean[1], extract=mean[2], total=mean[3],
                logQ=float(b.group(1)) if b else float("nan"), n=n)


def dev(ours, paper):
    return f"{100 * (ours - paper) / paper:+5.1f}%"


def cell(r, key):
    if r is None:
        return "   --  "
    if "error" in r:
        return " FAIL  "
    return f"{r[key]:7.2f}"


def tab5(logs):
    print("\nTable 5  m=2^16, digit extraction (s), 80-bit point")
    print(f"{'p':>6} {'logQ':>6} {'OFF':>7} {'ON':>7} {'speedup':>8} | "
          f"{'paper':>6} {'dev':>7}   n")
    for p, ref in PAPER["tab5_po2"].items():
        if p.startswith("_"):
            continue
        d = logs / "tab5_po2"
        off = parse(d / f"R_p{p}_OFF.out") or parse(d / f"V_p{p}_OFF.out")
        on = parse(d / f"R_p{p}_ON.out") or parse(d / f"V_p{p}_ON.out")
        line = f"{p:>6} "
        if off and on and "error" not in off and "error" not in on:
            s = off["extract"] / on["extract"]
            line += (f"{off['logQ']:6.0f} {off['extract']:7.2f} {on['extract']:7.2f} "
                     f"{s:7.2f}x | {ref['speedup']:5.2f}x {dev(s, ref['speedup'])}   "
                     f"{min(off['n'], on['n'])}")
        else:
            line += f"{'':6} {cell(off, 'extract')} {cell(on, 'extract')} {'':8} | {ref['speedup']:5.2f}x"
        print(line)


def tab6(logs):
    print("\nTable 6  general cyclotomic rings, digit extraction and whole recryption (s)")
    print(f"{'set':>5} {'logQ':>7} {'OFF':>7} {'ON':>7} {'speedup':>8} {'tot OFF':>8} "
          f"{'tot ON':>7} | {'paper':>6} {'dev':>7}   n")
    for s, ref in PAPER["tab6_general"].items():
        if s.startswith("_"):
            continue
        d = logs / "tab6_general"
        off, on = parse(d / f"{s}_off.log"), parse(d / f"{s}_on.log")
        if off and on and "error" not in off and "error" not in on:
            sp = off["extract"] / on["extract"]
            print(f"{s:>5} {off['logQ']:7.1f} {off['extract']:7.2f} {on['extract']:7.2f} "
                  f"{sp:7.2f}x {off['total']:8.2f} {on['total']:7.2f} | "
                  f"{ref['speedup']:5.2f}x {dev(sp, ref['speedup'])}   {min(off['n'], on['n'])}")
        else:
            print(f"{s:>5} {'':7} {cell(off, 'extract')} {cell(on, 'extract')} {'':8} "
                  f"{'':8} {'':7} | {ref['speedup']:5.2f}x")


def tab7(logs):
    print("\nTable 7  p=65537, m=50731 (s); speedup of total vs Ma et al.")
    d = logs / "tab7_m50731"
    rows = {k: parse(d / f"{f}.log") for k, f in
            (("ma", "ma_largep"), ("order4", "order4"), ("combined", "combined"))}
    base = rows["ma"]
    print(f"{'row':>9} {'linear':>7} {'extract':>8} {'total':>7} {'speedup':>8} | "
          f"{'paper':>6} {'dev':>7}   n")
    for k, r in rows.items():
        ref = PAPER["tab7_m50731"][k]
        if r and "error" not in r and base and "error" not in base:
            sp = base["total"] / r["total"]
            print(f"{k:>9} {r['linear1'] + r['linear2']:7.2f} {r['extract']:8.2f} "
                  f"{r['total']:7.2f} {sp:7.2f}x | {ref['speedup']:5.2f}x "
                  f"{dev(sp, ref['speedup'])}   {r['n']}")
        else:
            print(f"{k:>9} {'--':>7}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--logs", type=Path, default=ROOT / "results" / "paper",
                    help="directory holding tab5_po2/, tab6_general/, tab7_m50731/")
    ap.add_argument("tables", nargs="*", default=["5", "6", "7"])
    a = ap.parse_args()
    for t in a.tables:
        {"5": tab5, "6": tab6, "7": tab7}[t](a.logs.resolve())


if __name__ == "__main__":
    main()
