#!/usr/bin/env python3
"""Report of 3_experiments/compare_baselines.sh: control groups vs treatment.

    python3 3_experiments/compare_report.py --dir out/compare
    python3 3_experiments/compare_report.py --dir 5_results/baseline_comparison   # shipped run

Reads <dir>/<experiment>/<group>.log and writes
    <dir>/summary.csv      one row per group: stage means, chain, radix, trials
    <dir>/report.md        the groups, the ratios that answer each question, PASS/FAIL
    <dir>/figures/*.png    stage breakdown per experiment (needs matplotlib)
"""
import argparse
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_tables as mt  # noqa: E402

AUX = re.compile(r"log2\(aux\) = ([\d.]+)")
ORDER4 = "HELIB_AUX_ORDER4_EVAL enabled"

GROUPS = {
    "m50731": [
        ("G0_original_ma_largep", "G0", "Ma et al. Large-p, original code (2_baselines/)"),
        ("G1_ours_off_radix35", "G1", "ours, evaluator OFF, baseline radix A = 35"),
        ("G2_ours_off_radix256", "G2", "ours, evaluator OFF, order-four radix A = 256"),
        ("G3_ours_on_radix256", "G3", "ours, evaluator ON, A = 256 (treatment)"),
        ("G4_ours_on_parallel", "G4", "G3 + parallel CoeffToSlot (2 cores)"),
    ],
    "po2": [
        ("G0_original_ma_ntt", "G0", "Ma et al. NTT, original code, own chain (1314 bits built)"),
        ("G1_ours_off_radix27", "G1", "ours, OFF, baseline radix A = 27, baseline chain (1314)"),
        ("G2_ours_off_radix256", "G2", "ours, OFF, order-four radix A = 256, baseline chain (1314)"),
        ("G3_ours_on_radix256", "G3", "ours, ON, A = 256, baseline chain (1314) (treatment)"),
        ("G4_ours_off_radix27_rec", "G4", "ours, OFF, A = 27, recommended 1066-bit chain (Table 5 OFF)"),
        ("G5_ours_off_radix256_rec", "G5", "ours, OFF, A = 256, recommended 1066-bit chain"),
        ("G6_ours_on_radix256_rec", "G6", "ours, ON, A = 256, recommended 1066-bit chain (Table 5 ON)"),
    ],
}

# (experiment, numerator, denominator, metric, question, expected (lo, hi) or None, paper value)
RATIOS = [
    ("m50731", "G0", "G1", "total", "Do our changes, switched off, leave the baseline unchanged?", (0.95, 1.05), None),
    ("m50731", "G1", "G2", "extract", "Effect of the radix alone (evaluator OFF)", None, None),
    ("m50731", "G2", "G3", "extract", "Effect of the evaluator alone (same radix)", (1.5, None), None),
    ("m50731", "G0", "G3", "total", "Ours vs original baseline, whole recryption (Table 8)", (1.35, None), 1.62),
    ("m50731", "G0", "G4", "total", "Combined vs original baseline (Table 8, 'combined')", (1.35, None), 1.80),
    ("po2", "G0", "G1", "total", "Do our changes, switched off, leave the baseline unchanged?", (0.95, 1.05), None),
    ("po2", "G1", "G2", "extract", "Effect of the radix alone, baseline chain (1314)", None, None),
    ("po2", "G2", "G3", "extract", "Effect of the evaluator alone, baseline chain (1314) (Table 5, last column, p = 65537)", (1.4, None), 1.59),
    ("po2", "G4", "G6", "extract", "Table 5 as measured: baseline radix OFF vs ON, 1066-bit chain", (1.5, None), 1.85),
    ("po2", "G5", "G6", "extract", "Effect of the evaluator alone, 1066-bit chain", (1.4, None), None),
    ("po2", "G4", "G6", "total", "Whole recryption, 1066-bit chain (paper: 42.34 s -> 33.25 s)", (1.1, None), 1.27),
]


def load(d):
    rows = {}
    for exp, groups in GROUPS.items():
        for fname, gid, desc in groups:
            path = d / exp / f"{fname}.log"
            r = mt.parse(path)
            rec = dict(experiment=exp, group=gid, description=desc, log=str(path.relative_to(d)))
            if r is None:
                rec["status"] = "not run"
            elif "error" in r:
                rec["status"] = "FAILED: " + r["error"]
            else:
                text = path.read_text(errors="replace")
                a = AUX.search(text)
                rec.update(status="ok", trials=r["n"], log2Q=r["logQ"],
                           radix=round(2 ** float(a.group(1))) if a else None,
                           order4=ORDER4 in text,
                           linear=r["linear1"] + r["linear2"], extract=r["extract"],
                           other=r["total"] - r["linear1"] - r["linear2"] - r["extract"],
                           total=r["total"])
            rows[(exp, gid)] = rec
    return rows


def write_csv(rows, path):
    cols = ["experiment", "group", "description", "status", "trials", "log2Q", "radix", "order4",
            "linear", "extract", "other", "total", "log"]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows.values():
            w.writerow({c: (f"{r[c]:.3f}" if isinstance(r.get(c), float) else r.get(c, ""))
                        for c in cols})


def figures(rows, outdir):
    try:
        import plotting as pl
    except ImportError:
        return []
    made = []
    for exp, groups in GROUPS.items():
        ok = [(gid, rows[(exp, gid)]) for _, gid, _ in groups if rows[(exp, gid)]["status"] == "ok"]
        if not ok:
            continue
        fig, ax = pl.plt.subplots(figsize=(8, 0.55 * len(ok) + 1.2))
        labels = [f"{gid}  {rows[(exp, gid)]['description']}" for gid, _ in ok]
        parts = [[r["linear"] for _, r in ok], [r["extract"] for _, r in ok],
                 [max(r["other"], 0) for _, r in ok]]
        pl.stacked_hbars(ax, labels, parts, ["linear transforms", "digit extraction", "other"])
        ax.set_xlabel("whole recryption (s), mean of trials")
        title = {"m50731": "Control groups, p = 65537, m = 50731 (vs Ma et al. Large-p)",
                 "po2": "Control groups, p = 65537, m = 2^16 (vs Ma et al. NTT)"}[exp]
        fig.suptitle(title, x=0.01, ha="left", fontsize=10)
        made.append(pl.save(fig, outdir / f"compare_{exp}"))
    return made


def report(rows, d):
    L = ["# Baseline comparison with control groups", "",
         "Generated by `3_experiments/compare_report.py` from the logs in this directory.",
         "Every run decrypts and checks all slots; times are means over the trials of one run.", ""]
    verdicts = []
    for exp, groups in GROUPS.items():
        L += [f"## {exp}", "", "| group | description | status | chain (bits) | radix A | evaluator | linear (s) | extract (s) | total (s) | trials |",
              "|---|---|---|---|---|---|---|---|---|---|"]
        for _, gid, desc in groups:
            r = rows[(exp, gid)]
            if r["status"] == "ok":
                L.append(f"| {gid} | {desc} | ok | {r['log2Q']:.0f} | {r['radix']} | "
                         f"{'ON' if r['order4'] else 'OFF'} | {r['linear']:.2f} | {r['extract']:.2f} | "
                         f"{r['total']:.2f} | {r['trials']} |")
            else:
                L.append(f"| {gid} | {desc} | {r['status']} | | | | | | | |")
        L += ["", "| question | ratio | value | expected | paper | verdict |", "|---|---|---|---|---|---|"]
        for e, num, den, metric, q, exp_rng, paper in RATIOS:
            if e != exp:
                continue
            a, b = rows[(e, num)], rows[(e, den)]
            if a["status"] != "ok" or b["status"] != "ok":
                L.append(f"| {q} | {num}/{den} ({metric}) | – | | {paper or ''} | not run |")
                continue
            v = a[metric] / b[metric]
            if exp_rng is None:
                verdict, rng = "informative", "–"
            else:
                lo, hi = exp_rng
                good = v >= lo and (hi is None or v <= hi)
                verdict = "PASS" if good else "FAIL"
                rng = f"{lo:.2f}–{hi:.2f}" if hi else f"≥ {lo:.2f}"
                verdicts.append(good)
            L.append(f"| {q} | {num}/{den} ({metric}) | **{v:.2f}x** | {rng} | "
                     f"{f'{paper:.2f}x' if paper else ''} | {verdict} |")
        L.append("")
    n_pass, n_all = sum(verdicts), len(verdicts)
    L += [f"**Checked ratios: {n_pass} of {n_all} pass.** "
          "Expected ranges are deliberately wide: a ratio of two timings on another machine "
          "moves by a few per cent, and more on a loaded one.", ""]
    (d / "report.md").write_text("\n".join(L) + "\n")
    return n_pass, n_all


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", type=Path, required=True)
    a = ap.parse_args()
    d = a.dir.resolve()
    rows = load(d)
    write_csv(rows, d / "summary.csv")
    n_pass, n_all = report(rows, d)
    made = figures(rows, d / "figures")
    print((d / "report.md").read_text())
    for m in made:
        print(f"wrote {m}")
    return 0 if n_pass == n_all and n_all > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
