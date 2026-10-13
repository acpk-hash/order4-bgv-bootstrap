#!/usr/bin/env python3
"""Figures of the ciphertext experiments and of the security estimates.

    python3 3_experiments/make_figures.py                 # from 5_results/paper_runs
    python3 3_experiments/make_figures.py --logs out      # from your runs (missing rows are skipped)

Writes PNG and PDF to 5_results/figures/ (or --out):

    table4_po2.png        camera-ready Table 4 [5]: extraction time OFF vs ON, ten primes on m = 2^16
    table5_general.png    camera-ready Table 5 [6]: the same on nine general cyclotomic sets
    table6_m50731.png     camera-ready Table 6 [8]: whole recryption by stage, Ma et al. vs ours
    speedup_rerun.png     speedups of the paper next to our Docker reruns (5_results/reruns)
    security_lambda.png   camera-ready Table 3: lambda* at the built chain vs the 80/128-bit targets
"""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_tables as mt  # noqa: E402
import plotting as pl     # noqa: E402

ROOT = HERE.parent
PAPER = mt.PAPER


def fig_po2(logs, out):
    labels, off, on, sp = [], [], [], []
    for p, ref, a, b, s in mt.rows_po2(logs):
        if s is None:
            continue
        labels.append(p); off.append(a["extract"]); on.append(b["extract"]); sp.append(f"{s:.2f}x")
    if not labels:
        return None
    fig, ax = pl.plt.subplots(figsize=(8, 3.4))
    pl.grouped_bars(ax, labels, [off, on], ["OFF: Ma et al. evaluator, radix A = 27", "ON: X/2 + X^3 Q(X^4), order-four radix A"],
                    above=sp)
    ax.set_xlabel("plaintext prime p  (m = 2^16, 80-bit chain, h = 26)")
    ax.set_ylabel("digit extraction (s)")
    ax.set_title("Table 4: digit extraction on the power-of-two ring, mean of 3 trials")
    return pl.save(fig, out / "table4_po2")


def fig_general(logs, out):
    labels, off, on, sp = [], [], [], []
    for s_, ref, a, b, s in mt.rows_general(logs):
        if s is None:
            continue
        labels.append(s_); off.append(a["extract"]); on.append(b["extract"]); sp.append(f"{s:.2f}x")
    if not labels:
        return None
    fig, ax = pl.plt.subplots(figsize=(8, 3.4))
    pl.grouped_bars(ax, labels, [off, on], ["OFF: generic evaluation of P_A", "ON: X/2 + X^3 Q(X^4)"], above=sp)
    ax.set_xlabel("parameter set (80-x: h = 12, 128-x: h = 26)")
    ax.set_ylabel("digit extraction (s)")
    ax.set_title("Table 5: digit extraction on general cyclotomic rings (same radix A for OFF and ON), mean of 3 trials")
    return pl.save(fig, out / "table5_general")


def stages(r):
    lin = r["linear1"] + r["linear2"]
    return lin, r["extract"], max(r["total"] - lin - r["extract"], 0.0)


def fig_m50731(logs, out):
    rows = [(k, r) for k, ref, base, r, sp in mt.rows_m50731(logs) if mt.ok(r)]
    if not rows:
        return None
    names = {"ma": "Ma et al. (A = 35)", "order4": "ours, order-four (A = 256)",
             "combined": "ours + parallel CoeffToSlot"}
    parts = list(zip(*[stages(r) for _, r in rows]))
    fig, ax = pl.plt.subplots(figsize=(7.5, 2.6))
    pl.stacked_hbars(ax, [names[k] for k, _ in rows], parts,
                     ["linear transforms", "digit extraction", "other (key switching, mod-down)"])
    ax.set_xlabel("whole recryption (s)")
    fig.suptitle("Table 6: p = 65537, m = 50731, 2784 slots, by stage", x=0.01, ha="left", fontsize=10)
    return pl.save(fig, out / "table6_m50731")


def fig_rerun(out):
    rr = ROOT / "5_results" / "reruns"
    labels, paper, ours = [], [], []
    for p, ref, a, b, s in mt.rows_po2(rr):
        if s is not None:
            labels.append(f"T4 p={p}"); paper.append(ref["speedup"]); ours.append(s)
    for s_, ref, a, b, s in mt.rows_general(rr):
        if s is not None:
            labels.append(f"T5 {s_}"); paper.append(ref["speedup"]); ours.append(s)
    for k, ref, base, r, s in mt.rows_m50731(rr):
        if s is not None and k != "ma":
            labels.append(f"T6 {k}"); paper.append(ref["speedup"]); ours.append(s)
    if not labels:
        return None
    fig, ax = pl.plt.subplots(figsize=(7, 3.2))
    pl.grouped_bars(ax, labels, [paper, ours], ["paper", "Docker rerun (5_results/reruns)"],
                    value_fmt="{:.2f}x")
    ax.set_ylabel("speedup")
    ax.set_title("Speedups: paper vs our rerun of the artifact")
    return pl.save(fig, out / "speedup_rerun")


def fig_security(out):
    d = ROOT / "4_parameters" / "transcripts" / "built_chain"
    labels, vals, binding = [], [], []
    for s in PAPER["security_built_chain"]:
        if s.startswith("_"):
            continue
        f = d / f"{s}.json"
        if f.exists() and f.read_text().strip():
            j = json.loads(f.read_text())
            labels.append(s[:-1] + "-" + s[-1]); vals.append(j["lambda_star"])
    if not labels:
        return None
    fig, ax = pl.plt.subplots(figsize=(7, 3))
    bars = ax.bar(range(len(labels)), vals, 0.6, color=pl.SERIES[0], linewidth=0,
                  label="lambda* at the built chain (main key, h' = 120)")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.2f}", ha="center", va="bottom",
                fontsize=7, color=pl.INK_2)
    for i, t in enumerate((80, 128)):
        ax.axhline(t, color=pl.INK_2, linewidth=0.8, linestyle=(0, (4, 3)),
                   label="80- and 128-bit targets" if i == 0 else None)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels)
    ax.set_ylabel("bits"); ax.set_ylim(0, 175)
    ax.legend(loc="upper left", fontsize=8, ncol=2)
    ax.set_title("Table 3: security at the chain HElib built (Lattice Estimator 3e48ef4)")
    return pl.save(fig, out / "security_lambda")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--logs", type=Path, default=ROOT / "5_results" / "paper_runs")
    ap.add_argument("--out", type=Path, default=ROOT / "5_results" / "figures")
    a = ap.parse_args()
    logs = a.logs.resolve()
    made = [fig_po2(logs, a.out), fig_general(logs, a.out), fig_m50731(logs, a.out)]
    if logs == (ROOT / "5_results" / "paper_runs").resolve():
        made += [fig_rerun(a.out), fig_security(a.out)]
    for m in made:
        if m:
            print(f"wrote {m.relative_to(ROOT) if ROOT in m.parents else m} (+ .pdf)")


if __name__ == "__main__":
    main()
