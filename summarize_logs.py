#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path


TIME_RE = re.compile(
    r"time for linear1 = ([0-9.]+), linear2 = ([0-9.]+), extract = ([0-9.]+), total = ([0-9.]+)"
)
QKS_RE = re.compile(r"final log2\(qks\) = ([0-9.]+)")
OK_RE = re.compile(r"### bts finished, everything ok ###")


def parse(path: Path) -> dict:
    text = path.read_text(errors="replace")
    times = TIME_RE.findall(text)
    if not times:
        raise SystemExit(f"no timing summary found in {path}")
    linear1, linear2, extract, total = map(float, times[-1])
    qks = QKS_RE.findall(text)
    return {
        "path": str(path),
        "qks": float(qks[-1]) if qks else None,
        "linear1": linear1,
        "linear2": linear2,
        "extract": extract,
        "total": total,
        "ok": bool(OK_RE.search(text)),
    }


def pct(base: float, opt: float) -> float:
    return (base - opt) / base * 100.0


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit("usage: summarize_logs.py LOG [LOG...]")
    rows = [parse(Path(arg)) for arg in sys.argv[1:]]
    print("| log | qks | linear1 | linear2 | extract | total | ok |")
    print("|---|---:|---:|---:|---:|---:|---|")
    for row in rows:
        qks = "-" if row["qks"] is None else f"{row['qks']:.6f}"
        print(
            f"| {Path(row['path']).name} | {qks} | {row['linear1']:.6f}s | "
            f"{row['linear2']:.6f}s | {row['extract']:.6f}s | {row['total']:.6f}s | {row['ok']} |"
        )
    if len(rows) >= 2:
        base, opt = rows[0], rows[1]
        print()
        print("| metric | speedup |")
        print("|---|---:|")
        for key in ("linear1", "linear2", "extract", "total"):
            print(f"| {key} | {pct(base[key], opt[key]):.2f}% |")


if __name__ == "__main__":
    main()
