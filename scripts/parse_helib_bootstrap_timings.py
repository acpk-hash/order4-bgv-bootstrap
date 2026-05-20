#!/usr/bin/env python3
"""Parse HElib fatboot timing lines into a compact JSON summary."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re


TIMING_RE = re.compile(
    r"time for linear1 = (?P<linear1>[0-9.]+), "
    r"linear2 = (?P<linear2>[0-9.]+), "
    r"extract = (?P<extract>[0-9.]+), "
    r"total = (?P<total>[0-9.]+)"
)


def mean(xs: list[float]) -> float:
    return sum(xs) / len(xs)


def sample_std(xs: list[float]) -> float:
    if len(xs) < 2:
        return 0.0
    mu = mean(xs)
    return math.sqrt(sum((x - mu) ** 2 for x in xs) / (len(xs) - 1))


def parse_log(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    rows = [
        {k: float(v) for k, v in match.groupdict().items()}
        for match in TIMING_RE.finditer(text)
    ]
    if not rows:
        raise ValueError(f"no timing lines found in {path}")

    passed = "### bts finished, everything ok ###" in text
    failed = "what():" in text or "thin results not match" in text
    runs = rows[:-1] if len(rows) > 1 else rows
    reported_average = rows[-1] if len(rows) > 1 else None

    stats = {}
    for key in ["linear1", "linear2", "extract", "total"]:
        values = [row[key] for row in runs]
        stats[key] = {
            "mean": mean(values),
            "sample_std": sample_std(values),
            "min": min(values),
            "max": max(values),
        }

    return {
        "path": str(path),
        "passed": passed,
        "failed": failed,
        "num_run_timings": len(runs),
        "runs": runs,
        "reported_average": reported_average,
        "stats": stats,
    }


def pct_improvement(base: float, opt: float) -> float:
    return (base - opt) / base * 100.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--optimized", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    baseline = parse_log(args.baseline)
    optimized = parse_log(args.optimized)
    improvements = {}
    for key in ["linear1", "linear2", "extract", "total"]:
        improvements[key] = pct_improvement(
            baseline["stats"][key]["mean"], optimized["stats"][key]["mean"]
        )

    result = {
        "baseline": baseline,
        "optimized": optimized,
        "improvement_percent_mean": improvements,
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    if args.json:
        print(text)
        return

    print("metric baseline_mean optimized_mean improvement_percent")
    for key in ["linear1", "linear2", "extract", "total"]:
        print(
            f"{key} {baseline['stats'][key]['mean']:.6f} "
            f"{optimized['stats'][key]['mean']:.6f} {improvements[key]:.2f}"
        )


if __name__ == "__main__":
    main()
