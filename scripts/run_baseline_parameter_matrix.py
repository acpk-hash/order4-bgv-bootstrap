#!/usr/bin/env python3
"""Run same-parameter cleaner and ciphertext comparisons against baselines.

The script is intentionally conservative.  It does not modify the baseline
repositories.  For each parameter set it:

1. runs the baseline fatboot variants that are available locally;
2. parses the actual overflow bound B from the HElib log;
3. searches exact bounded-support cleaners for the same (p, B);
4. runs the isolated optimized artifact with the selected aux radix when the
   non-power-of-p path is valid.

The output is a JSON summary plus a Markdown report under
results/baseline_matrix/.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.append(str(SCRIPT_DIR))

from character_projected_cleaner_search import (  # noqa: E402
    candidate_decompositions,
    default_candidates,
    interpolate_cleaner,
    strip_terms,
)
from parse_helib_bootstrap_timings import parse_log  # noqa: E402


BOUND_RE = re.compile(r"bound on I = (?P<bound>[0-9.]+)")
DEG_RE = re.compile(r"deg of poly for non-power-of-p aux is (?P<deg>\d+)")
AUX_RE = re.compile(r"log2\(aux\) = (?P<log2aux>[0-9.]+)")
EXPLICIT_AUX_RE = re.compile(r"using HELIB_EXPLICIT_AUX=(?P<aux>\d+)")


MA_PARAMS = [
    {"article": "Ma-large-p", "label": "A", "i": 0, "p": 17, "h": 14},
    {"article": "Ma-large-p", "label": "B", "i": 1, "p": 127, "h": 12},
    {"article": "Ma-large-p", "label": "C", "i": 2, "p": 257, "h": 12},
    {"article": "Ma-large-p", "label": "D", "i": 3, "p": 8191, "h": 12},
    {"article": "Ma-large-p", "label": "E", "i": 4, "p": 65537, "h": 12},
]

NTT_PARAMS = [
    {"article": "Ma-NTT", "label": "I", "i": 0, "p": 65537, "h": 26},
    {"article": "Ma-NTT", "label": "II", "i": 1, "p": 8191, "h": 24},
    {"article": "Ma-NTT", "label": "III", "i": 2, "p": 131071, "h": 26},
]


def run_command(
    cmd: list[str],
    log_path: Path,
    timeout_s: int,
    env: dict[str, str] | None = None,
    force: bool = False,
) -> dict[str, Any]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if log_path.exists() and not force:
        return {"status": "cached", "log": str(log_path), "returncode": 0}
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    try:
        with log_path.open("w", encoding="utf-8", errors="replace") as fh:
            proc = subprocess.run(
                cmd,
                cwd=str(ROOT),
                env=merged_env,
                stdout=fh,
                stderr=subprocess.STDOUT,
                timeout=timeout_s,
                check=False,
            )
        return {
            "status": "ok" if proc.returncode == 0 else "failed",
            "log": str(log_path),
            "returncode": proc.returncode,
        }
    except subprocess.TimeoutExpired:
        with log_path.open("a", encoding="utf-8", errors="replace") as fh:
            fh.write(f"\n### TIMEOUT after {timeout_s}s ###\n")
        return {"status": "timeout", "log": str(log_path), "returncode": None}


def parse_extra_log_fields(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    text = path.read_text(encoding="utf-8", errors="replace")
    out: dict[str, Any] = {
        "passed_marker": "### bts finished, everything ok ###" in text,
        "timeout_marker": "### TIMEOUT" in text,
    }
    if m := BOUND_RE.search(text):
        bound = float(m.group("bound"))
        out["bound_I"] = bound
        out["B"] = int(math.ceil(bound))
    if m := DEG_RE.search(text):
        out["extract_poly_degree_from_log"] = int(m.group("deg"))
    if m := AUX_RE.search(text):
        out["log2_aux"] = float(m.group("log2aux"))
    if m := EXPLICIT_AUX_RE.search(text):
        out["explicit_aux"] = int(m.group("aux"))
    try:
        out["timing"] = parse_log(path)
    except Exception as exc:  # noqa: BLE001 - report parser failure in artifact
        out["timing_parse_error"] = str(exc)
    return out


def support_type_b_feasible(p: int, B: int) -> bool:
    return (2 * B + 1) * (2 * B + 1) < p


def search_cleaners(p: int, B: int, top: int = 8) -> dict[str, Any]:
    candidates = set(default_candidates(p, B, max_order=64))
    candidates.add(2 * B + 1)
    candidates.add((-(2 * B + 1)) % p)
    rows = []
    for aux in sorted(a % p for a in candidates if a % p not in {0, 1}):
        row = interpolate_cleaner(p, B, aux)
        if row.get("valid"):
            decomps = candidate_decompositions(row, 16)
            row["best_decomposition"] = decomps[0] if decomps else None
            row["decompositions"] = decomps[:5]
            row["score"] = (
                row["best_decomposition"]["estimated_mults"]
                if row["best_decomposition"]
                else row["generic_ps"]["mults"],
                row["nonzero_terms"],
                abs(row["centered_aux"]),
            )
        rows.append(row)
    valid = [r for r in rows if r.get("valid")]
    invalid = [r for r in rows if not r.get("valid")]
    valid.sort(key=lambda r: r["score"])
    generic_aux = 2 * B + 1
    generic = strip_terms(interpolate_cleaner(p, B, generic_aux % p))
    best = [strip_terms(r) for r in valid[:top]]
    chosen = None
    for row in valid:
        centered = row.get("centered_aux", row["aux"])
        if centered > 1:
            chosen = strip_terms(row)
            break
    if chosen is None and valid:
        chosen = strip_terms(valid[0])
    return {
        "p": p,
        "B": B,
        "support_size": (2 * B + 1) ** 2,
        "num_candidates": len(rows),
        "num_valid": len(valid),
        "num_invalid": len(invalid),
        "generic_aux": generic_aux,
        "generic_aux_result": generic,
        "best": best,
        "chosen_practical": chosen,
        "type_b_feasible": support_type_b_feasible(p, B),
    }


def mean_metric(row: dict[str, Any], key: str) -> float | None:
    timing = row.get("timing")
    if not timing:
        return None
    try:
        return float(timing["stats"][key]["mean"])
    except Exception:
        return None


def fmt(x: Any, digits: int = 3) -> str:
    if x is None:
        return "--"
    if isinstance(x, float):
        return f"{x:.{digits}f}"
    return str(x)


def pct(base: float | None, opt: float | None) -> float | None:
    if base is None or opt is None or base == 0:
        return None
    return (base - opt) / base * 100.0


def make_markdown(summary: dict[str, Any]) -> str:
    lines: list[str] = []
    lines.append("# Baseline Parameter Matrix Report")
    lines.append("")
    lines.append("This report is generated from local logs.  Timings are real ciphertext")
    lines.append("fatboot runs when the corresponding status is `ok` or `cached` with a")
    lines.append("passing correctness marker.")
    lines.append("")

    for family in ["Ma-large-p", "Ma-NTT"]:
        rows = [r for r in summary["parameters"] if r["article"] == family]
        if not rows:
            continue
        lines.append(f"## {family}")
        lines.append("")
        lines.append(
            "| set | p | h | B | type-B? | generic poly | best/chosen poly | variants |"
        )
        lines.append("|---|---:|---:|---:|---|---|---|---|")
        for r in rows:
            search = r.get("cleaner_search") or {}
            generic = search.get("generic_aux_result") or {}
            chosen = search.get("chosen_practical") or {}
            best_decomp = (chosen.get("best_decomposition") or {})
            generic_s = (
                f"A={search.get('generic_aux')}, deg={generic.get('degree')}, "
                f"terms={generic.get('nonzero_terms')}"
            )
            chosen_s = (
                f"A={chosen.get('aux')}, deg={chosen.get('degree')}, "
                f"terms={chosen.get('nonzero_terms')}, "
                f"d={best_decomp.get('modulus')}, est={best_decomp.get('estimated_mults')}"
            )
            variant_s = ", ".join(
                f"{name}:{v.get('status')}" for name, v in r.get("variants", {}).items()
            )
            lines.append(
                f"| {r['label']} | {r['p']} | {r['h']} | {fmt(r.get('B'), 0)} | "
                f"{search.get('type_b_feasible')} | {generic_s} | {chosen_s} | {variant_s} |"
            )
        lines.append("")

        lines.append("### Ciphertext Timings")
        lines.append("")
        lines.append(
            "| set | variant | status | pass | linear1 | linear2 | extract | total | "
            "extract speedup vs Ma-typeB | total speedup vs Ma-typeB |"
        )
        lines.append("|---|---|---|---|---:|---:|---:|---:|---:|---:|")
        for r in rows:
            base = r.get("variants", {}).get("ma_typeB")
            base_extract = mean_metric(base or {}, "extract")
            base_total = mean_metric(base or {}, "total")
            for name, v in r.get("variants", {}).items():
                extract = mean_metric(v, "extract")
                total = mean_metric(v, "total")
                lines.append(
                    f"| {r['label']} | {name} | {v.get('status')} | "
                    f"{v.get('passed_marker')} | {fmt(mean_metric(v, 'linear1'))} | "
                    f"{fmt(mean_metric(v, 'linear2'))} | {fmt(extract)} | {fmt(total)} | "
                    f"{fmt(pct(base_extract, extract), 2)} | {fmt(pct(base_total, total), 2)} |"
                )
        lines.append("")

    lines.append("## Non-Comparable Local Baselines")
    lines.append("")
    for item in summary.get("non_comparable", []):
        lines.append(f"- `{item['directory']}`: {item['reason']}")
    lines.append("")
    return "\n".join(lines)


def copy_existing(src: Path, dst: Path, force: bool) -> dict[str, Any] | None:
    if not src.exists():
        return None
    dst.parent.mkdir(parents=True, exist_ok=True)
    if force or not dst.exists():
        shutil.copy2(src, dst)
    return {"status": "cached", "log": str(dst), "returncode": 0}


def run_ma_large_p(args: argparse.Namespace, out_dir: Path) -> list[dict[str, Any]]:
    base_bin = Path("/home/user/experiments/baselines/BGV-Boot-for-Large-p/build/fatboot")
    opt_bin = ROOT / "src/BGV-Boot-auxradix-opt/build/fatboot"
    old_logs = ROOT / "results"
    rows = []
    for param in MA_PARAMS:
        row = dict(param)
        variants: dict[str, Any] = {}
        common = [
            f"i={param['i']}",
            f"h={param['h']}",
            "newks=1",
            "thick=0",
            f"repeat={args.repeat}",
        ]
        env = {"HELIB_ZZX_CACHE_DIR": str(ROOT / "cache/saved_ZZX")}

        # Ma type-A, the power-of-p auxiliary path.
        log_type_a = out_dir / f"ma_{param['label']}_typeA_repeat{args.repeat}.log"
        variants["ma_typeA"] = run_command(
            [str(base_bin), *common, "t=0", "newbts=1"],
            log_type_a,
            args.timeout,
            env=env,
            force=args.force,
        )
        variants["ma_typeA"].update(parse_extra_log_fields(log_type_a))

        B = variants["ma_typeA"].get("B")
        row["B"] = B
        if B is not None:
            search = search_cleaners(param["p"], int(B), top=args.top)
            row["cleaner_search"] = search
            search_path = out_dir / f"ma_{param['label']}_cleaner_search.json"
            search_path.write_text(json.dumps(search, indent=2, sort_keys=True) + "\n")

        # Ma type-B/default non-power-of-p path.  Reuse the already measured
        # repeat-three log for the largest parameter when possible.
        log_type_b = out_dir / f"ma_{param['label']}_typeB_repeat{args.repeat}.log"
        reused = None
        if param["label"] == "E" and args.repeat == 3:
            reused = copy_existing(
                old_logs / "baseline_BGV_Boot_for_Large_p_i4_repeat3.log",
                log_type_b,
                args.force,
            )
        if reused is not None:
            variants["ma_typeB"] = reused
        elif B is not None and support_type_b_feasible(param["p"], int(B)):
            variants["ma_typeB"] = run_command(
                [str(base_bin), *common, "t=-1", "newbts=1"],
                log_type_b,
                args.timeout,
                env=env,
                force=args.force,
            )
        else:
            variants["ma_typeB"] = {
                "status": "skipped",
                "reason": "non-power-of-p condition (2B+1)^2 < p is false or B unknown",
                "log": str(log_type_b),
            }
        variants["ma_typeB"].update(parse_extra_log_fields(log_type_b))

        # Native HElib/Chen-Han style run.  This can be extremely slow for
        # large-p, so it is opt-in.
        if args.include_native:
            log_native = out_dir / f"ma_{param['label']}_native_repeat{args.repeat}.log"
            variants["helib_native"] = run_command(
                [str(base_bin), *common, "t=0", "newbts=0"],
                log_native,
                args.native_timeout,
                env=env,
                force=args.force,
            )
            variants["helib_native"].update(parse_extra_log_fields(log_native))

        # Our isolated artifact on the same type-B support, if legal.
        search = row.get("cleaner_search") or {}
        chosen = search.get("chosen_practical") or {}
        chosen_aux = chosen.get("aux")
        can_run_ours = B is not None and chosen_aux and support_type_b_feasible(param["p"], int(B))
        if can_run_ours:
            log_sparse = out_dir / f"ma_{param['label']}_ours_aux{chosen_aux}_repeat{args.repeat}.log"
            variants["ours_sparse"] = run_command(
                [str(opt_bin), *common, "t=-1", "newbts=1"],
                log_sparse,
                args.timeout,
                env={**env, "HELIB_EXPLICIT_AUX": str(chosen_aux)},
                force=args.force,
            )
            variants["ours_sparse"].update(parse_extra_log_fields(log_sparse))

            decomp = chosen.get("best_decomposition") or {}
            is_order4 = chosen.get("sqrt_minus_one") and decomp.get("modulus") == 4
            if is_order4:
                log_order4 = out_dir / f"ma_{param['label']}_ours_order4_aux{chosen_aux}_repeat{args.repeat}.log"
                if param["label"] == "E" and chosen_aux == 256 and args.repeat == 3:
                    reused = copy_existing(old_logs / "order4_aux256_repeat3.log", log_order4, args.force)
                    variants["ours_order4"] = reused or {}
                else:
                    variants["ours_order4"] = run_command(
                        [str(opt_bin), *common, "t=-1", "newbts=1"],
                        log_order4,
                        args.timeout,
                        env={
                            **env,
                            "HELIB_EXPLICIT_AUX": str(chosen_aux),
                            "HELIB_AUX_ORDER4_EVAL": "1",
                        },
                        force=args.force,
                    )
                variants["ours_order4"].update(parse_extra_log_fields(log_order4))
            else:
                variants["ours_order4"] = {
                    "status": "skipped",
                    "reason": "chosen cleaner is not an order-four sqrt(-1) decomposition",
                }
        else:
            variants["ours_sparse"] = {
                "status": "skipped",
                "reason": "same-support non-power-of-p path is not valid for this parameter",
            }
            variants["ours_order4"] = {
                "status": "skipped",
                "reason": "same-support non-power-of-p path is not valid for this parameter",
            }

        row["variants"] = variants
        rows.append(row)
    return rows


def run_ntt(args: argparse.Namespace, out_dir: Path) -> list[dict[str, Any]]:
    if not args.include_ntt:
        return [
            {
                **p,
                "status": "not_run",
                "reason": "NTT baseline uses a different power-of-two linear-transform executor; enable with --include-ntt",
            }
            for p in NTT_PARAMS
        ]
    bin_path = Path("/home/user/experiments/baselines/bgv-bootstrapping-with-homomorphic-NTT/build/fatboot")
    rows = []
    for param in NTT_PARAMS:
        row = dict(param)
        common = [
            f"i={param['i']}",
            f"h={param['h']}",
            "t=0",
            "newbts=1",
            "newks=1",
            "thick=0",
            f"repeat={args.repeat}",
        ]
        variants: dict[str, Any] = {}
        log_opt = out_dir / f"ntt_{param['label']}_optimized_repeat{args.repeat}.log"
        variants["ntt_article_optimized"] = run_command(
            [str(bin_path), *common, "baseline=0"],
            log_opt,
            args.timeout,
            force=args.force,
        )
        variants["ntt_article_optimized"].update(parse_extra_log_fields(log_opt))
        B = variants["ntt_article_optimized"].get("B")
        row["B"] = B
        if B is not None:
            search = search_cleaners(param["p"], int(B), top=args.top)
            row["cleaner_search"] = search
            (out_dir / f"ntt_{param['label']}_cleaner_search.json").write_text(
                json.dumps(search, indent=2, sort_keys=True) + "\n"
            )
        log_base = out_dir / f"ntt_{param['label']}_baseline_repeat{args.repeat}.log"
        variants["ntt_article_baseline"] = run_command(
            [str(bin_path), *common, "baseline=1"],
            log_base,
            args.timeout,
            force=args.force,
        )
        variants["ntt_article_baseline"].update(parse_extra_log_fields(log_base))
        variants["ours_ciphertext"] = {
            "status": "not_integrated",
            "reason": (
                "Our order-four evaluator is implemented in the isolated Ma-large-p HElib tree; "
                "the NTT repository uses a different patched HElib/executor, so same-parameter "
                "ciphertext comparison would require a separate port."
            ),
        }
        row["variants"] = variants
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--native-timeout", type=int, default=900)
    parser.add_argument("--top", type=int, default=8)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--include-native", action="store_true")
    parser.add_argument("--include-ntt", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/baseline_matrix")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (ROOT / "cache/saved_ZZX").mkdir(parents=True, exist_ok=True)

    rows = []
    rows.extend(run_ma_large_p(args, args.output_dir))
    rows.extend(run_ntt(args, args.output_dir))

    summary = {
        "repeat": args.repeat,
        "timeout": args.timeout,
        "parameters": rows,
        "non_comparable": [
            {
                "directory": "Bootstrapping_Polyfunctions",
                "reason": "requires Magma; no `magma` binary is available in this environment, so polynomial regeneration cannot be run locally.",
            },
            {
                "directory": "fheanor/openfhe-development/lattigo",
                "reason": "different scheme/library benchmarks; no drop-in HElib BGV digit-extraction pipeline for same ciphertext timing.",
            },
        ],
    }
    json_path = args.output_dir / f"summary_repeat{args.repeat}.json"
    md_path = args.output_dir / f"report_repeat{args.repeat}.md"
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    md_path.write_text(make_markdown(summary) + "\n")
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
