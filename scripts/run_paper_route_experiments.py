#!/usr/bin/env python3
"""Run the paper-oriented route experiment matrix.

The script intentionally writes only under upload/results/paper_routes.  It
combines the existing HElib baseline logs, the SAT/sparse cleaner experiment,
the structural linear-transform models, and a small set of Fheanor tests that
exercise recent alternative bootstrapping code.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "paper_routes"
FHEANOR = ROOT / "external" / "fheanor"


@dataclass
class CommandResult:
    label: str
    command: list[str]
    returncode: int
    stdout: str
    stderr: str

    @property
    def ok(self) -> bool:
        return self.returncode == 0


def run(label: str, command: list[str], cwd: Path) -> CommandResult:
    proc = subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)
    result = CommandResult(label, command, proc.returncode, proc.stdout, proc.stderr)
    (RESULTS / f"{label}.log").write_text(
        "$ " + " ".join(command) + "\n\n" + proc.stdout + proc.stderr,
        encoding="utf-8",
    )
    return result


def parse_time_line(text: str) -> dict[str, float | str] | None:
    pattern = re.compile(
        r"time for linear1 = ([0-9.]+), linear2 = ([0-9.]+), "
        r"extract = ([0-9.]+), total = ([0-9.]+)"
    )
    matches = pattern.findall(text)
    if not matches:
        return None
    linear1, linear2, extract, total = matches[-1]
    return {
        "linear1": float(linear1),
        "linear2": float(linear2),
        "extract": float(extract),
        "total": float(total),
        "correctness": (
            "pass"
            if "### bts finished, everything ok ###" in text
            or "thin results match" in text
            or "results match" in text
            else "unknown"
        ),
    }


def existing_baselines() -> dict:
    logs = {
        "aux35_direct_ma_baseline_repeat3": ROOT / "results" / "baseline_BGV_Boot_for_Large_p_i4_repeat3.log",
        "aux35_internal_reproduction_repeat3": ROOT / "results" / "upload_aux_default_repeat3.log",
        "aux256_sparse_repeat3": ROOT / "results" / "upload_aux_256_repeat3.log",
        "aux256_parallel_coeff2slot_repeat3": ROOT / "results" / "parallel_coeff2slot_aux256_repeat3.log",
        "aux35_original_repeat1": ROOT / "results" / "upload_aux_default_repeat1.log",
        "aux256_sparse_repeat1": ROOT / "results" / "upload_aux_256_repeat1.log",
    }
    out = {}
    for label, path in logs.items():
        if path.exists():
            out[label] = parse_time_line(path.read_text(encoding="utf-8", errors="replace"))
    return out


def parse_json_from_stdout(result: CommandResult) -> dict:
    return json.loads(result.stdout)


def cargo_summary(result: CommandResult) -> dict:
    text = result.stdout + result.stderr
    passed = re.search(r"test result: ok\. ([0-9]+) passed; 0 failed", text)
    failed = "FAILED" in text or result.returncode != 0
    done_ms = [int(x) for x in re.findall(r"done in ([0-9]+) ms", text)]
    return {
        "ok": result.ok and not failed,
        "passed": int(passed.group(1)) if passed else 0,
        "done_ms_values": done_ms,
        "max_done_ms": max(done_ms) if done_ms else None,
    }


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    commands: list[CommandResult] = []

    commands.append(
        run(
            "sat_aux35_vs_aux256",
            [
                sys.executable,
                "scripts/aux_cleaning_global_sweep.py",
                "--p",
                "65537",
                "--B",
                "17",
                "--mode",
                "list",
                "--aux",
                "35",
                "256",
                "--json",
            ],
            ROOT,
        )
    )

    commands.append(
        run(
            "sat_aux256_feature_search",
            [
                sys.executable,
                "scripts/sat_cleaning_circuit_search.py",
                "--p",
                "65537",
                "--B",
                "17",
                "--aux",
                "256",
                "--degree-limit",
                "64",
                "--feature-set",
                "mixed",
                "--max-features",
                "160",
                "--json",
            ],
            ROOT,
        )
    )

    commands.append(
        run(
            "route_models",
            [sys.executable, "scripts/structural_route_models.py", "--json"],
            ROOT,
        )
    )

    commands.append(
        run(
            "routeB_plaintext_rader97",
            [sys.executable, "scripts/plaintext_rader_evalmap.py", "--ell", "97", "--seed", "7"],
            ROOT,
        )
    )

    cargo_tests = [
        (
            "routeC_fheanor_slots_to_coeffs_thin",
            ["cargo", "test", "lin_transform::pow2::test_slots_to_coeffs_thin", "--", "--nocapture"],
        ),
        (
            "routeC_fheanor_coeffs_to_slots_thin",
            ["cargo", "test", "lin_transform::pow2::test_coeffs_to_slots_thin", "--", "--nocapture"],
        ),
        (
            "routeC_fheanor_bgv_bootstrap17",
            ["cargo", "test", "bgv::bootstrap::test_pow2_bgv_thin_bootstrapping_17", "--", "--nocapture"],
        ),
        (
            "routeC_fheanor_bfv_bootstrap17",
            ["cargo", "test", "bfv::bootstrap::test_pow2_bfv_thin_bootstrapping_17", "--", "--nocapture"],
        ),
        (
            "routeD_fheanor_bounded_digit_retain",
            ["cargo", "test", "poly_eval::digit_extract::test_bounded_digit_retain_poly", "--", "--nocapture"],
        ),
    ]

    if FHEANOR.exists():
        for label, command in cargo_tests:
            commands.append(run(label, command, FHEANOR))

    summary = {
        "existing_helib_baselines": existing_baselines(),
        "sat_sparse_cleaner": parse_json_from_stdout(commands[0]) if commands[0].ok else {"ok": False},
        "sat_feature_search": parse_json_from_stdout(commands[1]) if commands[1].ok else {"ok": False},
        "structural_route_models": parse_json_from_stdout(commands[2]) if commands[2].ok else {"ok": False},
        "routeB_plaintext_rader97": {
            "ok": commands[3].ok,
            "stdout": commands[3].stdout.strip().splitlines(),
        },
        "fheanor_tests": {
            r.label: cargo_summary(r)
            for r in commands[4:]
            if r.label.startswith("routeC_") or r.label.startswith("routeD_")
        },
        "claim_boundaries": {
            "paper_main_claim": "HElib ciphertext result for aux=256 sparse exact cleaner",
            "implementation_choice": "parallel coeffToSlot is wall-clock scheduling, not algebraic complexity",
            "route_A": "model-level only; needs exact EvalMap matrix extraction before ciphertext claim",
            "route_B": "plaintext operation count positive; previous unfused ciphertext executor failed",
            "route_C": "external Fheanor ciphertext tests pass at small p=17; not apples-to-apples with HElib p=65537",
            "route_D": "bounded cleaner tests pass; no CKKS ciphertext replacement implemented locally",
        },
    }

    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))

    failed = [c.label for c in commands if not c.ok]
    if failed:
        print("failed commands: " + ", ".join(failed), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
