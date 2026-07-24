#!/usr/bin/env bash
# run_p257_boundary.sh - p=257 boundary demonstration (h=7) in the power-of-two pipeline.
#
# OFF = generic bounded-support Paterson-Stockmeyer; ON = order-four evaluator.
# Radix A=16 (16^2 = 256 = -1 mod 257). Lowering the key weight to h=7 makes the
# support condition hold at p=257: I_bound=6.61 -> B=7 -> (2B+1)^2 = 225 < 257.
#
# Requires the Ma et al. homomorphic-NTT pipeline (ASIACRYPT 2024) with
# patches/power_of_two_order4_port.patch applied; the patch adds parameter
# index i=14 = {m=2^16, p=257, c=3, bits=980, nslots=128}. Point PO2_FATBOOT
# at the resulting fatboot binary if it is not in the default location.
#
# Expected (single thread): OFF extract ~12.6 s, ON ~8.0 s (1.58x);
# ON log shows "HELIB_AUX_ORDER4_EVAL enabled: deg(P)=223, terms(P)=57, deg(Q)=55";
# both verdicts "### bts finished, everything ok ###" (all 128 slots correct).
# Reference logs: results/p257_order4/h7_{OFF,ON}.log, results/p257_order4/SUMMARY.md
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${PO2_FATBOOT:-$ROOT/baselines/bgv-bootstrapping-with-homomorphic-NTT/build/fatboot}"
RESULTS="$ROOT/results/p257_order4"
mkdir -p "$RESULTS"
export HELIB_ZZX_CACHE_DIR="${HELIB_ZZX_CACHE_DIR:-$ROOT/cache/p257_order4}"
mkdir -p "$HELIB_ZZX_CACHE_DIR"

H="${H:-7}"

echo "[p257] OFF (generic bounded-support PS), h=$H"
HELIB_EXPLICIT_AUX=16 "$BIN" i=14 h="$H" t=-1 newbts=1 newks=1 thick=0 repeat=3 \
  > "$RESULTS/h${H}_OFF.log" 2>&1
grep -m1 "extract" "$RESULTS/h${H}_OFF.log" || true

echo "[p257] ON (order-four evaluator), h=$H"
HELIB_EXPLICIT_AUX=16 HELIB_AUX_ORDER4_EVAL=1 "$BIN" i=14 h="$H" t=-1 newbts=1 newks=1 thick=0 repeat=3 \
  > "$RESULTS/h${H}_ON.log" 2>&1
grep -m1 "ORDER4_EVAL enabled" "$RESULTS/h${H}_ON.log" || true

echo "[p257] done: logs in $RESULTS (expected 12.6s -> 8.0s, 1.58x)"
