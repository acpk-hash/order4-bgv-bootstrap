#!/usr/bin/env bash
# run_po2_prime_sweep.sh - ten-prime OFF/ON sweep on the power-of-two ring
# m = 2^16 at key weight h=26 (matching Set I of the homomorphic-NTT pipeline).
#
# OFF = generic bounded-support PS; ON = order-four evaluator. Every prime
# shares B=13 and digit-polynomial degree 727; expected speedups 1.59-1.74x.
# Reference logs: results/power_of_two/expand/ (and grid/ for the h and
# ring-size scans: i0 h in {12,18,26,40}; i3 = m=2^15, i4 = m=2^14 at h=26).
#
# Requires the patched homomorphic-NTT pipeline (see run_p257_boundary.sh).
# Entries: prime : parameter index i : minimal order-four radix A.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${PO2_FATBOOT:-$ROOT/baselines/bgv-bootstrapping-with-homomorphic-NTT/build/fatboot}"
RESULTS="$ROOT/results/power_of_two/expand"
mkdir -p "$RESULTS"
export HELIB_ZZX_CACHE_DIR="${HELIB_ZZX_CACHE_DIR:-$ROOT/cache/po2_expand}"
mkdir -p "$HELIB_ZZX_CACHE_DIR"

H="${H:-26}"
SWEEP="65537:0:256 1601:5:40 3137:6:56 14401:7:120 4513:8:95 1297:9:36 7057:10:84 13457:11:116 15377:12:124 2521:13:71"

for e in $SWEEP; do
  IFS=: read -r p i A <<< "$e"
  echo "[po2] p=$p (i=$i, A=$A) OFF"
  HELIB_EXPLICIT_AUX="$A" "$BIN" i="$i" h="$H" t=-1 newbts=1 newks=1 thick=0 repeat=3 \
    > "$RESULTS/${p}_OFF.log" 2>&1
  echo "[po2] p=$p (i=$i, A=$A) ON"
  HELIB_EXPLICIT_AUX="$A" HELIB_AUX_ORDER4_EVAL=1 "$BIN" i="$i" h="$H" t=-1 newbts=1 newks=1 thick=0 repeat=3 \
    > "$RESULTS/${p}_ON.log" 2>&1
  echo "[po2] p=$p done"
done

echo "[po2] sweep complete: logs in $RESULTS"
