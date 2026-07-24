#!/usr/bin/env bash
# run_crossover_trend.sh - encrypted butterfly-vs-BSGS crossover across block
# dimensions D = 192, 288, 576, 1152 (Contribution 2 trend experiment).
#
# For each D a ring is used whose native hypercube dimension has size exactly D
# (see encrypted_linear_transform/crossover/results_trend/TREND_SUMMARY.md for
# the ring-construction recipe). Both paths run in the same session on the same
# input ciphertext: the mixed-radix DIF butterfly stage chain (stages<D>.txt,
# plaintext-verified against the row-permuted DFT matrix by derive_trend.py)
# versus HElib's monolithic BSGS MatMul1D over the same dense DxD block; both
# outputs are verified by decryption (0 slot mismatches).
#
# Expected within-row time ratios butterfly/BSGS: 1.25 (192), 1.03 (288),
# 0.61 (576), 0.59 (1152) -- crossover just above D=288. Absolute seconds are
# not comparable across rows (rings differ); the per-row ratio is the metric.
# Reference logs: encrypted_linear_transform/crossover/results_trend/
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CD="$ROOT/encrypted_linear_transform/crossover"
BIN="${ENC_CROSSOVER:-$CD/enc_crossover}"
BITS="${BITS:-600}"
RESULTS="$CD/results_trend"
mkdir -p "$RESULTS"

if [ ! -x "$BIN" ]; then
  echo "building enc_crossover (HElib fork from local/helib_auxradix_opt must be installed; see README Build)"
  g++ -O2 -std=c++17 "$CD/enc_crossover.cpp" -o "$BIN" \
    $(pkg-config --cflags --libs helib 2>/dev/null || echo "-lhelib -lntl -lgmp -lpthread")
fi

run_one() {
  local D=$1 m=$2
  echo "[trend] D=$D m=$m bits=$BITS"
  "$BIN" "$CD/stages$D.txt" "$m" "$BITS" 1 1 | tee "$RESULTS/run_D${D}_m${m}_bits${BITS}.log"
}

run_one 192  17177
run_one 288  14119
run_one 576  51353
run_one 1152 19601

echo "[trend] done: logs in $RESULTS"
