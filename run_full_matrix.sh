#!/usr/bin/env bash
# run_full_matrix.sh - Run ALL benchmark configurations
# Including: baseline, sparse, order-4, parallel, combined, AKS variants
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$ROOT/src/BGV-Boot-auxradix-opt/build/fatboot"
BASELINE_BIN="$ROOT/baselines/BGV-Boot-for-Large-p/build/fatboot"
RESULTS="$ROOT/results"
mkdir -p "$RESULTS"

REPEAT="${REPEAT:-3}"
TIMEOUT="${TIMEOUT:-1800s}"
COMMON=(i=4 h=12 t=-1 newbts=1 newks=1 thick=0 repeat="$REPEAT")

export HELIB_ZZX_CACHE_DIR="$ROOT/cache/saved_ZZX"
mkdir -p "$HELIB_ZZX_CACHE_DIR"

echo "================================================================"
echo "  Full Benchmark Matrix: Order-Four BGV Bootstrap"
echo "  p=65537, m=50731, h=12, thin bootstrap, repeat=$REPEAT"
echo "================================================================"

run_config() {
    local name="$1"; shift
    local logfile="$RESULTS/${name}_repeat${REPEAT}.log"
    echo -n "  [$name] "
    if timeout "$TIMEOUT" env "$@" "$BIN" "${COMMON[@]}" > "$logfile" 2>&1; then
        local timing=$(grep "time for linear" "$logfile" | tail -1)
        echo "$timing"
    else
        echo "FAILED (see $logfile)"
    fi
}

# Config 1: Baseline (Ma et al.)
echo ""
echo "[1/6] Baseline (Ma et al., aux=35)..."
if [[ -x "$BASELINE_BIN" ]]; then
    timeout "$TIMEOUT" "$BASELINE_BIN" "${COMMON[@]}" > "$RESULTS/baseline_repeat${REPEAT}.log" 2>&1
    grep "time for linear" "$RESULTS/baseline_repeat${REPEAT}.log" | tail -1
else
    echo "  SKIPPED (baseline binary not found)"
fi

# Config 2: Sparse cleaner only (aux=256)
echo ""
echo "[2/6] Sparse cleaner (aux=256, generic PS)..."
run_config "sparse_aux256" HELIB_EXPLICIT_AUX=256

# Config 3: Order-4 evaluator
echo ""
echo "[3/6] Order-4 evaluator (aux=256 + specialized)..."
run_config "order4_aux256" HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1

# Config 4: Parallel CoeffToSlot
echo ""
echo "[4/6] Parallel CoeffToSlot (aux=256)..."
run_config "parallel_c2s" HELIB_EXPLICIT_AUX=256 HELIB_AUX_PARALLEL_COEFF2SLOT=1

# Config 5: Combined (Order-4 + Parallel) -- BEST
echo ""
echo "[5/6] Combined: Order-4 + Parallel (BEST)..."
run_config "combined_best" HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1 HELIB_AUX_PARALLEL_COEFF2SLOT=1

# Config 6: Combined + AKS BSGS (experimental)
echo ""
echo "[6/6] Combined + AKS BSGS (experimental)..."
run_config "combined_aks" HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1 HELIB_AUX_PARALLEL_COEFF2SLOT=1 HELIB_AUX_LC_AKS_BSGS=1

# Summary
echo ""
echo "================================================================"
echo "  RESULTS SUMMARY"
echo "================================================================"
echo ""
printf "  %-30s %8s %8s %8s %8s\n" "Config" "linear1" "linear2" "extract" "total"
printf "  %-30s %8s %8s %8s %8s\n" "------" "-------" "-------" "-------" "-----"

for logfile in "$RESULTS"/*_repeat${REPEAT}.log; do
    name=$(basename "$logfile" _repeat${REPEAT}.log)
    timing=$(grep "time for linear" "$logfile" 2>/dev/null | tail -1)
    if [[ -n "$timing" ]]; then
        l1=$(echo "$timing" | grep -oP 'linear1 = \K[\d.]+')
        l2=$(echo "$timing" | grep -oP 'linear2 = \K[\d.]+')
        ex=$(echo "$timing" | grep -oP 'extract = \K[\d.]+')
        tot=$(echo "$timing" | grep -oP 'total = \K[\d.]+')
        printf "  %-30s %8s %8s %8s %8s\n" "$name" "$l1" "$l2" "$ex" "$tot"
    fi
done

echo ""
echo "  Logs saved to: $RESULTS/"
echo "================================================================"
