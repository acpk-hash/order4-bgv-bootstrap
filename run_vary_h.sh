#!/usr/bin/env bash
# run_vary_h.sh - Vary encapsulated key weight h on p=65537
# Tests how B (noise bound) affects the order-4 speedup
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$ROOT/src/BGV-Boot-auxradix-opt/build/fatboot"
RESULTS="$ROOT/results/vary_h"
mkdir -p "$RESULTS"

export HELIB_ZZX_CACHE_DIR="$ROOT/cache/saved_ZZX"
mkdir -p "$HELIB_ZZX_CACHE_DIR"

REPEAT="${REPEAT:-1}"
TIMEOUT="${TIMEOUT:-1800s}"

echo "================================================================"
echo "  Varying h (key weight) on p=65537, m=50731"
echo "  h affects B (noise bound), which affects polynomial degree"
echo "================================================================"
echo ""

# h values to test: 12 (default), 16, 20, 24
for h in 12 16 20 24; do
    echo "--- h=$h ---"

    # Compute expected B
    # B = ceil(I_bound) where I_bound ≈ 0.5 + beta * sqrt(phi(m)/m * h * 2^k/3) * 0.5
    # For h=12: B=17, for h=16: B≈20, for h=20: B≈22, for h=24: B≈24

    # Baseline
    echo -n "  [baseline] "
    timeout "$TIMEOUT" "$BIN" i=4 h="$h" t=-1 newbts=1 newks=1 thick=0 repeat="$REPEAT" \
        > "$RESULTS/h${h}_baseline.log" 2>&1
    grep "time for linear" "$RESULTS/h${h}_baseline.log" | tail -1 || echo "FAILED"

    # Order-4
    echo -n "  [order-4]  "
    HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1 \
    timeout "$TIMEOUT" "$BIN" i=4 h="$h" t=-1 newbts=1 newks=1 thick=0 repeat="$REPEAT" \
        > "$RESULTS/h${h}_order4.log" 2>&1
    grep "time for linear" "$RESULTS/h${h}_order4.log" | tail -1 || echo "FAILED"

    echo ""
done

echo "================================================================"
echo "  SUMMARY: Effect of h on Order-4 Speedup"
echo "================================================================"
printf "  %-5s %-12s %-12s %-8s\n" "h" "Baseline" "Order-4" "Speedup"
printf "  %-5s %-12s %-12s %-8s\n" "---" "--------" "-------" "-------"

for h in 12 16 20 24; do
    bl=$(grep "time for linear" "$RESULTS/h${h}_baseline.log" 2>/dev/null | tail -1 | grep -oP 'total = \K[\d.]+' || echo "—")
    o4=$(grep "time for linear" "$RESULTS/h${h}_order4.log" 2>/dev/null | tail -1 | grep -oP 'total = \K[\d.]+' || echo "—")
    if [[ "$bl" != "—" && "$o4" != "—" ]]; then
        sp=$(python3 -c "print(f'{float(\"$bl\")/float(\"$o4\"):.2f}x')")
    else
        sp="—"
    fi
    printf "  %-5s %-12s %-12s %-8s\n" "$h" "${bl}s" "${o4}s" "$sp"
done
echo ""
echo "  Note: larger h → larger B → higher degree → more speedup from order-4"
echo "================================================================"
