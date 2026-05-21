#!/usr/bin/env bash
# run_multi_p.sh - Run benchmarks across multiple plaintext primes
# Tests order-4 on applicable primes and baseline on all
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$ROOT/src/BGV-Boot-auxradix-opt/build/fatboot"
BASELINE_BIN="$ROOT/baselines/BGV-Boot-for-Large-p/build/fatboot"
RESULTS="$ROOT/results/multi_p"
mkdir -p "$RESULTS"

export HELIB_ZZX_CACHE_DIR="$ROOT/cache/saved_ZZX"
mkdir -p "$HELIB_ZZX_CACHE_DIR"

REPEAT="${REPEAT:-1}"
TIMEOUT="${TIMEOUT:-1800s}"

echo "================================================================"
echo "  Multi-Prime Benchmark: Order-Four Applicability"
echo "  Testing across p=17, 127, 257, 8191, 65537"
echo "  Repeat=$REPEAT"
echo "================================================================"
echo ""

# Parameter sets:
# i=3: p=17, m=38309 (Set A, order-4 applicable, A=4)
# i=4: p=65537, m=50731 (Set E, order-4 applicable, A=256) [main result]
# i=5: p=127, m=56647 (Set B, NOT applicable)
# i=6: p=257, m=55427 (Set C, order-4 applicable, A=16)
# i=7: p=8191, m=45193 (Set D, NOT applicable)

declare -A PARAM_NAMES
PARAM_NAMES[3]="p=17_m=38309"
PARAM_NAMES[4]="p=65537_m=50731"
PARAM_NAMES[5]="p=127_m=56647"
PARAM_NAMES[6]="p=257_m=55427"
PARAM_NAMES[7]="p=8191_m=45193"

declare -A ORDER4_AUX
ORDER4_AUX[3]="4"      # 4^2=16=-1 mod 17
ORDER4_AUX[4]="256"    # 256^2=65536=-1 mod 65537
ORDER4_AUX[5]=""       # NOT applicable (127≡3 mod 4)
ORDER4_AUX[6]="16"     # 16^2=256=-1 mod 257
ORDER4_AUX[7]=""       # NOT applicable (8191≡3 mod 4)

run_one() {
    local idx="$1"
    local config="$2"
    local logname="$3"
    shift 3
    local logfile="$RESULTS/${logname}.log"
    echo -n "  [$config] "
    if timeout "$TIMEOUT" env "$@" "$BIN" i="$idx" h=12 t=-1 newbts=1 newks=1 thick=0 repeat="$REPEAT" > "$logfile" 2>&1; then
        local timing=$(grep "time for linear" "$logfile" | tail -1)
        if [[ -n "$timing" ]]; then
            echo "$timing"
        else
            echo "completed (no timing line)"
        fi
    else
        echo "FAILED or TIMEOUT"
    fi
}

for idx in 3 4 5 6 7; do
    name="${PARAM_NAMES[$idx]}"
    aux="${ORDER4_AUX[$idx]}"

    echo "--- ${name} (index $idx) ---"

    # Baseline
    run_one "$idx" "baseline" "${name}_baseline"

    # Order-4 (if applicable)
    if [[ -n "$aux" ]]; then
        run_one "$idx" "order-4 (A=$aux)" "${name}_order4" \
            HELIB_EXPLICIT_AUX="$aux" HELIB_AUX_ORDER4_EVAL=1
    else
        echo "  [order-4] NOT APPLICABLE (p≡3 mod 4)"
    fi
    echo ""
done

# Summary
echo "================================================================"
echo "  SUMMARY"
echo "================================================================"
echo ""
printf "  %-25s %-10s %-12s %-12s %-8s\n" "Parameter Set" "Order-4?" "Baseline" "Order-4" "Speedup"
printf "  %-25s %-10s %-12s %-12s %-8s\n" "-------------" "--------" "--------" "-------" "-------"

for idx in 3 4 5 6 7; do
    name="${PARAM_NAMES[$idx]}"
    aux="${ORDER4_AUX[$idx]}"

    bl_total=$(grep "time for linear" "$RESULTS/${name}_baseline.log" 2>/dev/null | tail -1 | grep -oP 'total = \K[\d.]+' || echo "—")

    if [[ -n "$aux" ]]; then
        o4_total=$(grep "time for linear" "$RESULTS/${name}_order4.log" 2>/dev/null | tail -1 | grep -oP 'total = \K[\d.]+' || echo "—")
        applicable="YES (A=$aux)"
        if [[ "$bl_total" != "—" && "$o4_total" != "—" ]]; then
            speedup=$(python3 -c "print(f'{float(\"$bl_total\")/float(\"$o4_total\"):.2f}x')" 2>/dev/null || echo "—")
        else
            speedup="—"
        fi
    else
        o4_total="N/A"
        applicable="NO"
        speedup="N/A"
    fi

    printf "  %-25s %-10s %-12s %-12s %-8s\n" "$name" "$applicable" "${bl_total}s" "${o4_total}s" "$speedup"
done

echo ""
echo "================================================================"
