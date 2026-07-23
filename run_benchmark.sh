#!/usr/bin/env bash
# run_benchmark.sh - Run complete benchmark suite
# Compares: baseline (Ma et al.) vs order-four evaluator
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUR_BIN="$ROOT/src/BGV-Boot-auxradix-opt/build/fatboot"
BASELINE_BIN="$ROOT/baselines/BGV-Boot-for-Large-p/build/fatboot"
RESULTS="$ROOT/results"
mkdir -p "$RESULTS"

REPEAT="${REPEAT:-3}"
TIMEOUT="${TIMEOUT:-1800s}"
COMMON_ARGS=(i=4 h=12 t=-1 newbts=1 newks=1 thick=0 repeat="$REPEAT")

# Check binaries
if [[ ! -x "$OUR_BIN" ]]; then
    echo "ERROR: Our binary not found. Run ./build.sh first."
    exit 1
fi
if [[ ! -x "$BASELINE_BIN" ]]; then
    echo "WARNING: Baseline binary not found at $BASELINE_BIN"
    echo "  Build it with: cd baselines/BGV-Boot-for-Large-p && mkdir -p build && cd build && cmake .. && make"
fi

echo "================================================================"
echo "  Order-Four BGV Bootstrap Benchmark"
echo "  Parameters: p=65537, m=50731, h=12, thin bootstrap"
echo "  Repeat: $REPEAT"
echo "================================================================"
echo ""

# Step 1: SAT search for optimal polynomial
echo "[Step 1] Generating order-four cleaner polynomial via SAT search..."
source "$ROOT/../.venv_sat/bin/activate" 2>/dev/null || true
python3 "$ROOT/scripts/character_projected_cleaner_search.py" \
    --p 65537 --B 17 --aux 256 \
    --output "$RESULTS/cleaner_poly.json" 2>&1 | tail -5
echo ""

# Step 2: Run baseline
if [[ -x "$BASELINE_BIN" ]]; then
    echo "[Step 2] Running baseline (Ma et al., aux=35)..."
    timeout "$TIMEOUT" "$BASELINE_BIN" "${COMMON_ARGS[@]}" \
        > "$RESULTS/baseline_repeat${REPEAT}.log" 2>&1
    echo "  Done. Log: $RESULTS/baseline_repeat${REPEAT}.log"
    grep "time for linear" "$RESULTS/baseline_repeat${REPEAT}.log" | tail -1
else
    echo "[Step 2] SKIPPED (baseline binary not available)"
fi
echo ""

# Step 3: Run sparse cleaner (aux=256, no specialized evaluator)
echo "[Step 3] Running sparse cleaner (aux=256, generic PS)..."
HELIB_EXPLICIT_AUX=256 timeout "$TIMEOUT" "$OUR_BIN" "${COMMON_ARGS[@]}" \
    > "$RESULTS/sparse_aux256_repeat${REPEAT}.log" 2>&1
echo "  Done. Log: $RESULTS/sparse_aux256_repeat${REPEAT}.log"
grep "time for linear" "$RESULTS/sparse_aux256_repeat${REPEAT}.log" | tail -1
echo ""

# Step 4: Run order-four evaluator (aux=256 + specialized evaluator)
echo "[Step 4] Running order-four evaluator (aux=256 + specialized)..."
HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1 timeout "$TIMEOUT" "$OUR_BIN" "${COMMON_ARGS[@]}" \
    > "$RESULTS/order4_aux256_repeat${REPEAT}.log" 2>&1
echo "  Done. Log: $RESULTS/order4_aux256_repeat${REPEAT}.log"
grep "time for linear" "$RESULTS/order4_aux256_repeat${REPEAT}.log" | tail -1
echo ""

# Step 5: Parse and compare
echo "================================================================"
echo "  RESULTS COMPARISON"
echo "================================================================"
python3 "$ROOT/scripts/parse_helib_bootstrap_timings.py" \
    "$RESULTS/baseline_repeat${REPEAT}.log" \
    "$RESULTS/sparse_aux256_repeat${REPEAT}.log" \
    "$RESULTS/order4_aux256_repeat${REPEAT}.log" 2>/dev/null || \
python3 -c "
import re, sys

def parse_log(path):
    try:
        with open(path) as f:
            lines = [l for l in f if 'time for linear' in l]
            if not lines: return None
            last = lines[-1]
            m = re.search(r'linear1 = ([\d.]+).*linear2 = ([\d.]+).*extract = ([\d.]+).*total = ([\d.]+)', last)
            if m: return tuple(float(x) for x in m.groups())
    except: pass
    return None

results = {}
for name, path in [
    ('Baseline (Ma)', '$RESULTS/baseline_repeat${REPEAT}.log'),
    ('Sparse (A=256)', '$RESULTS/sparse_aux256_repeat${REPEAT}.log'),
    ('Order-4 (A=256)', '$RESULTS/order4_aux256_repeat${REPEAT}.log'),
]:
    r = parse_log(path)
    if r:
        results[name] = r
        print(f'  {name:20s}: linear1={r[0]:.2f}s  linear2={r[1]:.2f}s  extract={r[2]:.2f}s  total={r[3]:.2f}s')

if 'Baseline (Ma)' in results and 'Order-4 (A=256)' in results:
    b, o = results['Baseline (Ma)'], results['Order-4 (A=256)']
    print(f'')
    print(f'  Extract speedup: {(b[2]-o[2])/b[2]*100:.1f}%')
    print(f'  Total speedup:   {(b[3]-o[3])/b[3]*100:.1f}%')
"
echo ""
echo "================================================================"
echo "  Benchmark complete. Results in: $RESULTS/"
echo "================================================================"
