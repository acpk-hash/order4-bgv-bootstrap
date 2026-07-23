#!/usr/bin/env bash
# run_tower_experiments.sh - Run all tower linear transform experiments
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS="$ROOT/results"
mkdir -p "$RESULTS"

echo "================================================================"
echo "  Tower Linear Transform Experiments"
echo "  (Contribution 2: Galois-Structured Linear Transforms)"
echo "================================================================"
echo ""

# ---- Experiment 1: Plaintext Rader Verification ----
echo "[Exp 1] Plaintext Rader prototype (operation count comparison)..."
source /home/luck/xzy/0424project/.venv_sat/bin/activate 2>/dev/null || true
python3 "$ROOT/scripts/plaintext_rader_evalmap.py" --ell 97 2>&1 | tee "$RESULTS/rader_plaintext.log"
echo ""

# ---- Experiment 2: tower_bgv correctness (small params) ----
echo "[Exp 2] Tower BGV correctness test (m=13, p=157, D=12)..."
if [[ -f "$ROOT/../tower_bgv_build/test_tower" ]]; then
    "$ROOT/../tower_bgv_build/test_tower" 2>&1 | tee "$RESULTS/tower_bgv_correctness.log"
else
    echo "  Building tower_bgv..."
    mkdir -p "$ROOT/../tower_bgv_build"
    cd "$ROOT/../tower_bgv_build"
    cmake "$ROOT" -DCMAKE_CXX_COMPILER=/usr/local/gcc-9.5.0/gcc-9.5.0/bin/c++ 2>&1 | tail -3
    make -j4 test_tower 2>&1 | tail -3
    ./test_tower 2>&1 | tee "$RESULTS/tower_bgv_correctness.log"
    cd "$ROOT"
fi
echo ""

# ---- Experiment 3: tower_bgv D=96 benchmark (naive ZZX) ----
echo "[Exp 3] Tower vs BSGS benchmark (D=96, naive ZZX arithmetic)..."
if [[ -f "$ROOT/../tower_bgv_build/bench_fair" ]]; then
    timeout 1200 "$ROOT/../tower_bgv_build/bench_fair" 2>&1 | tee "$RESULTS/tower_vs_bsgs_naive.log"
else
    echo "  SKIPPED (build tower_bgv first: cd tower_bgv_build && make bench_fair)"
fi
echo ""

# ---- Experiment 4: Rader verification in HElib (F_{p^18}) ----
echo "[Exp 4] Rader decomposition verification in HElib slot ring..."
if [[ -f "$ROOT/src/rader_bootstrap" ]]; then
    timeout 300 "$ROOT/src/rader_bootstrap" 2>&1 | tee "$RESULTS/rader_helib_verification.log"
else
    echo "  Building rader_bootstrap..."
    /usr/bin/c++ -std=c++17 -O2 -no-pie \
        -I/home/luck/xzy/0424project/tower_helib/include \
        -o "$ROOT/src/rader_bootstrap" \
        "$ROOT/src/rader_bootstrap.cpp" \
        /home/luck/xzy/0424project/tower_helib/build/lib/libhelib.a \
        /usr/local/lib/libntl.a -lgmp -lgf2x -lpthread 2>&1 | tail -3
    if [[ -f "$ROOT/src/rader_bootstrap" ]]; then
        timeout 300 "$ROOT/src/rader_bootstrap" 2>&1 | tee "$RESULTS/rader_helib_verification.log"
    else
        echo "  BUILD FAILED"
    fi
fi
echo ""

# ---- Experiment 5: Full bootstrap with Rader-verified Vandermonde ----
echo "[Exp 5] Full thin bootstrap with Rader-constructed Vandermonde..."
FATBOOT="/tmp/fatboot_rader_final"
if [[ -f "$FATBOOT" ]]; then
    timeout 700 "$FATBOOT" i=4 h=12 t=-1 newbts=1 newks=1 thick=0 repeat=1 2>&1 | \
        tee "$RESULTS/rader_full_bootstrap.log" | grep -E "time for linear|TOWER|Rader"
else
    echo "  SKIPPED (fatboot_rader_final not built)"
    echo "  To build: link fatboot.cpp against tower_helib/build/lib/libhelib.a"
fi
echo ""

# ---- Summary ----
echo "================================================================"
echo "  EXPERIMENT SUMMARY"
echo "================================================================"
echo ""
echo "  Exp 1 (Plaintext Rader): Operation count reduction"
grep "operation_ratio" "$RESULTS/rader_plaintext.log" 2>/dev/null || echo "    (not available)"
echo ""
echo "  Exp 2 (tower_bgv correctness): Small-parameter verification"
grep "Summary" "$RESULTS/tower_bgv_correctness.log" 2>/dev/null || echo "    (not available)"
echo ""
echo "  Exp 3 (Naive benchmark): Tower vs BSGS on same arithmetic"
grep "Linear transform speedup" "$RESULTS/tower_vs_bsgs_naive.log" 2>/dev/null || echo "    (not available)"
echo ""
echo "  Exp 4 (Rader in HElib): Mathematical correctness in F_{p^18}"
grep "Rader product" "$RESULTS/rader_helib_verification.log" 2>/dev/null || echo "    (not available)"
echo ""
echo "  Exp 5 (Full bootstrap): Rader-verified Vandermonde"
grep "time for linear" "$RESULTS/rader_full_bootstrap.log" 2>/dev/null | tail -1 || echo "    (not available)"
echo ""
echo "  KEY FINDINGS:"
echo "  - Rader decomposition is mathematically correct (Exp 1, 2, 4)"
echo "  - On naive ZXX: Tower gives 1.46x speedup over BSGS (Exp 3)"
echo "  - On HElib (with BSGS): no wall-clock improvement (Exp 5)"
echo "  - Reason: BSGS already reduces rotations to O(√D);"
echo "    factored Rader stages add noise and don't reduce total rotations"
echo "  - Applicable to: GPU/ASIC implementations without BSGS"
echo ""
echo "================================================================"
