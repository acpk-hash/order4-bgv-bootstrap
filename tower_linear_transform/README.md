# Tower Linear Transform (Contribution 2)

## Status: Theoretically Correct, Practically Limited in BSGS Framework

This directory contains the implementation and experiments for the Galois-structured
linear transform (Section 6 of the paper). The butterfly decomposition is
**mathematically verified** but does not yield wall-clock improvement in HElib's
BSGS framework.

## Experiments

| # | Experiment | Status | Result |
|---|---|---|---|
| 1 | Plaintext Rader prototype | ✅ PASS | 30% fewer operations (ratio 0.7031) |
| 2 | tower_bgv correctness (m=13, p=157) | ✅ PASS | Round-trip identity verified |
| 3 | Tower vs BSGS (naive ZZX, D=96) | ✅ PASS | 1.46x speedup (same arithmetic) |
| 4 | Rader verification in HElib (F_{p^18}) | ✅ PASS | Product == Vandermonde |
| 5 | Full bootstrap (Rader Vandermonde) | ✅ PASS | 216s (≈ baseline 208s) |
| 6 | Factored Rader (4 independent stages) | ❌ FAIL | Noise overflow |

## Key Findings

1. **The Rader decomposition is mathematically correct**: the product of
   conj_NTT × diag(B) × INTT × P_out equals the Step2Matrix Vandermonde,
   verified numerically in HElib's actual slot ring F_{p^18}.

2. **On naive arithmetic (no BSGS)**: the butterfly gives 1.46x speedup
   because it reduces rotations from D=96 to ~13.

3. **On HElib (with BSGS)**: no improvement because:
   - BSGS already reduces the single-matrix application to ~18 rotations
   - The factored form (4 dense stages) needs ~40 rotations total
   - Sequential dense stages cause noise accumulation exceeding modulus chain

4. **Applicable scenarios**: GPU/ASIC implementations where:
   - All rotations are equally expensive (no BSGS benefit)
   - Native sparse-diagonal evaluation is available
   - The O(log D) per-layer count translates to real savings

## Directory Structure

```
tower_linear_transform/
├── tower_bgv_src/          # Standalone implementation (NTL-based)
├── test/                   # Correctness tests (m=13, p=157)
├── bench/                  # Benchmarks (D=96 naive, D=96 fair comparison)
├── src/
│   └── rader_bootstrap.cpp # Rader verification in HElib
├── helib_patch/
│   └── EvalMap_tower.cpp   # Modified HElib EvalMap with Rader integration
├── scripts/
│   └── plaintext_rader_evalmap.py  # Python Rader prototype
├── results/                # Experiment logs
└── run_tower_experiments.sh
```

## Running

```bash
chmod +x run_tower_experiments.sh
./run_tower_experiments.sh
```

## Relation to Peikert-Pepin (2025)

Our work provides the **first concrete instantiation** of the Peikert-Pepin
structured linear transform framework on a practical HElib parameter set.
While Peikert-Pepin gives the theoretical O(log D) framework, they explicitly
leave implementation to future work. We:

1. Identify the condition D | p^d - 1 as the concrete enabler
2. Compute omega_96 in F_{p^18} and verify the Rader product
3. Demonstrate that in BSGS-based libraries, the practical benefit is limited
4. Characterize the applicable scenarios (GPU/ASIC)

This honest analysis is itself a contribution: it clarifies the gap between
the theoretical O(log D) promise and the practical O(√D) reality in current
FHE libraries.
