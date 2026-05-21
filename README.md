# Order-Four Character Filter for Large-Prime BGV Bootstrapping

Artifact for the paper "Efficient Large-Prime BGV Bootstrapping via Algebraic Structure Exploitation" (AsiaCrypt 2026 submission).

## Main Result

For HElib thin bootstrapping with `p=65537, m=50731, 2784 slots`:

| Configuration | Extract | Total | Speedup |
|---|---:|---:|---:|
| Ma et al. baseline (aux=35) | 160.8s | 213.1s | — |
| Sparse cleaner (aux=256) | 147.9s | 199.4s | 6.4% |
| **Order-four evaluator (aux=256)** | **80.8s** | **131.2s** | **38.4%** |

## Key Idea

When the auxiliary radix `A` satisfies `A² ≡ -1 (mod p)` (e.g., `A=256` for `p=65537`):
- The cleaner polynomial satisfies `P_A(AX) + A·P_A(X) = AX`
- Nonzero monomials are restricted to `k=1` and `k ≡ 3 (mod 4)`
- Term count drops from 612 → 307 (49.8% reduction)
- Factored form: `P_A(X) = (1/2)X + X³·Q(X⁴)` with `deg(Q) = 305`
- Evaluation cost: `O(√(d/4))` instead of `O(√d)`

## Directory Structure

```
order4bgv/
├── src/
│   ├── BGV-Boot-auxradix-opt/    # Modified fatboot driver
│   └── HElib_auxradix_opt/       # Modified HElib with order-4 evaluator
├── tower_linear_transform/        # Contribution 2: Galois-structured LT
│   ├── tower_bgv_src/            # Standalone butterfly NTT implementation
│   ├── test/                     # Correctness tests (m=13, p=157)
│   ├── bench/                    # Naive ZZX benchmarks (D=96)
│   ├── src/rader_bootstrap.cpp   # Rader verification in HElib F_{p^18}
│   ├── helib_patch/              # Modified EvalMap.cpp with Rader
│   ├── scripts/                  # Plaintext Rader prototype
│   └── run_tower_experiments.sh  # All tower experiments
├── baselines/                     # Symlinks to baseline implementations
│   ├── BGV-Boot-for-Large-p/     # Ma et al. (2024)
│   ├── HElib/                    # Upstream HElib
│   ├── subring_bts/              # Ma et al. subring (2025)
│   └── Bootstrapping_Polyfunctions/  # Geelen et al. (2022)
├── scripts/
│   ├── verify_order4.py          # Quick polynomial verification
│   ├── character_projected_cleaner_search.py  # Full SAT search
│   └── ...
├── results/                       # Benchmark logs and comparisons
├── build.sh                       # Build script
└── run_benchmark.sh               # Full benchmark suite (Contribution 1)
```

## Quick Start

```bash
# 1. Build (requires NTL, GMP, pthreads)
./build.sh

# 2. Verify polynomial properties
python3 scripts/verify_order4.py

# 3. Run full benchmark (takes ~15 minutes)
./run_benchmark.sh
```

## Environment Variables

| Variable | Description |
|---|---|
| `HELIB_EXPLICIT_AUX=256` | Use A=256 as auxiliary radix |
| `HELIB_AUX_ORDER4_EVAL=1` | Enable specialized order-four evaluator |
| `HELIB_ZZX_CACHE_DIR=path` | Cache directory for precomputed polynomials |

## SAT Search

The polynomial search uses Z3 (via `.venv_sat`):

```bash
source ../.venv_sat/bin/activate
python3 scripts/character_projected_cleaner_search.py --p 65537 --B 17 --aux 256
```

## Baselines

| Paper | Directory | Parameters |
|---|---|---|
| Ma et al. "Large-p" (2024) | `baselines/BGV-Boot-for-Large-p` | Same (p=65537, m=50731) |
| Geelen et al. "Polyfunctions" (2022) | `baselines/Bootstrapping_Polyfunctions` | Different p |
| Ma et al. "Subring" (2025) | `baselines/subring_bts` | Power-of-two cyclotomics |
| HElib upstream | `baselines/HElib` | Reference implementation |

## Reproducing Paper Results

```bash
# Contribution 1: Order-four character filter (digit extraction)
REPEAT=3 ./run_benchmark.sh

# Contribution 2: Tower linear transform (Rader decomposition)
cd tower_linear_transform && ./run_tower_experiments.sh
```

### Contribution 2 Results Summary

| Experiment | Status | Key Result |
|---|---|---|
| Plaintext Rader (ell=97) | ✅ | 30% fewer ops (ratio 0.7031) |
| Butterfly correctness (m=13, p=157) | ✅ | Round-trip PASS (4/4 tests) |
| Tower vs BSGS (naive ZZX, D=96) | ✅ | 1.46x speedup |
| Rader in HElib (F_{p^18}) | ✅ | Product == Vandermonde |
| Full bootstrap (Rader matrix) | ✅ | 216s (≈ baseline 208s) |
| Factored 4-stage Rader | ❌ | Noise overflow (inherent limitation) |

See `tower_linear_transform/README.md` for detailed analysis.

## License

Apache-2.0 (same as HElib)
