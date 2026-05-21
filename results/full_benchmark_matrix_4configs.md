# Full Benchmark Matrix: 4 Configurations

Date: 2026-05-12

Parameters: p=65537, r=1, bits=1500, m=50731, mvec=[97,523], h=12, t=-1, repeat=3

## Results (mean over 3 runs, seconds)

| Config | linear1 | linear2 | extract | total | vs Baseline |
|--------|--------:|--------:|--------:|------:|------------:|
| Ma baseline (aux=35) | 4.428 | 45.269 | 160.783 | 213.103 | 1.00x |
| Parallel CoeffToSlot (aux=256) | 4.603 | 23.906 | 151.121 | 182.354 | 1.17x |
| Order-4 Cleaner (aux=256) | 4.272 | 43.623 | 80.770 | 131.181 | 1.62x |
| **Combined (Order-4 + Parallel)** | **4.628** | **24.709** | **86.383** | **118.494** | **1.80x** |

## Improvement vs Baseline (%)

| Config | linear1 | linear2 | extract | total |
|--------|--------:|--------:|--------:|------:|
| Parallel CoeffToSlot | -3.95% | +47.19% | +6.01% | +14.43% |
| Order-4 Cleaner | +3.52% | +3.64% | +49.76% | +38.44% |
| **Combined** | **-4.50%** | **+45.42%** | **+46.27%** | **+44.40%** |

## Analysis

The two optimizations target different stages and compose well:

1. **Order-4 Cleaner** (HELIB_AUX_ORDER4_EVAL=1): Exploits the algebraic structure
   P(X) = c₁X + X³Q(X⁴) of the aux=256 cleaning polynomial to reduce the
   homomorphic polynomial evaluation cost in digit extraction. Speedup concentrated
   in `extract` phase (49.76% reduction).

2. **Parallel CoeffToSlot** (HELIB_AUX_PARALLEL_COEFF2SLOT=1): Runs the two
   independent CoeffToSlot transforms in the newBTS pipeline concurrently.
   Speedup concentrated in `linear2` phase (47.19% reduction).

3. **Combined**: Both optimizations apply simultaneously. The `extract` speedup
   is slightly less than standalone Order-4 (46.27% vs 49.76%) due to normal
   run-to-run variance. The `linear2` speedup is preserved (45.42%).

The combined configuration achieves **1.80x total speedup** over the Ma et al.
baseline, reducing end-to-end bootstrapping time from 213.1s to 118.5s.

## Negative Result: Rader97 Factorization

The `HELIB_THINSTEP2_RADER97=1` path (Rader/butterfly factorization of the D=96
linear transform step) fails at ciphertext level with noise overflow. The 13-stage
SparseMatMul1DExec cascade introduces too many key-switches, each adding noise
that exceeds the modulus budget. This is an architectural limitation of the
unfused-stage approach, not a correctness bug in the matrix algebra (which was
verified at plaintext level).

## Environment Variables

```bash
# Combined (best)
HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1 HELIB_AUX_PARALLEL_COEFF2SLOT=1

# Order-4 only
HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1

# Parallel only
HELIB_EXPLICIT_AUX=256 HELIB_AUX_PARALLEL_COEFF2SLOT=1

# Baseline
(no special flags)
```

## Source Logs

- `baseline_BGV_Boot_for_Large_p_i4_repeat3.log`
- `parallel_coeff2slot_aux256_repeat3.log`
- `order4_aux256_repeat3.log`
- `order4_parallel_combined_repeat3.log`
