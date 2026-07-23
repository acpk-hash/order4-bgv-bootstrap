# Real Ciphertext Bootstrapping Comparison

Date: 2026-04-27

Parameter set:

```text
fatboot i=4 h=12 t=-1 newbts=1 newks=1 thick=0
p=65537, r=1, m=50731, mvec=[97,523], ords=[96,29]
thin bootstrapping, nslots=2784
I_bound=16.402088, B=17
```

## Claim 1: Sparse Auxiliary-Radix Cleaning

This optimization is integrated in HElib through `HELIB_EXPLICIT_AUX`.

Commands:

```bash
./build/fatboot i=4 h=12 t=-1 newbts=1 newks=1 thick=0 repeat=3

HELIB_EXPLICIT_AUX=256 \
  ./build/fatboot i=4 h=12 t=-1 newbts=1 newks=1 thick=0 repeat=3
```

Logs:

```text
/home/user/experiments/results/aux_radix_helib/aux_default35_repeat3.log
/home/user/experiments/results/aux_radix_helib/aux_256_repeat3.log
/home/user/experiments/results/aux_radix_helib/aux_64_repeat1.log
```

Repeat-3 results:

| aux | algebraic tag | qks bits | poly nonzero terms | linear1 | linear2 | extract | total | correctness |
|---:|---|---:|---:|---:|---:|---:|---:|---|
| 35 | default | 38.373361 | 612 | 4.377296s | 45.154748s | 158.456112s | 210.573619s | pass |
| 256 | `A^2=-1 mod 65537` | 41.244073 | 307 | 4.205747s | 43.463309s | 141.778187s | 191.931511s | pass |

Improvement:

| metric | baseline | optimized | speedup |
|---|---:|---:|---:|
| linear1 | 4.377296s | 4.205747s | 3.92% |
| linear2 | 45.154748s | 43.463309s | 3.75% |
| extract | 158.456112s | 141.778187s | 10.53% |
| total | 210.573619s | 191.931511s | 8.85% |

Control:

```text
A=64 is also a power of two but is not algebraic sqrt(-1).
repeat=1 result: extract=160.716592s, total=211.733621s.
```

Interpretation:

```text
The real ciphertext bootstrapping result supports claim 1.  The speedup tracks
the exact cleaning polynomial sparsity change 612 -> 307.  The qks bound
increases by 2.870712 bits, but correctness and final capacity remain fine.
```

## Claim 2: Structured Linear Transform

Current status:

```text
Plaintext exact Rader prototype: done.
HElib ciphertext executor: not integrated yet.
Real optimized ciphertext timing: not available yet.
```

Plaintext prototype:

```text
/home/user/experiments/scripts/plaintext_rader_evalmap.py
```

For the target `ell=97` block:

| method | finite-field adds | finite-field muls | total ops |
|---|---:|---:|---:|
| dense primitive Vandermonde | 9216 | 9216 | 18432 |
| Rader + mixed-radix NTT | 3840 | 9120 | 12960 |

Ratio:

```text
12960 / 18432 = 0.7031
```

Real ciphertext baseline for future integration:

```bash
HELIB_EXPLICIT_AUX=256 \
  ./build/fatboot i=4 h=12 t=-1 newbts=1 newks=1 thick=0 repeat=1
```

Log:

```text
/home/user/experiments/results/aux_radix_helib/structured_linear_target_i4_aux256_repeat1.log
```

Result:

| setting | linear1 | linear2 | extract | total | correctness |
|---|---:|---:|---:|---:|---|
| current HElib + `A=256` | 4.443499s | 43.832781s | 143.722755s | 194.269403s | pass |

Relevant timers:

| timer | value |
|---|---:|
| `AAA_slotToCoeff` | 4.44344s |
| `AAA_coeffToSlot` | 21.6985s |
| `mul_MatMul1DExec` | 41.6792s / 6 calls |
| `automorph` in `matmul.cpp` | 15.2621s / 151 calls |
| `smartAutomorph` | 13.731s / 27 calls |

Interpretation:

```text
Claim 2 is not yet a real ciphertext speedup claim.  It is currently a
validated exact plaintext algorithm plus a real ciphertext baseline.  The next
required step is to integrate a guarded ell=97 Rader executor into the
ThinStep2Matrix/Step2Matrix path and rerun this same full bootstrapping test.
```
