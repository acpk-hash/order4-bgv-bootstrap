# Encrypted prototype: O(log D) mixed-radix/Rader butterfly CoeffToSlot block (Contribution 2)

Date: 2026-07-22. All results measured on this machine, single thread, HElib fork
`github_order4_cleaner/local/helib_auxradix_opt` (same install that `build_o4/fatboot` links).

## What was built

An end-to-end ENCRYPTED evaluation of the ThinEvalMap Step2 (CoeffToSlot) dimension-0
block for the paper's main ring: p=65537, m=50731=97x523 (Set E), D=96, slot ring
F_p[X]/G(X) with deg G=18, hypercube dims {96,29}, both native ("good").

Pipeline (derived and plaintext-verified in `derive_butterfly.py` against the
instrumented HElib dump `github_order4_cleaner/step2_matrix_dim0_sz96.txt`, which is
the exact matrix HElib applies; dump model verified: S[i][j] = zeta^(i*38*5^(-j) mod 97)):

```
w   = INTT3.INTT2.INTT1( Dmid * NTT3.NTT2.NTT1( v ) )      6 butterfly stages
out = w with dim0-slot 48 := DC = sum_j v[j]               Rader DC fixup
y_true[i] = out[dlog_5(i)] (i>=1),  y_true[0] = out[48]    fixed output permutation
```

- NTT chain: mixed radices [6,4,4] (96 = 6*4*4), decimation-in-frequency,
  wrap-correct diagonal supports:
  - stage 1 (radix 6, stride 16): offsets {0,16,32,48,64,80} - 6 diagonals (subgroup, wrap-free)
  - stage 2 (radix 4, stride 4):  offsets {0,+-4,+-8,+-12}  - 7 diagonals
  - stage 3 (radix 4, stride 1):  offsets {0,+-1,+-2,+-3}   - 7 diagonals
- INTT chain: same supports in reverse order. Middle diagonal Dmid (Rader kernel
  NTT'd, kernel = row 1 of the true Step2 matrix) folded into NTT stage 3 for free.
- Rotation budget: 34 nontrivial rotate1D in the 6 stages + 7 for the DC total-sum
  = 41 rotations, over 17 distinct rotation amounts
  ({16,32,48,64,80} u {4,8,12,84,88,92} u {1,2,3,93,94,95}); DC amounts are a subset.
- Depth: 7 constant-mult levels (6 stages + DC masking).
- The residual Rader permutation (t -> 5^t mod 97 on dim-0 coordinates) is fixed and
  known; it can be absorbed into the constants of the adjacent dim-1 Step2 stage /
  downstream maps (permutations commute with dim-1 MatMul and compose into any
  adjacent linear transform for free), so it is not evaluated homomorphically.
  Verification is done against the TRUE block output composed with this fixed permutation.

## Encrypted results (m=50731, real Set E ring)

`enc_butterfly.cpp` builds the Set E context (gens {48117,5239}, ords {96,29}, c=3),
encrypts 2784 random slots of F_p^18, evaluates the chain, decrypts, and compares
every slot against an independent NTL zz_pE matvec of the dump matrix.

| run | keygen | chain wall | rotations | capacity in -> out (consumed) | verified |
|---|---|---|---|---|---|
| bits=600  | add1DMatrices (all 1D KS) | 27.1 s | 41 | 569.3 -> 377.1 (192.1 bits) | PASS 0/2784 mismatches |
| bits=600  | addSome1DMatrices (HElib default) | 38.7 s | 41 | 569.3 -> 376.3 (193.0 bits) | PASS 0/2784 |
| bits=1500 | add1DMatrices | 59.1 s | 41 | 1470.4 -> 1276.4 (194.0 bits) | PASS 0/2784 |

Per-stage capacity consumption is ~21-27 bits/stage (constant-mult level + rotations),
uniform across stages; see run logs.

## Baseline (same ring, same block, same input ciphertext)

Monolithic BSGS MatMul1D (HElib MatMul1DExec) of the identical 96x96 block:

| run | precomp | mul wall | capacity consumed | verified |
|---|---|---|---|---|
| bits=600  | 19 s | 6.4 s  | 55.5 bits | YES (transposed orientation) |
| bits=1500 | 27 s | 14.9 s | 57.2 bits | YES |

Reference from full bootstrapping logs (`results/multi_p/p65537_off.log`): Set E
linear2 (CoeffToSlot, BOTH dims) = 43.3 s, 62.3 bits.

## Honest comparison / findings

1. The encrypted butterfly chain is CORRECT on the real ring: first encrypted
   execution of the paper's Contribution-2 decomposition, exact (0 slot errors).
2. At D=96 in HElib it does NOT beat BSGS: 59.1 s vs 14.9 s (bits=1500) and
   194 vs 57 capacity bits. Reasons, now with numbers: 41 vs ~18 automorphisms, and
   7 vs 1 constant-mult levels (this quantifies the paper's "cumulative noise"
   concern: ~3.4x capacity of the monolithic map, but far from overflow -
   the chain fits comfortably even at bits=600).
3. Key-switching gap: with HElib's default `addSomeIDMatrices` keygen the chain still
   runs (smartAutomorph multi-hop) but costs +43% wall time (38.7 s vs 27.1 s at
   bits=600) at equal noise. Generating all 1D KS matrices (add1DMatrices) removes
   the gap; no non-1D Galois elements are needed - all 17 rotation amounts are native
   dim-0 rotations because the Rader reindexing makes HElib's dim-0 coordinate the
   dlog coordinate.
4. Data provenance finding: the precomputed `tower_stage_coefficients.json`
   (stage1 {0,16,...,80}, stage2 {0..15}) is NOT usable for a cyclic (rotate1D)
   evaluation: exhaustive convention testing (32 composition conventions, plus
   permutation/twist/diagonal-multiset analysis) shows its 2-stage product agrees
   with a twisted Vandermonde only on the diagonal band offset in [0,16] and is
   corrupted on wrapped diagonals - the generation ignored block-boundary wrap
   (correct cyclic stages need the +- offset pairs used here). The stage data used
   by this prototype was therefore re-derived from scratch and re-verified.

## Files

- `derive_butterfly.py`   - derivation + plaintext verification vs HElib dump; writes `butterfly_stages.json`
- `make_stages_txt.py`    - flattens JSON to `stages.txt` for C++
- `enc_butterfly.cpp`     - the encrypted prototype + BSGS baseline (build: g++ -O2 -std=c++17 -march=native enc_butterfly.cpp -I$L/include -L$L/lib -lhelib -lntl -lgmp -lpthread, L=github_order4_cleaner/local/helib_auxradix_opt)
- `verify_chain_plaintext.py` - initial (negative) check of the old JSON stage data
- `run_bits600.log`, `run_bits600_someks.log`, `run_bits1500.log` - raw measurements
