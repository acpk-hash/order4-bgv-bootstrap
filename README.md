# Artifact: Order-Four Symmetry in BGV Bootstrapping

Artifact for the Asiacrypt 2026 paper

> **Order-Four Symmetry in BGV Bootstrapping: Faster Digit Extraction for Large Primes**
> Zhenyu Xiong, Mingsheng Wang.

Badges requested: **Available**, **Functional**, **Reproduced**.

Paper references use the numbering of the **camera-ready version** (the PDF
submitted with this artifact). The IACR ePrint full version numbers Tables 1–8
and Sections 1–7 identically and adds appendices with Tables 9 and 10.

---

## Contents

1. [Overview](#1-overview)
2. [Repository layout](#2-repository-layout)
3. [Requirements](#3-requirements)
4. [Getting started (30 minutes)](#4-getting-started-30-minutes)
5. [Claims and how they are checked](#5-claims-and-how-they-are-checked)
6. [One-click comparison with the original baselines](#6-one-click-comparison-with-the-original-baselines)
7. [Step-by-step reproduction](#7-step-by-step-reproduction)
8. [Parameter selection and security](#8-parameter-selection-and-security)
9. [Reading the output](#9-reading-the-output)
10. [Reference reruns and figures](#10-reference-reruns-and-figures)
11. [Known deviations between the paper text and the measurements](#11-known-deviations-between-the-paper-text-and-the-measurements)
12. [Implementation and modifications to public code](#12-implementation-and-modifications-to-public-code)
13. [Troubleshooting](#13-troubleshooting)
14. [Licence](#14-licence)

---

## 1. Overview

BGV bootstrapping for a large plaintext prime `p` spends most of its time in
**digit extraction**, the homomorphic evaluation of a polynomial `P_A` of degree
`(2B+1)² − 2` (Ma et al., EUROCRYPT 2024). The paper shows that if the auxiliary
radix satisfies `A² ≡ −1 (mod p)`, then

```
P_A(AX) + A·P_A(X) = AX,   hence   P_A(X) = ½·X + X³·Q(X⁴),   deg Q = (deg P_A − 3)/4,
```

so only the monomials `X¹` and `X^k`, `k ≡ 3 (mod 4)`, survive, and the
evaluator works on `Q` instead of `P_A`. Our HElib switches between the generic
evaluation of `P_A` (**OFF**, the baseline's evaluator) and the factored one
(**ON**) through the environment variable `HELIB_AUX_ORDER4_EVAL`, so the same
binary measures control and treatment. The original baselines are also built
unmodified from their published sources, so that every comparison can be traced
back to code we did not write.

## 2. Repository layout

The artifact is grouped by role; each group has its own README.

| folder | contents | README |
|---|---|---|
| `1_implementation/` | our code: modified HElib, driver, patches, Section 5 prototype | `1_implementation/README.md` |
| `2_baselines/` | **unmodified** sources of every baseline: HElib `3e337a6`, Ma et al. Large-p `83b54d5`, Ma et al. NTT `800bae0`, and NTL 11.5.1 | `2_baselines/README.md` |
| `3_experiments/` | build, run, one-click baseline comparison, tests, tables, figures | this file, Sec. 5–7 |
| `4_parameters/` | parameter selection: Lattice Estimator `3e48ef4` (pinned copy), scripts, transcripts, **recommended-parameter tables** | `4_parameters/README.md` |
| `5_results/` | logs of the reported runs, our reruns, baseline-comparison output, **figures** | `5_results/README.md` |
| `6_proofs/` | Lean 4 + Mathlib formalisation of the algebraic core (no `sorry`) | `6_proofs/lean/README.md` |

```
3_experiments/
  build_all.sh           builds everything from this repository, offline
  run.sh                 one target per experiment of the paper
  compare_baselines.sh   one-click comparison: original baselines + control groups (Sec. 6)
  verify.sh              automated tests, PASS/FAIL per claim: quick | smoke | compare | full
  make_tables.py         Tables 5, 6, 8 from logs, next to the paper, with --check
  make_figures.py        figures of Tables 4, 5, 6, 8 (5_results/figures/)
  compare_report.py      report, CSV and figures of compare_baselines.sh
  paper_values.json      the paper's numbers and every test tolerance
  tests/                 test_plaintext.py, test_security.py
```

Run-time output goes to `out/` (git-ignored). Every script prints its usage with `-h` or without arguments.

## 3. Requirements

**Hardware.** x86-64 with AVX2. The reported runs used Intel Xeon processors
(Tables 5 and 6: Xeon Gold 6230R, 2.1 GHz).

| level | what runs | wall clock | peak RAM | cores |
|---|---|---|---|---|
| quick | build checks, plaintext algebra, shipped logs, security transcripts, admissibility | 2 min (+ 10 min Docker build) | 2 GB | 1 |
| smoke | quick + one Table 5 row + Section 5 prototype | 15 min | 16 GB | 1 |
| compare | quick + one-click baseline comparison (Sec. 6) | 2.5 h | 32 GB | 2 |
| full | smoke + compare + every ciphertext run of the paper | 13 h | **128 GB** (one row needs 118 GB) | 2 |

Memory per run grows as the slot count falls, not with the ring: on `m = 2^16`,
`p = 2521` (4 slots) needs 118 GB, `p = 65537` (32768 slots) 29 GB (Sec. 7).

**Software.** Docker (recommended). Without Docker: Linux, GCC 11.4, CMake ≥ 3.16,
Python 3.8+ (matplotlib for figures), GMP 6.2.1, and NTL 11.5.1 built from
`2_baselines/ntl-11.5.1.tar.gz` with `./configure SHARED=on` (**default
configuration**, Sec. 13). No network access and no commercial software are
needed. Optional: SageMath 9+ (security recomputation), elan / Lean 4 (proofs).

**Timing.** Results are ratios of two runs on one machine. They are sensitive to
concurrent load, including memory-bandwidth load from other users: run one
experiment at a time on an otherwise idle machine (Sec. 10).

## 4. Getting started (30 minutes)

```sh
unzip order4-bgv-bootstrap-ac26-artifact.zip && cd order4-bgv-bootstrap
#   or: git clone https://github.com/acpk-hash/order4-bgv-bootstrap.git && git checkout ac26-artifact-v1

docker build -t order4-bgv .                                   # ~10 min, offline
mkdir -p out
docker run --rm -v "$PWD/out:/artifact/out" order4-bgv 3_experiments/verify.sh          # ~2 min
docker run --rm -v "$PWD/out:/artifact/out" order4-bgv 3_experiments/verify.sh smoke    # ~15 min, 16 GB
```

The build prints `NTL 11.5.1, NTL_SP_NBITS=60: ok` and ends with a listing of
`build/bin/` (five executables). `verify.sh smoke` should end with:

```
PASS E1    environment: NTL configuration, executables built
PASS P1    plaintext: P_A = X/2 + X^3 Q(X^4) on every (p,A,B) used
PASS L1    shipped logs regenerate Tables 5,6,8 within 0.5%
PASS S1    security at built chain: transcripts equal paper (0.02 bit)
PASS A1    admissibility failures: p=257 abort, p=8191 no radix
PASS C1    ciphertext smoke: po2 p=14401 OFF/ON, decrypts, >= 1.5x
PASS B1    Sec. 5 butterfly: 41 rotations, 2784 slots, 0 mismatches

summary: 7 passed, 0 failed, 0 skipped
```

The report is in `out/verify_report.txt`, each test's output in
`out/tests/<ID>.log`. To see the smoke run next to the paper:
`python3 3_experiments/make_tables.py --logs out/smoke po2` (expect `p = 14401`
near 1.9x with a 1066-bit chain; the paper's 2.13x is discussed in Sec. 11, item 2).

## 5. Claims and how they are checked

Level = the `verify.sh` level that runs the check (Sec. 3).

| ID  | Claim | Paper | Check | Expected | Level |
|-----|---------------------------------------|--------------|---------------------------|--------------------|--------|
| P1  | `P_A = X/2 + X^3 Q(X^4)`: only `k = 1` and `k = 3 mod 4` occur, at most `(N-1)/4 + 1` terms for `N = (2B+1)^2`; at `p = 65537, B = 17`: 612 → 307 terms, `deg Q = 305` | Sec. 4.1–4.2, Thms 2, 4 | `tests/test_plaintext.py` | `PASSED` | quick |
| L1  | Tables 5, 6 and 8 are the means of the shipped logs | Sec. 6.2–6.4 | `make_tables.py --check shipped` | 22 rows within 0.5% | quick |
| S1  | λ* at the built chain, per attack | Table 4, full version Table 10 | `tests/test_security.py` | 51 values within 0.02 bit | quick |
| A1  | Both admissibility conditions fail visibly: `p = 257, h = 12` aborts; `p = 8191` (3 mod 4) has no order-four radix | Sec. 6.3, Thm 3 | `run.sh admissibility` | abort message; no `A` | quick |
| C1  | Order-four digit extraction is about 2x faster on `m = 2^16` | Table 5 | `run.sh smoke` | decrypts; 1066-bit chain; ≥ 1.5x | smoke |
| B1  | Section 5 butterfly on ciphertexts: 41 rotations, correct on all 2784 slots | Sec. 5 | `run.sh butterfly` | 0 mismatches | smoke |
| K1  | Against the original baseline code, with equivalence, radix and evaluator controls | Tables 5 and 8 | `compare_baselines.sh` | all checked ratios PASS (Sec. 6) | compare |
| C2a | Ten primes on `m = 2^16`: 1.84x–2.13x | Table 5, Sec. 6.2 | `run.sh po2` | each row within 15% and ≥ 1.5x | full |
| C2b | Nine general cyclotomic sets: 1.72x–1.88x | Table 6, Sec. 6.3 | `run.sh general` | each row within 15% and ≥ 1.5x | full |
| C2c | `m = 50731`: Ma et al. 213.1 s → ours 131.2 s (1.62x); combined 118.5 s (1.80x) | Table 8, Sec. 6.4 | `run.sh m50731` | each ratio within 15% | full |
| S2  | λ* recomputed with the Lattice Estimator | as S1 | `verify.sh --security-rerun 80B` | within 0.02 bit | optional |
| F1  | The algebraic core is machine-checked | Sec. 6 | `verify.sh --lean` | builds, no `sorry` | optional |

**P1 in detail.** `tests/test_plaintext.py` reuses no code of the implementation.
It builds `P_A` by interpolation over F_p from its definition `P_A(hA + l) = l`
for `|h|, |l| ≤ B`, for the paper's main point and for every `(p, A, B)` that
appears in a shipped log (20 points), and checks B-injectivity of the support,
`deg P_A = N − 2`, the functional equation, the monomial filter, the term bound,
the coefficient `1/2`, the factored form on every support point, equality with
the polynomial files in `5_results/supplementary/polynomials/`, and that HElib
printed the same degree in every run. (At `p = 2521, B = 13` one coefficient
vanishes by accident mod `p`; the count is one below the bound, reported as a note.)

**Not covered by scripts** (supporting data, `5_results/supplementary/`): the
composition with Okada et al. (Sec. 6.4) needs their SEAL code; the full
degree-scaling table (Table 7); the HElib-native and Geelen et al.
rows of Table 8, which the paper does not measure (`2_baselines/README.md`).

## 6. One-click comparison with the original baselines

```sh
3_experiments/compare_baselines.sh           # both experiments, ~2.5 h, 32 GB, 2 cores
3_experiments/compare_baselines.sh m50731    # ~1 h, 17 GB
3_experiments/compare_baselines.sh po2       # ~1.5 h, 32 GB
```

Each experiment is a fixed set of groups at `p = 65537`; every group is one run
of three trials, decrypting and checking every slot.

**`m50731`** (Ma et al. Large-p, Table 8):

| group | executable | radix A | evaluator | role |
|---|---|---|---|---|
| G0 | `baseline_ma_largep` (original code from `2_baselines/`) | 35 | – | baseline |
| G1 | `fatboot_order4` | 35 | OFF | **equivalence control**: our build with the evaluator off must match G0 |
| G2 | `fatboot_order4` | 256 | OFF | **radix control**: the order-four radix alone |
| G3 | `fatboot_order4` | 256 | ON | **treatment** |
| G4 | `fatboot_order4` + parallel CoeffToSlot | 256 | ON | the paper's "combined" row (needs 2 cores) |

**`po2`** (Ma et al. homomorphic NTT, `m = 2^16`, Table 5):

| group | executable | chain | radix A | evaluator | role |
|---|---|---|---|---|---|
| G0 | `baseline_ma_ntt` (original code) | baseline's own, 1314 bits built | 27 | – | baseline |
| G1 | `fatboot_po2` | 1314 | 27 | OFF | **equivalence control** |
| G2 | `fatboot_po2` | 1314 | 256 | OFF | **radix control** |
| G3 | `fatboot_po2` | 1314 | 256 | ON | treatment (Table 5, last column) |
| G4 | `fatboot_po2` | recommended, 1066 | 27 | OFF | Table 5 OFF as measured |
| G5 | `fatboot_po2` | 1066 | 256 | OFF | radix control at the recommended chain |
| G6 | `fatboot_po2` | 1066 | 256 | ON | Table 5 ON |

`3_experiments/compare_report.py` then answers one question per ratio and marks
it PASS/FAIL against a stated range:

| ratio | question | expected |
|---|---|---|
| G0/G1 (total) | do our changes, switched off, leave the baseline unchanged? | 0.95–1.05 |
| G1/G2 (extract) | what does the radix alone do? | informative |
| G2/G3 (extract) | what does the evaluator alone do, at the same radix? | ≥ 1.5 (m50731), ≥ 1.4 (po2) |
| G0/G3, G0/G4 (total) | ours vs the original baseline (Table 8) | ≥ 1.35; paper 1.62, 1.80 |
| G4/G6 (extract, total) | Table 5 as measured | ≥ 1.5 / ≥ 1.1; paper 1.85 / 1.27 |
| G5/G6 (extract) | evaluator alone at the recommended chain | ≥ 1.4 |

Output: `out/compare/{m50731,po2}/<group>.log`, `out/compare/summary.csv`,
`out/compare/report.md` and `out/compare/figures/compare_{m50731,po2}.png`.
Our own run is in `5_results/baseline_comparison/` (Sec. 10).

## 7. Step-by-step reproduction

All commands run inside the container
(`docker run --rm -it -v "$PWD/out:/artifact/out" order4-bgv`) or after a native
`JOBS=8 3_experiments/build_all.sh`. A run without `### bts finished, everything ok ###`
is a failure. `3_experiments/run.sh list` prints the cost of every target.

### 7.1 Build

```sh
JOBS=8 3_experiments/build_all.sh       # or: docker build -t order4-bgv .
```

Five executables in `build/bin/` (Sec. 12): `fatboot_order4`, `fatboot_po2` and
`enc_butterfly` (ours), `baseline_ma_largep` and `baseline_ma_ntt` (original
code), plus HElib installs in `local/`. Sources are copied from `1_implementation/`
and `2_baselines/` into `build/src/` before patching; the repository is never modified.

### 7.2 Plaintext, logs and security (no ciphertexts)

```sh
python3 3_experiments/tests/test_plaintext.py        # P1, 30 s
python3 3_experiments/make_tables.py --check shipped # L1, 1 s
python3 3_experiments/tests/test_security.py         # S1, 1 s
```

### 7.3 Table 5: power-of-two ring `m = 2^16`

```sh
3_experiments/run.sh po2                         # all ten primes, ~2.7 h
3_experiments/run.sh po2 1601 3137 14401         # cheapest three, ~25 min, 15 GB
python3 3_experiments/make_tables.py --logs out --check rerun po2
```

| p | 1297 | 1601 | 2521 | 3137 | 4513 | 7057 | 13457 | 14401 | 15377 | 65537 |
|---|---|---|---|---|---|---|---|---|---|---|
| peak RAM (GB) | 59 | 15 | **118** | 15 | 30 | 59 | 59 | 15 | 59 | 29 |
| wall, OFF + ON (min) | 19 | 8 | 37 | 7 | 11 | 19 | 18 | 8 | 19 | 20 |

### 7.4 Table 6: general cyclotomic rings

```sh
3_experiments/run.sh general                     # all nine sets, ~5.6 h
3_experiments/run.sh general 80B 80A 128D        # cheapest three, ~1.2 h, 13 GB
python3 3_experiments/make_tables.py --logs out --check rerun general
```

| set | 80A | 80B | 80C | 80D | 128A | 128B | 128C | 128D | 128E |
|---|---|---|---|---|---|---|---|---|---|
| peak RAM (GB) | 12 | 8 | 12 | 10 | 61 | 38 | 22 | 13 | 16 |
| wall, OFF + ON (min) | 20 | 20 | 23 | 27 | 61 | 58 | 32 | 30 | 30 |

### 7.5 Table 8: `p = 65537, m = 50731`

```sh
3_experiments/run.sh m50731                      # ma, sparse, order4, combined: ~1 h, 17 GB
3_experiments/run.sh m50731 ma order4            # the headline ratio only, ~30 min
python3 3_experiments/make_tables.py --logs out --check rerun m50731
```

`ma` is the original Ma et al. code; `sparse` is ours at `A = 256` with the
evaluator OFF (not a row of the paper; it separates radix and evaluator);
`order4` is ON; `combined` adds parallel CoeffToSlot and **needs two cores**
(`docker run --cpuset-cpus=0,1 ...`).

### 7.6 Section 6.3 admissibility, Section 5 butterfly

```sh
3_experiments/run.sh admissibility               # <1 min
3_experiments/run.sh butterfly                   # ~3 min, 5 GB
```

### 7.7 Figures

```sh
python3 3_experiments/make_figures.py                                 # shipped logs -> 5_results/figures/
python3 3_experiments/make_figures.py --logs out --out out/figures    # your runs
```

### 7.8 Lean proofs (optional)

```sh
cd 6_proofs/lean && lake exe cache get && lake build   # Lean v4.31.0-rc2, Mathlib 24ef8a8
```

### 7.9 Everything

```sh
3_experiments/verify.sh full        # ~13 h, 128 GB
```

## 8. Parameter selection and security

`4_parameters/README.md` explains how a set is chosen (radix, noise bound,
ring, two LWE instances, minimum over attacks, recomputation at the built chain)
and gives the **recommended parameters** as tables and as
`4_parameters/recommended_parameters.csv`. In short (λ* in bits):

| set | ring | p | m | h | log2 Q built | λ* |
|---|---|---|---|---|---|---|
| po2-80 | power of two | 65537 | 2^16 | 26 | 1066 | 86.7 |
| po2-100 / po2-128 | power of two | 65537 | 2^17 | 26 | (recommended 1676 / 1277) | 106.0 / 131.4 |
| 80-A … 80-D | general | 4513 … 2521 | 37637 … 50731 | 12 | 1524 … 2026 | 80.5 … 82.75 |
| 128-A … 128-E | general | 13457 … 14401 | 65047 … 54359 | 26 | 1407 … 1133 | 129.81 … 131.0 |

The Lattice Estimator is included at the pinned commit
(`4_parameters/lattice-estimator-3e48ef4/`), so
`cd 4_parameters/scripts && sage -python built_chain.py 80B` recomputes a row
with no download (15–40 min); `3_experiments/tests/test_security.py --rerun 80B`
does the same and compares with the paper.

## 9. Reading the output

```
INFO: total bits = 1066, security level = 92.8902        # modulus chain HElib built
log2(aux) = 8.000000                                      # radix A = 256
HELIB_AUX_ORDER4_EVAL enabled: deg(P)=727, terms(P)=183, deg(Q)=181, terms(Q)=182
time for linear1 = 0.51, linear2 = 7.50, extract = 22.80, total = 31.26
### bts finished, everything ok ###                       # all slots decrypted correctly
```

`linear1`, `linear2` are SlotToCoeff and CoeffToSlot, `extract` is digit
extraction (the only stage the method changes), `total` the whole recryption.
`make_tables.py` prints, for each row, the speedup next to the paper's value:

```
po2: Table 5, Sec. 6.2
       p   logQ   OFF(s)    ON(s)  speedup   paper     dev  n
   14401   1066    22.80    12.07    1.89x   2.13x  -11.3%  3
```

`speedup` is OFF/ON of `extract` (Tables 5–6) or Ma/ours of `total` (Table 8);
`dev` its relative deviation from the paper; `n` the trials. **Tolerances** are
in `3_experiments/paper_values.json`: shipped logs must match the paper to 0.5%
(they are its data); reruns must be within 15% and above 1.5x (1.35x for
Table 8), because a ratio of two timings moves by a few per cent on another
machine and by more on a loaded one. Absolute seconds are never checked.

## 10. Reference reruns and figures

We rebuilt the image and reran a subset on the paper's host (one dedicated core
per run, two for "combined", three trials, idle host). Logs:
`5_results/reruns/` (`python3 3_experiments/make_tables.py --logs 5_results/reruns`).

| run | paper | Docker rerun |
|---|---|---|
| Table 5, `p = 14401`, extract OFF → ON | 24.39 → 11.48 s, 2.13x | 22.80 → 12.07 s, 1.89x (Sec. 11, item 2) |
| Table 6, set 80B, extract OFF → ON | 100.05 → 57.40 s, 1.74x | 100.86 → 58.59 s, 1.72x |
| Table 8, Ma et al., total | 213.10 s | 211.75 s |
| Table 8, ours, total | 131.18 s, 1.62x | 134.83 s, 1.57x |
| Table 8, combined, total | 118.49 s, 1.80x | 113.18 s, 1.87x |
| Sec. 6.3, `p = 257, h = 12` | aborts | aborts, same message |
| Sec. 5 butterfly | 41 rotations, 0 mismatches | 41 rotations, 0 mismatches |
| Table 4 λ*, full version Table 10 | nine sets | all 51 values within 0.01 bit |

**Baseline comparison (Sec. 6) on the same host**, Docker image, two
dedicated cores, three trials per group; full report in
`5_results/baseline_comparison/report.md`, per-group logs and `summary.csv`
beside it. All 9 checked ratios pass:

| experiment | ratio | question | measured | paper |
|---|---|---|---|---|
| m50731 | G0/G1, total | original code vs our build with the evaluator OFF | 212.7 / 214.3 s = 0.99x | – |
| m50731 | G1/G2, extract | radix alone (35 → 256) | 1.10x | – |
| m50731 | G2/G3, extract | evaluator alone | 1.76x | – |
| m50731 | G0/G3, total | ours vs original Ma et al. | 1.57x | 1.62x |
| m50731 | G0/G4, total | combined vs original Ma et al. | 1.87x | 1.80x |
| po2 | G0/G1, total | original code vs our build with the evaluator OFF | 54.1 / 53.6 s = 1.01x | – |
| po2 | G1/G2, extract | radix alone (27 → 256), baseline chain | 1.14x | – |
| po2 | G2/G3, extract | evaluator alone, baseline chain | 1.66x | 1.59x (first trial) |
| po2 | G4/G6, extract | Table 5 as measured, 1066-bit chain | 1.94x | 1.85x |
| po2 | G5/G6, extract | evaluator alone, 1066-bit chain | 1.69x | – |
| po2 | G4/G6, total | whole recryption, 1066-bit chain | 1.32x | 1.27x |

The equivalence controls (G0/G1 = 0.99x and 1.01x) show that our executables
with the evaluator OFF behave as the original code. The radix controls
separate the two parts of the gain: the order-four radix alone speeds digit
extraction up by 1.10–1.14x, the evaluator alone by 1.66–1.76x. `5_results/reruns/m50731/under_load/`
keeps two attempts made while another user loaded the host (load average ≈ 56):
single trials of the order-four row ranged from 81.6 s to 132.0 s, and
`linear2`, which the method does not touch, slowed by the same factor.

Figures (PNG and PDF) are in `5_results/figures/` and
`5_results/baseline_comparison/figures/`; `5_results/README.md` lists them.

## 11. Known deviations between the paper text and the measurements

1. **Table 5, OFF radix.** Section 6 states that OFF and ON use the same
   auxiliary radix `A`. This holds for Tables 6 and 8. In Table 5,
   OFF ran Ma et al.'s driver at its default radix `A = 2B+1 = 27` and ON at the
   order-four radix; both build the same chain. Its ratio therefore includes the
   choice of radix as well as the evaluator. The comparison of Sec. 6 measures
   both parts separately (po2: G4/G6 as in the paper, G5/G6 the evaluator alone).
2. **Table 5 runs were concurrent.** OFF and ON of a prime ran at the same
   time, with two to four runs on the host. Alone, `p = 14401` gives 1.89x instead
   of 2.13x; the other nine primes lie at 1.84x–1.95x.
3. **Trials.** The "at 1314–1332" column of Table 5 and the power-of-two rows
   of Table 7 are first trials, not means of three; the means
   differ by up to 0.05x (`5_results/supplementary/power_of_two_first_trials/`).
   Set 80-E and the general-ring rows of Table 7 are single trials.
4. **Compiler flags.** The paper says `-O2`; HElib's `Release` build, used for
   every run, compiles with `-O3 -march=native`.

## 12. Implementation and modifications to public code

| executable | used for | HElib | driver |
|---|---|---|---|
| `fatboot_order4` | Tables 6 and 8 | `1_implementation/helib_order4` | `1_implementation/driver_order4` |
| `fatboot_po2` | Table 5 | HElib `3e337a6` + Ma et al. NTT patch + `patches/helib-ntt-order4.patch` | Ma et al. NTT driver + `patches/ma-ntt-fatboot-order4.patch` |
| `enc_butterfly` | Section 5 | as `fatboot_order4` | `1_implementation/butterfly/enc_butterfly.cpp` |
| `baseline_ma_largep` | Table 8, Ma et al.; G0 | HElib `3e337a6` + `bgvboot-final.patch` | original `83b54d5` |
| `baseline_ma_ntt` | G0 of the po2 comparison | HElib `3e337a6` + `github_20241003.patch` | original `800bae0` |

What we changed in public code (`1_implementation/README.md` has details):

* `1_implementation/patches/helib-auxradix-order4.REFERENCE.patch`: complete
  diff of `helib_order4/` against HElib `3e337a6` + Ma et al.'s patch. The
  evaluator is in `src/extractDigits.cpp` (`HELIB_AUX_ORDER4_EVAL`), the radix
  selection in `src/recryption.cpp` (`HELIB_EXPLICIT_AUX`); the rest is
  experimental code behind variables that no experiment sets, except
  `HELIB_AUX_PARALLEL_COEFF2SLOT` ("combined" row).
* `patches/helib-ntt-order4.patch` (2 files, +122 lines): the evaluator in Ma et
  al.'s power-of-two HElib.
* `patches/ma-ntt-fatboot-order4.patch`: Table 5's parameter sets and the
  options `bits=`, `mainhwt=` in their driver.
* `driver_order4/fatboot.cpp`: sets J–S (Table 6) and `bits=` in the Large-p driver.

| variable | meaning |
|---|---|
| `HELIB_EXPLICIT_AUX=A` | auxiliary radix `A` instead of the driver default `2B+1` |
| `HELIB_AUX_ORDER4_EVAL=1` | factored evaluation; ignored with a message if `A² ≢ −1 (mod p)` |
| `HELIB_AUX_PARALLEL_COEFF2SLOT=1` | parallel CoeffToSlot (Table 8, "combined") |
| `HELIB_ZZX_CACHE_DIR=dir` | cache of digit-extraction polynomials (set by the scripts) |

Driver arguments: `i` parameter-set index, `h` encapsulated-key weight, `bits`
requested chain length (HElib builds a longer one and prints it), `mainhwt`
main-key weight, `repeat` trials. The values of every run are in
`3_experiments/run.sh` and `compare_baselines.sh`.

## 13. Troubleshooting

| symptom | cause and fix |
|---|---|
| `INFO: total bits = 1086` instead of 1066 on `m = 2^16`, then a segmentation fault | NTL built with `NTL_ENABLE_AVX_FFT=on` (50-bit primes). Build NTL from `2_baselines/` with `./configure SHARED=on` only; `3_experiments/check_ntl.cpp` tests this. |
| killed, `std::bad_alloc`, exit code 137 | out of memory; see Sec. 7. Docker Desktop limits container memory by default. |
| a speedup far below the paper, and `linear2` slower than the paper as well | the machine is loaded; `linear2` does not depend on the method. Rerun on an idle machine. |
| `combined` / G4 no faster than `order4` / G3 | it got one core; give it two. |
| `HELIB_AUX_ORDER4_EVAL skipped` | `A` is not an order-four element for this `p`; use the scripts, which know the radix of every prime. |
| G0/G1 outside 0.95–1.05 | the two runs were not on comparable cores (frequency scaling, hyper-threads, load); pin both to the same physical core and rerun with `FORCE=1`. |

## 14. Licence

See `LICENSE`. Our code is MIT. `1_implementation/helib_order4` is a modified
HElib (Apache-2.0, modifications in its `NOTICE`); `1_implementation/driver_order4`
derives from Ma et al. (MIT). `2_baselines/` and `4_parameters/lattice-estimator-3e48ef4/`
keep their own licences; the Ma et al. NTT code has none published and is
included only for reproduction (`2_baselines/bgv-bootstrapping-with-homomorphic-NTT-800bae0/NOTICE.md`).
