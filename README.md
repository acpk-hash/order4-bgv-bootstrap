# Artifact: Order-Four Symmetry in BGV Bootstrapping

Artifact for the Asiacrypt 2026 paper

> **Order-Four Symmetry in BGV Bootstrapping: Faster Digit Extraction for Large Primes**
> Zhenyu Xiong, Mingsheng Wang.

Badges requested: **Available**, **Functional**, **Reproduced**.

References to the paper use the numbering of the **camera-ready version** (the PDF
submitted with this artifact). The IACR ePrint full version has an extra table
at the front and an appendix, so its numbers are given in brackets, e.g.
"Table 4 [5]".

---

## Contents

1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Getting started (30 minutes)](#3-getting-started-30-minutes)
4. [Claims and how they are checked](#4-claims-and-how-they-are-checked)
5. [Step-by-step reproduction](#5-step-by-step-reproduction)
6. [Reading the output](#6-reading-the-output)
7. [Reference reruns](#7-reference-reruns)
8. [Known deviations between the paper text and the measurements](#8-known-deviations-between-the-paper-text-and-the-measurements)
9. [Implementation and modifications to public code](#9-implementation-and-modifications-to-public-code)
10. [Repository layout](#10-repository-layout)
11. [Troubleshooting](#11-troubleshooting)
12. [Licence](#12-licence)

---

## 1. Overview

BGV bootstrapping for a large plaintext prime `p` spends most of its time in
**digit extraction**, the homomorphic evaluation of a polynomial `P_A` of degree
`|S_A| − 2 = (2B+1)² − 2` (Ma et al., EUROCRYPT 2024). The paper shows that if
the auxiliary radix satisfies `A² ≡ −1 (mod p)`, then

```
P_A(AX) + A·P_A(X) = AX,   hence   P_A(X) = ½·X + X³·Q(X⁴),   deg Q = (deg P_A − 3)/4,
```

so only the monomials `X¹` and `X^k` with `k ≡ 3 (mod 4)` survive and the
evaluator can work on `Q` instead of `P_A`. The implementation is a modified
HElib in which the environment variable `HELIB_AUX_ORDER4_EVAL=1` switches the
digit-extraction evaluator from the generic Paterson–Stockmeyer evaluation of
`P_A` (OFF) to the factored evaluation (ON). Everything else in the bootstrapping
pipeline is unchanged, so the OFF/ON ratio isolates the contribution.

The artifact contains:

| component | what it supports |
|---|---|
| modified HElib, four executables, Dockerfile | the ciphertext measurements of Sections 5 and 6 |
| `ae/run.sh`: one target per experiment, with the exact command lines of the reported runs | Tables 4, 5, 6 [5, 6, 8]; Sec. 5; Sec. 6.3 |
| `ae/verify.sh`: automated tests with PASS/FAIL per claim | all claims in Section 4 below |
| `results/paper/`: raw logs of every reported ciphertext run | Tables 4–6 [5, 6, 8] |
| `security/`: Lattice Estimator scripts, pinned version, transcripts | Tables 2–3 [3–4], App. D [Table 10] |
| `lean/`: Lean 4 + Mathlib formalisation, zero `sorry` | Theorems 2–4 |

## 2. Requirements

**Hardware.** x86-64 with AVX2. The reported runs used Intel Xeon processors
(Tables 4 and 5: Intel Xeon Gold 6230R, 2.1 GHz). Time and memory per level:

| level | what runs | wall clock | peak RAM | cores |
|---|---|---|---|---|
| quick | build checks, plaintext, logs, security transcripts, admissibility | 2 min (+10 min Docker build) | 2 GB | 1 |
| smoke | quick + one Table 4 row + Section 5 prototype | 15 min | 16 GB | 1 |
| subset | three rows of each table (Sec. 5.3–5.5) | 3 h | 60 GB | 1 (2 for one row) |
| full | every ciphertext run in the paper | 10 h | **128 GB** (one row needs 118 GB) | 1 (2 for one row) |

Memory per run grows as the slot count falls, not with the ring: on `m = 2^16`
the prime `p = 2521` (4 slots) needs 118 GB and `p = 65537` (32768 slots) 29 GB.
Section 5 lists the requirement of every run.

**Software.** Docker (recommended), or: Linux, GCC 11.4, CMake ≥ 3.16, Python 3.8+,
GMP 6.2.1, NTL 11.5.1 in its **default configuration** (Sec. 11), git, and
network access to github.com during the build. Optional: SageMath 9+ (security
recomputation), elan / Lean 4 (formal proofs). No commercial software is needed.

**Timing.** Results are ratios of two runs on the same machine. They are
sensitive to concurrent load, including memory-bandwidth load from other
users: run one experiment at a time on an otherwise idle machine (Sec. 7).

## 3. Getting started (30 minutes)

```sh
git clone https://github.com/acpk-hash/order4-bgv-bootstrap.git
cd order4-bgv-bootstrap && git checkout ac26-artifact-v1

docker build -t order4-bgv .                 # ~10 min: NTL, 3 x HElib, 4 executables
mkdir -p ae/out
docker run --rm -v "$PWD/ae/out:/artifact/ae/out" order4-bgv ae/verify.sh          # ~2 min
docker run --rm -v "$PWD/ae/out:/artifact/ae/out" order4-bgv ae/verify.sh smoke    # ~15 min, 16 GB
```

The build prints `NTL 11.5.1, NTL_SP_NBITS=60: ok` and ends with a listing of
`build/bin/`. `ae/verify.sh smoke` should print:

```
PASS E1    environment: NTL configuration, executables built
PASS P1    plaintext: P_A = X/2 + X^3 Q(X^4) on every (p,A,B) used
PASS L1    shipped logs regenerate Tables 4-6 within 0.5%
PASS S1    security at built chain: transcripts equal paper (0.02 bit)
PASS A1    admissibility failures: p=257 abort, p=8191 no radix
PASS C1    ciphertext smoke: po2 p=14401 OFF/ON, decrypts, >= 1.5x
PASS B1    Sec. 5 butterfly: 41 rotations, 2784 slots, 0 mismatches

summary: 7 passed, 0 failed, 0 skipped
```

The report is written to `ae/out/verify_report.txt` and the output of each test
to `ae/out/tests/<ID>.log`. To see the smoke run next to the paper:

```sh
python3 ae/make_tables.py --logs ae/out/smoke po2
```

Expect `p = 14401` at about 1.9x (paper 2.13x; see Sec. 8, item 2), with the
chain built at 1066 bits.

## 4. Claims and how they are checked

Test scripts are in `ae/` (`run.sh`, `make_tables.py`) and `ae/tests/`. Level
is the `ae/verify.sh` level that runs the test (Sec. 2).

| ID  | Claim | Paper | Check | Expected | Level |
|-----|---------------------------------------|--------------|---------------------------|--------------------|--------|
| P1  | `P_A = X/2 + X^3 Q(X^4)`: only `k = 1` and `k = 3 mod 4` occur, at most `(N-1)/4 + 1` terms for `N = (2B+1)^2`; at `p = 65537, B = 17`: 612 → 307 terms, `deg Q = 305` | Sec. 4.1–4.2, Thms 2, 4 | `test_plaintext.py` | `PASSED` | quick |
| L1  | Tables 4–6 [5, 6, 8] are the means of the shipped logs | Sec. 6.2–6.4 | `make_tables.py --check shipped` | 22 rows within 0.5% | quick |
| S1  | λ* at the built chain, per attack | Table 3 [4], App. D [Tab. 10] | `test_security.py` | 51 values within 0.02 bit | quick |
| A1  | Both admissibility conditions fail visibly: `p = 257, h = 12` aborts; `p = 8191` (3 mod 4) has no order-four radix | Sec. 6.3, Thm 3 | `run.sh admissibility` | abort message; no `A` | quick |
| C1  | Order-four digit extraction is about 2x faster on `m = 2^16` | Table 4 [5] | `run.sh smoke` | decrypts; 1066-bit chain; ≥ 1.5x | smoke |
| B1  | Section 5 butterfly on ciphertexts: 41 rotations, correct on all 2784 slots | Sec. 5 | `run.sh butterfly` | 0 mismatches | smoke |
| C2a | Ten primes on `m = 2^16`: 1.84x–2.13x | Table 4 [5], Sec. 6.2 | `run.sh po2` | each row within 15% and ≥ 1.5x | full |
| C2b | Nine general cyclotomic sets: 1.72x–1.88x | Table 5 [6], Sec. 6.3 | `run.sh general` | each row within 15% and ≥ 1.5x | full |
| C2c | `m = 50731`: Ma et al. 213.1 s → ours 131.2 s (1.62x); combined 118.5 s (1.80x) | Table 6 [8], Sec. 6.4 | `run.sh m50731` | each ratio within 15% | full |
| S2  | λ* recomputed with the Lattice Estimator | as S1 | `verify.sh --security-rerun 80B` | within 0.02 bit | optional |
| F1  | The algebraic core is machine-checked | Sec. 6 | `verify.sh --lean` | builds, no `sorry` | optional |

**P1 in detail.** `ae/tests/test_plaintext.py` does not reuse any code from the
implementation. It builds `P_A` by interpolation over F_p from its definition
`P_A(hA + l) = l` for `|h|, |l| ≤ B`, for the paper's main point and for every
`(p, A, B)` that appears in a shipped log (20 points), and checks:
B-injectivity of the support; `deg P_A = N − 2`; the functional equation
`P_A(AX) + A·P_A(X) = AX`; the monomial filter; the term bound; the leading
coefficient `1/2`; agreement of the factored form on every support point; equality,
coefficient by coefficient, with the polynomial files in `results/`; and that
HElib printed the same degree in every run. At `p = 2521, B = 13` one
coefficient vanishes by accident mod `p`, so the count is one below the bound;
the test reports this as a note.

`ae/verify.sh full` runs C2a–C2c on top of the smoke level and checks every row
(test `C2` in the report).

**Not covered by scripts** (supporting data, kept for completeness): the
composition with Okada et al. (Sec. 6.4) needs their SEAL implementation
(`okada_comparison/` holds our harness and logs); the degree-scaling table of
the full version (Table 7) is in `results/degree_scaling/` and
`results/paper/degscale_general/`; the HElib-native and Geelen et al. rows of
Table 6 [8] were not measured by us, as the paper states.

## 5. Step-by-step reproduction

All commands run inside the container
(`docker run --rm -it -v "$PWD/ae/out:/artifact/ae/out" order4-bgv`) or in a
native build after `JOBS=8 ae/build_all.sh`. Every run decrypts and compares
all slots; a run without `### bts finished, everything ok ###` is a failure.
`ae/run.sh list` prints the cost of every target.

### 5.1 Build

```sh
JOBS=8 ae/build_all.sh            # or: docker build -t order4-bgv .
```

This produces `build/bin/{fatboot_order4, fatboot_po2, fatboot_ma_largep, enc_butterfly}`
(Sec. 9). Upstream code is fetched at pinned commits; transient GitHub errors
are retried up to eight times.

### 5.2 Plaintext, logs and security (no ciphertexts)

```sh
python3 ae/tests/test_plaintext.py        # P1, 30 s
python3 ae/make_tables.py --check shipped # L1, 1 s
python3 ae/tests/test_security.py         # S1, 1 s
```

### 5.3 Table 4 [5]: power-of-two ring `m = 2^16`

```sh
ae/run.sh po2                          # all ten primes, ~2.7 h
ae/run.sh po2 1601 3137 14401          # cheapest three, ~25 min, 15 GB
python3 ae/make_tables.py --logs ae/out --check rerun po2
```

| p | 1297 | 1601 | 2521 | 3137 | 4513 | 7057 | 13457 | 14401 | 15377 | 65537 |
|---|---|---|---|---|---|---|---|---|---|---|
| peak RAM (GB) | 59 | 15 | **118** | 15 | 30 | 59 | 59 | 15 | 59 | 29 |
| wall, OFF + ON (min) | 19 | 8 | 37 | 7 | 11 | 19 | 18 | 8 | 19 | 20 |

### 5.4 Table 5 [6]: general cyclotomic rings

```sh
ae/run.sh general                      # all nine sets, ~5.6 h
ae/run.sh general 80B 80A 128D         # cheapest three, ~1.2 h, 13 GB
python3 ae/make_tables.py --logs ae/out --check rerun general
```

| set | 80A | 80B | 80C | 80D | 128A | 128B | 128C | 128D | 128E |
|---|---|---|---|---|---|---|---|---|---|
| peak RAM (GB) | 12 | 8 | 12 | 10 | 61 | 38 | 22 | 13 | 16 |
| wall, OFF + ON (min) | 20 | 20 | 23 | 27 | 61 | 58 | 32 | 30 | 30 |

### 5.5 Table 6 [8]: `p = 65537, m = 50731`

```sh
ae/run.sh m50731                       # ma, sparse, order4, combined: ~1 h, 17 GB
ae/run.sh m50731 ma order4             # the headline ratio only, ~30 min
python3 ae/make_tables.py --logs ae/out --check rerun m50731
```

`ma` is the unmodified Ma et al. baseline at its radix `A = 35`; `sparse` is ours
at `A = 256` with the evaluator OFF (not a row of the paper; it separates the
effect of the radix from that of the evaluator); `order4` is ON; `combined`
adds HElib's parallel CoeffToSlot and **needs two cores**
(e.g. `docker run --cpuset-cpus=0,1 ...`).

### 5.6 Section 6.3 admissibility and Section 5 butterfly

```sh
ae/run.sh admissibility                # <1 min
ae/run.sh butterfly                    # ~3 min, 5 GB
```

### 5.7 Security recomputation (optional, SageMath)

```sh
git clone https://github.com/malb/lattice-estimator third_party/src/lattice-estimator
git -C third_party/src/lattice-estimator checkout 3e48ef4
python3 ae/tests/test_security.py --rerun 80B 128C       # ~20 min per set
```

`security/README.md` describes the remaining security scripts (Tables 2–3
[3–4], including the rotated hybrid of ePrint 2026/279), which take several
hours.

### 5.8 Lean proofs (optional)

```sh
cd lean && lake exe cache get && lake build       # Lean v4.31.0-rc2, Mathlib 24ef8a8
```

### 5.9 Everything

```sh
ae/verify.sh full        # ~10 h, 128 GB; or run the subsets of 5.3-5.5
```

## 6. Reading the output

A bootstrapping log prints one line per trial and then the driver's mean:

```
INFO: total bits = 1066, security level = 92.8902        # modulus chain HElib built
HELIB_AUX_ORDER4_EVAL enabled: deg(P)=727, terms(P)=183, deg(Q)=181, terms(Q)=182
time for linear1 = 0.51, linear2 = 7.50, extract = 22.80, total = 31.26
### bts finished, everything ok ###                       # all slots decrypted correctly
```

`linear1` and `linear2` are SlotToCoeff and CoeffToSlot, `extract` is digit
extraction (the only stage the method changes), and `total` is the whole
recryption. `ae/make_tables.py` averages the per-trial lines and prints, for
each row, the speedup next to the paper's value:

```
po2: Table 4, Sec. 6.2 [full version: Table 5, Sec. 6.2]
       p   logQ   OFF(s)    ON(s)  speedup   paper     dev  n
   14401   1066    22.80    12.07    1.89x   2.13x  -11.3%  3
```

`speedup` is OFF/ON of `extract` for Tables 4–5 [5–6], and Ma/ours of `total`
for Table 6 [8]; `dev` is its relative deviation from the paper; `n` is the
number of trials. **Tolerances** (in `ae/paper_values.json`): the shipped logs
must match the paper to 0.5%, since they are the paper's data. Reruns must be
within 15% and above 1.5x (1.35x for Table 6 [8]), because a ratio of two
timings moves by a few per cent on another machine and by more on a loaded one.
Absolute seconds are never checked.

## 7. Reference reruns

We rebuilt the Docker image from scratch and reran a subset on the paper's host
(one dedicated core per run, two for "combined", three trials, idle host). The
logs are in `results/ae_rerun_docker/`; reproduce this table with
`python3 ae/make_tables.py --logs results/ae_rerun_docker`.

| run | paper | Docker rerun |
|---|---|---|
| Table 4 [5], `p = 14401`, extract OFF → ON | 24.39 → 11.48 s, 2.13x | 22.80 → 12.07 s, 1.89x (Sec. 8, item 2) |
| Table 5 [6], set 80B, extract OFF → ON | 100.05 → 57.40 s, 1.74x | 100.86 → 58.59 s, 1.72x |
| Table 6 [8], Ma et al., total | 213.10 s | 211.75 s |
| Table 6 [8], ours, total | 131.18 s, 1.62x | 134.83 s, 1.57x |
| Table 6 [8], combined, total | 118.49 s, 1.80x | 113.18 s, 1.87x |
| Sec. 6.3, `p = 257, h = 12` | aborts | aborts, same message |
| Sec. 5 butterfly | 41 rotations, 0 mismatches | 41 rotations, 0 mismatches |
| Table 3 [4] λ*, App. D [Table 10] | nine sets | all 51 values within 0.01 bit |

`results/ae_rerun_docker/m50731/under_load/` keeps two Table 6 [8] attempts
made while another user loaded the host (load average ≈ 56): single trials of
the order-four row ranged from 81.6 s to 132.0 s of extraction, and `linear2`,
which the method does not touch, slowed by the same factor.

## 8. Known deviations between the paper text and the measurements

1. **Table 4 [5], OFF radix.** Section 6 states that OFF and ON use the same
   auxiliary radix `A`. This holds for Tables 5 and 6 [6, 8]. In Table 4 [5],
   OFF ran Ma et al.'s driver at its default radix `A = 2B+1 = 27` and ON at
   the order-four radix of the prime; both build the same chain (`INFO: total
   bits` in `results/paper/po2/`). Table 4 [5] therefore compares against the
   baseline as published, and its ratio includes the choice of radix as well as
   the evaluator. Table 6 [8] separates the two (`sparse` vs `order4`).
   `ae/run.sh po2` reproduces exactly what was measured.
2. **Table 4 [5] runs were concurrent.** OFF and ON of a prime ran at the same
   time, with two to four runs on the host. Run alone, `p = 14401` gives 1.89x
   instead of 2.13x; the other nine primes lie at 1.84x–1.95x.
3. **Trials.** The "at 1314–1332" column of Table 4 [5] and the power-of-two
   rows of the full version's Table 7 are first trials, not means of three;
   the means differ by up to 0.05x (`results/power_of_two/`). Set 80-E and the
   general-ring rows of the full version's Table 7 are single trials.
4. **Compiler flags.** The paper says `-O2`; HElib's `Release` build, used for
   every run, compiles with `-O3 -march=native`.

## 9. Implementation and modifications to public code

| executable | used for | HElib | driver |
|---|---|---|---|
| `fatboot_order4` | Tables 5, 6 [6, 8] (ours) | `src/HElib_auxradix_opt` | `src/BGV-Boot-auxradix-opt` |
| `fatboot_po2` | Table 4 [5] | HElib `3e337a6` + Ma et al. NTT patch + `third_party/patches/helib-ntt-order4.patch` | Ma et al. NTT driver `800bae0` + `third_party/patches/ma-ntt-fatboot-order4.patch` |
| `fatboot_ma_largep` | Table 6 [8], Ma et al. row | HElib `3e337a6` + Ma et al. `bgvboot-final.patch` | upstream `83b54d5`, unmodified |
| `enc_butterfly` | Section 5 | as `fatboot_order4` | `enc_butterfly.cpp` in `encrypted_linear_transform/` |

Upstream sources (HElib `3e337a6`, `msh086/BGV-Boot-for-Large-p` `83b54d5`,
`msh086/bgv-bootstrapping-with-homomorphic-NTT` `800bae0`) are cloned at these
commits by `ae/build_all.sh`. The NTT repository has no licence, so only our
diff against it is distributed.

**Our changes:**

* `third_party/patches/helib-auxradix-order4.REFERENCE.patch` is the complete
  diff of `src/HElib_auxradix_opt` against HElib `3e337a6` with Ma et al.'s
  patch applied. The order-four evaluator is in `src/extractDigits.cpp`
  (`HELIB_AUX_ORDER4_EVAL`) and the radix selection in `src/recryption.cpp`
  (`HELIB_EXPLICIT_AUX`). The rest of that diff is experimental code behind
  environment variables that are off by default and used by no table.
* `third_party/patches/helib-ntt-order4.patch` (2 files, +122 lines) ports the
  evaluator to Ma et al.'s power-of-two HElib.
* `third_party/patches/ma-ntt-fatboot-order4.patch` adds the parameter sets of
  Table 4 [5] and the options `bits=` and `mainhwt=` to their driver.
* `src/BGV-Boot-auxradix-opt/fatboot.cpp` adds sets J–S (Table 5 [6]) and
  `bits=` to the Ma et al. Large-p driver.

**Configuration.** Environment variables read by the modified HElib:

| variable | meaning |
|---|---|
| `HELIB_EXPLICIT_AUX=A` | auxiliary radix `A` instead of the driver default `2B+1` |
| `HELIB_AUX_ORDER4_EVAL=1` | factored evaluation; ignored with a message if `A² ≢ −1 (mod p)` |
| `HELIB_AUX_PARALLEL_COEFF2SLOT=1` | parallel CoeffToSlot (Table 6 [8], "combined") |
| `HELIB_ZZX_CACHE_DIR=dir` | cache of digit-extraction polynomials (set by `ae/run.sh`) |

Driver arguments: `i` parameter-set index, `h` encapsulated-key weight, `bits`
requested chain length (HElib builds a longer chain and prints it), `mainhwt`
main-key weight, `repeat` number of trials. The values for every row are in
`ae/run.sh`.

## 10. Repository layout

```
README.md                    this file
Dockerfile, .dockerignore    reproducible environment
ae/
  build_all.sh               builds build/bin/*
  run.sh                     one target per experiment (ae/run.sh list)
  verify.sh                  automated tests: quick | smoke | full
  make_tables.py             Tables 4-6 [5, 6, 8] from logs, with --check
  paper_values.json          the paper's numbers and the test tolerances
  tests/                     test_plaintext.py, test_security.py
  check_ntl.cpp              NTL configuration check
results/
  paper/                     raw logs of the reported runs:
                             po2/ general/ m50731/ admissibility/ degscale_general/
  ae_rerun_docker/           our reruns in the Docker image (Sec. 7)
  *.txt, *.json, ...         polynomial files and earlier development logs
src/HElib_auxradix_opt/      modified HElib (Apache-2.0, see its NOTICE)
src/BGV-Boot-auxradix-opt/   driver for Tables 5-6 [6, 8] (MIT, Ma et al.)
third_party/patches/         our diffs against pinned upstream code
security/                    Lattice Estimator scripts and transcripts (security/README.md)
  built_chain.py             per-attack security at the built chain
lean/                        Lean 4 formalisation (lean/README.md)
encrypted_linear_transform/  Section 5 prototype (its README.md)
okada_comparison/            Section 6.4 harness and logs
```

The other top-level scripts and `scripts/` are development history and are not
needed for evaluation. Logs in `results/paper/` are unedited except that the
home directory of the measurement host is replaced by `<HOST>`.

## 11. Troubleshooting

| symptom | cause and fix |
|---|---|
| `INFO: total bits = 1086` instead of 1066 on `m = 2^16`, then a segmentation fault | NTL was built with `NTL_ENABLE_AVX_FFT=on` (50-bit primes). Rebuild NTL with `./configure SHARED=on` only; `ae/check_ntl.cpp` tests this. |
| killed, `std::bad_alloc`, or exit code 137 | out of memory; see the RAM tables in Sec. 5. Docker Desktop limits container memory by default. |
| `clone of https://github.com/... failed` during the build | transient network error; the script retries eight times. Rerun `ae/build_all.sh`, which reuses completed clones. |
| a speedup far below the paper, and `linear2` slower than the paper as well | the machine is loaded; `linear2` does not depend on the method. Rerun on an idle machine. |
| `combined` row no faster than `order4` | it was given one core; give it two. |
| `HELIB_AUX_ORDER4_EVAL skipped` | the radix is not an order-four element for this `p`; use `ae/run.sh`, which knows the radix of every prime. |

## 12. Licence

See `LICENSE`. Our code is MIT. `src/HElib_auxradix_opt` is a modified HElib
(Apache-2.0, modifications listed in its `NOTICE`); `src/BGV-Boot-auxradix-opt`
derives from Ma et al. (MIT). `scripts/install_magma_local.sh` and
`scripts/run_polyfunctions_magma_generation.sh` concern an unused baseline that
needs the commercial Magma system; no result depends on them.
