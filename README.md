# Order-Four Symmetry in BGV Bootstrapping — Artifact

Artifact for

> Zhenyu Xiong, Mingsheng Wang. **Order-Four Symmetry in BGV Bootstrapping:
> Faster Digit Extraction for Large Primes.** Asiacrypt 2026.

The paper replaces the bounded-support digit-extraction polynomial `P_A` of
Ma et al. (EUROCRYPT 2024) by its order-four factorisation
`P_A(X) = ½X + X³·Q(X⁴)` whenever the auxiliary radix satisfies `A² ≡ −1 (mod p)`,
and evaluates `Q` instead of `P_A` inside HElib thin bootstrapping. This
repository contains the implementation (a modified HElib), the drivers and exact
command lines behind every ciphertext table in Section 7, the raw logs of the
reported runs, the security-estimation scripts and transcripts, and a Lean 4
formalisation of the algebraic core.

**Everything needed for evaluation is reachable from this README.** Directories
not listed in the map below are development history and are not needed.

---

## 1. Quick start (≈ 25 min, 16 GB RAM)

```sh
docker build -t order4-bgv .                    # ~10 min: NTL, 3 x HElib, 3 drivers
docker run --rm -v $PWD/ae/out:/artifact/ae/out order4-bgv ae/run.sh smoke
docker run --rm -v $PWD/ae/out:/artifact/ae/out order4-bgv ae/run.sh admissibility
python3 ae/make_tables.py --logs ae/out/smoke 5     # your smoke run beside the paper
python3 ae/make_tables.py                       # all tables from the shipped logs, no build needed
```

`smoke` runs one Table 5 configuration (`p = 14401` on `m = 2^16`) with the
order-four evaluator off and on, one trial each, decrypts, and checks every
slot. Expected: both runs print `### bts finished, everything ok ###`, the chain
is built at 1066 bits, and digit extraction is roughly 2x faster with the
evaluator on. A single trial in the Docker image on our host gave
23.62 s → 12.24 s (1.93x); three trials gave 1.89x. The paper's 2.13x for this
row comes from concurrent runs (Sec. 7, item 4).

## 2. Claims and where they are checked

| # | Claim (paper location) | Command | Expected | Cost |
|---|---|---|---|---|
| C1 | Factored form: `P_A` has 307 instead of 612 terms, `deg Q = 305` at `p = 65537, B = 17` (Sec. 5, Thm. 2) | `python3 scripts/verify_order4.py` | `Verification PASSED`, 612 → 307 | 1 min |
| C2 | Algebraic core machine-checked (Sec. 7 preamble) | `cd lean && lake exe cache get && lake build` | build succeeds, no `sorry` | 10 min + Mathlib cache |
| C3 | Table 5: digit extraction on `m = 2^16`, ten primes, 1.84x–2.13x, mean 1.91x | `ae/run.sh tab5 [p ...]` then `ae/make_tables.py --logs ae/out 5` | speedup per prime about 1.85x–1.95x; see Sec. 7, item 4 for `p = 14401` | 2.7 h, up to 118 GB |
| C4 | Table 6: nine general cyclotomic sets, 1.72x–1.88x | `ae/run.sh tab6 [set ...]` then `ae/make_tables.py --logs ae/out 6` | speedup per set within ~5% (our Docker rerun of 80B: 1.72x vs 1.74x) | 5.6 h, up to 61 GB |
| C5 | Table 7: `m = 50731`, Ma et al. 213.1 s → ours 131.2 s (1.62x), combined 118.5 s (1.80x) | `ae/run.sh tab7` then `ae/make_tables.py --logs ae/out 7` | ratios ±5% | 1 h, 17 GB |
| C6 | Two witnessed admissibility failures (Sec. 7.3) | `ae/run.sh admissibility` | `p too small to hold the overflow part` at `p = 257, h = 12`; no order-four radix for `p = 8191` | 1 min |
| C7 | Table 2 λ* at the built chain and the per-attack table (Sec. 7.1, App. "Per-attack detail") | `security/built_chain.py <set>` | the paper's columns to 0.01 bit | 15–40 min per set |
| C8 | Tables 1–2 security figures, rotated hybrid of ePrint 2026/279 (Sec. 7.1) | see `security/README.md` | `security/results/*.jsonl` | hours |
| C9 | Okada et al. composition, p = 257 (Sec. 7.4) | see `okada_comparison/` | logs there | needs their SEAL code |
| C10 | Butterfly linear transform, 41 rotations, 0 mismatches (Sec. 6) | see `encrypted_linear_transform/README.md` | logs there | 1 min per run |

Absolute times depend on the machine; the paper's claims are the **ratios**
between runs on one machine. `ae/make_tables.py` prints each ratio beside the
paper's value and its relative deviation. Running it with no arguments
tabulates the shipped logs of the reported runs; every row agrees with the
paper to 0.3%.

All runs are single-threaded except the "combined" row of Table 7, which turns
on HElib's parallel CoeffToSlot (`HELIB_AUX_PARALLEL_COEFF2SLOT=1`) and runs the
two CoeffToSlot transforms concurrently; give that run at least two cores.
Timings are sensitive to other load on the machine (memory bandwidth as much as
CPU): run one experiment at a time on an otherwise idle host.

## 3. Hardware and software

* x86-64 with AVX2. The reported runs used Intel Xeon processors; Tables 5
  and 6 were measured on an Intel Xeon Gold 6230R (2.1 GHz, 256 GB).
* RAM: 16 GB for `smoke`, `tab7` and the cheap rows; 64 GB for all of Table 6;
  128 GB for Table 5 row `p = 2521` (peak 118 GB). The other Table 5 rows need
  15–59 GB. Memory grows as the slot count falls, not with the ring.
* Disk: 3 GB for the image.
* Docker image: Ubuntu 22.04, GCC 11.4, CMake 3.22, NTL 11.5.1 (default
  configuration, built from source), GMP 6.2.1. These are the versions on the
  measurement host.

**NTL must use its default configuration.** With `NTL_ENABLE_AVX_FFT=on`, NTL's
single-precision primes shrink from 60 to 50 bits, HElib builds a different
modulus chain (1086 instead of 1066 bits on `m = 2^16`), and the `m = 2^16`
driver crashes. The Dockerfile checks this (`ae/check_ntl.cpp`).

Without Docker: install the packages above, then `JOBS=8 ae/build_all.sh`.

## 4. What is built

`ae/build_all.sh` produces three executables in `build/bin/`:

| executable | used for | HElib | driver |
|---|---|---|---|
| `fatboot_order4` | Tables 6 and 7 (ours) | `src/HElib_auxradix_opt` | `src/BGV-Boot-auxradix-opt` |
| `fatboot_po2` | Table 5 | HElib `3e337a6` + Ma et al. NTT patch + `third_party/patches/helib-ntt-order4.patch` | Ma et al. NTT `fatboot.cpp` `800bae0` + `third_party/patches/ma-ntt-fatboot-order4.patch` |
| `fatboot_ma_largep` | Table 7, row "Ma et al." | HElib `3e337a6` + Ma et al. `bgvboot-final.patch` | upstream `83b54d5`, unmodified |

The two Ma et al. repositories are cloned at the pinned commits during the
build. `bgv-bootstrapping-with-homomorphic-NTT` carries no licence, so only our
diff against it is distributed.

**Modifications to public code**, as required by the call for artifacts:

* `third_party/patches/helib-auxradix-order4.REFERENCE.patch` is the full
  diff of `src/HElib_auxradix_opt` against HElib `3e337a6` with Ma et al.'s
  patch applied. The order-four evaluator is `src/extractDigits.cpp`
  (`HELIB_AUX_ORDER4_EVAL`) and the radix selection is `src/recryption.cpp`
  (`HELIB_EXPLICIT_AUX`). Everything else in that diff is behind environment
  variables that are off by default and is not used by any table.
* `third_party/patches/helib-ntt-order4.patch` (2 files, +122 lines) ports the
  same evaluator to Ma et al.'s power-of-two HElib.
* `third_party/patches/ma-ntt-fatboot-order4.patch` adds parameter sets for the
  ten primes of Table 5 and the `bits=` / `mainhwt=` options to their driver.
* `src/BGV-Boot-auxradix-opt/fatboot.cpp` adds parameter sets J–S (Table 6)
  and `bits=` to the Ma et al. Large-p driver.

## 5. Configuration

Environment variables read by the modified HElib:

| variable | meaning |
|---|---|
| `HELIB_EXPLICIT_AUX=A` | use auxiliary radix `A` instead of the driver default `2B+1` |
| `HELIB_AUX_ORDER4_EVAL=1` | evaluate `P_A` as `½X + X³Q(X⁴)`; ignored, with a message, when `A² ≢ −1` |
| `HELIB_AUX_PARALLEL_COEFF2SLOT=1` | parallel CoeffToSlot (Table 7, "combined" only) |
| `HELIB_ZZX_CACHE_DIR=dir` | cache for digit-extraction polynomials (set by `ae/run.sh`) |

Driver arguments: `i` parameter-set index, `h` encapsulated-key weight,
`bits` requested chain length, `mainhwt` main-key weight, `repeat` trials.
`ae/run.sh` contains the exact values for every row. Note that HElib builds a
longer chain than the `bits` requested; the built chain is printed as
`INFO: total bits = ...`.

## 6. Reading the output

Each run prints one line per trial and then the mean:

```
time for linear1 = 3.56, linear2 = 18.30, extract = 20.06, total = 42.34
### bts finished, everything ok ###
```

`linear1` + `linear2` are SlotToCoeff and CoeffToSlot, `extract` is digit
extraction (the stage the paper changes), and `total` is the whole recryption.
`everything ok` means every slot decrypted to the expected value; a run without
it is not counted. With `HELIB_AUX_ORDER4_EVAL=1` the run also prints, for
example, `HELIB_AUX_ORDER4_EVAL enabled: deg(P)=1223, terms(P)=307, deg(Q)=305`.

## 7. Known deviations between the paper text and the measurements

We list these so that evaluators are not surprised.

1. **Table 5, OFF column, auxiliary radix.** Section 7 says OFF and ON share
   the auxiliary radix `A`. That holds for Tables 6 and 7. In Table 5, OFF runs
   Ma et al.'s driver at its default radix `2B+1 = 27`, and ON runs the
   order-four radix; both build the same chain. Table 5 therefore compares
   against the baseline as published, and its ratio includes the effect of the
   radix as well as the evaluator. `ae/run.sh tab5` reproduces exactly what was
   measured. (With `A` forced on OFF the `m = 2^16` driver does not run.)
2. **Single-trial and first-trial values.** Table 5's last column ("at
   1314–1332") and the power-of-two rows of the degree-scaling table in the
   full version report the first of three trials; the mean of three differs by
   up to 0.05x (`results/power_of_two/`). Set 80-E and the general-ring rows
   of the degree-scaling table are single trials (`results/baseline_matrix/`,
   `results/degree_scaling/`, `results/paper/degscale_gen/`).
3. **Compiler flags.** The paper says `-O2`; HElib's `Release` build, used for
   all runs, compiles with `-O3 -march=native`.
4. **Table 5 runs were concurrent.** The reported Table 5 runs executed OFF and
   ON of a prime at the same time, with two to four runs on the machine
   (`results/paper/tab5_po2/`, timestamps in the logs). Run alone in the Docker
   image, `p = 14401` gives 22.80 s → 12.07 s (1.89x, three trials) against the
   reported 24.39 s → 11.48 s (2.13x); the other primes of the table lie at
   1.84x–1.95x. Expect isolated reruns of that row near 1.9x.

Our reruns in the Docker image, one physical core per run (two for
"combined"), three trials each, logs in `results/ae_rerun_docker/`
(`python3 ae/make_tables.py --logs results/ae_rerun_docker`):

| run | paper | Docker rerun |
|---|---|---|
| Table 5, `p = 14401`, extract OFF → ON | 24.39 → 11.48 s, 2.13x | 22.80 → 12.07 s, 1.89x (item 4) |
| Table 6, set 80B, extract OFF → ON | 100.05 → 57.40 s, 1.74x | 100.86 → 58.59 s, 1.72x |
| Table 7, Ma et al. total | 213.10 s | 211.75 s |
| Table 7, sparse `A = 256` total (not in the paper; shipped log) | 199.40 s | 197.97 s |
| Table 7, order-four total | 131.18 s, 1.62x | 134.83 s, 1.57x |
| Table 7, combined total | 118.49 s, 1.80x | 113.18 s, 1.87x |
| Sec. 7.3, `p = 257, h = 12` | abort | abort, same message |
| Per-attack security table, nine sets | | equal to 0.01 bit |

`results/ae_rerun_docker/tab7_m50731/under_load/` keeps two Table 7 attempts
made while another user loaded the host (load average ≈ 56): single trials of
the order-four row ranged from 81.6 s to 132.0 s of extraction, and `linear2`,
which the evaluator does not touch, slowed by the same factor. Run on an idle
machine.

## 8. Repository map

```
ae/                        build, run and tabulation scripts for evaluation
  build_all.sh             builds build/bin/{fatboot_order4,fatboot_po2,fatboot_ma_largep}
  run.sh                   one target per table (ae/run.sh list)
  make_tables.py           Tables 5-7 from logs, next to the paper's values
  paper_values.json        the numbers printed in the paper
Dockerfile                 reproducible environment
results/paper/             raw logs of the reported runs, one file per run
  tab5_po2/  tab6_general/  tab7_m50731/  admissibility/  degscale_gen/
src/HElib_auxradix_opt/    modified HElib (Apache-2.0, see its NOTICE)
src/BGV-Boot-auxradix-opt/ driver for Tables 6 and 7 (MIT, Ma et al.)
third_party/patches/       our diffs against pinned upstream code
security/                  Lattice Estimator scripts, transcripts and summaries
  built_chain.py           per-attack lambda* at the built chain (Table 2, appendix)
  results/built_chain/     its output for the nine sets
lean/                      Lean 4 + Mathlib formalisation (zero sorry)
scripts/verify_order4.py   plaintext check of the factored form
encrypted_linear_transform/  Section 6 butterfly prototype and logs
okada_comparison/          Section 7.4 composition with Okada et al.
```

`results/paper/` logs are unedited except that the absolute path of the
measurement host's home directory is replaced by `<HOST>`.

## 9. Licence

See `LICENSE`. Our code is MIT. `src/HElib_auxradix_opt` is a modified HElib
(Apache-2.0); `src/BGV-Boot-auxradix-opt` derives from Ma et al. (MIT).
`scripts/install_magma_local.sh` and `scripts/run_polyfunctions_magma_generation.sh`
concern the Geelen et al. baseline, which needs the proprietary Magma system;
they are optional and no table depends on them.
