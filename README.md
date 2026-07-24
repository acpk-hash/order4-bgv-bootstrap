# Character-Filtered Exact Cleaners for Large-Prime BGV Bootstrapping

This repository is a GitHub-ready artifact for a conservative but broader paper
story built around three layers:

1. a **positive implemented ciphertext result** for the order-four exact cleaner
   and its specialized evaluator in large-prime BGV bootstrapping;
2. an **order-adaptive cleaner framework** supported by exact offline search and
   threshold sweeps;
3. a **Galois-structured linear-transform design route** supported by planner
   output, exact plaintext verification, and explicit negative/guardrail
   ciphertext evidence.

The artifact is intentionally self-contained but lightweight: it keeps the
modified HElib source, the large-prime fatboot driver, scripts, logs, and paper
sources, while excluding build directories, local installs, caches, and large
external exploration repositories.

## Main Positive Result

For the HElib thin-bootstrapping parameter set

```text
p = 65537, B = 17, m = 50731, mvec = [97, 523], ords = [96, 29]
```

the auxiliary radix `A=256` satisfies `A^2 = -1 mod p`. The exact cleaner has
the order-four form

```text
P(X) = c_1 X + X^3 Q(X^4).
```

The artifact implements two layers:

1. `HELIB_EXPLICIT_AUX=256`: choose the sparse exact cleaner.
2. `HELIB_AUX_ORDER4_EVAL=1`: evaluate the cleaner through the order-four
   decomposition above.

Repeat-three real ciphertext results. The first row is a direct run of
`/home/user/experiments/baselines/BGV-Boot-for-Large-p`; the other rows are
the isolated implementation in this artifact.

| Variant | linear1 | linear2 | extract | total |
|---|---:|---:|---:|---:|
| Ma baseline aux=35 | 4.428140s | 45.268850s | 160.783267s | 213.102573s |
| sparse aux=256 | 4.380695s | 44.512004s | 147.909934s | 199.397017s |
| aux=256 + order-four evaluator | 4.272161s | 43.623283s | 80.770389s | 131.180978s |

Speedup over sparse aux=256:

```text
extract = 45.39%
total   = 34.21%
```

Speedup over Ma baseline aux=35:

```text
extract = 49.76%
total   = 38.44%
```

## Broader Conservative Paper Scope

The broader paper claim is **not** that all of the following are implemented in
ciphertext form. Instead:

- the order-four cleaner/evaluator is the **implemented positive path**;
- higher-order subgroup cleaners are supported by **exact offline planner
  evidence**;
- the normal-basis / factorized transform route is supported by **planner output
  and exact plaintext verification**, while current ciphertext integrations are
  deliberately reported as negative or engineering-only evidence.

This avoids overclaiming while still documenting the larger Galois-structured
framework.

## Repository Layout

```text
src/HElib_auxradix_opt/        modified HElib source
src/BGV-Boot-auxradix-opt/     large-prime BGV fatboot driver
scripts/                       cleaner search, planner, and transform scripts
results/                       selected benchmark logs and summaries
lean/                          machine-checked proofs (Lean 4, zero sorry)
okada_comparison/              AC'23 Galois-evaluator reproduction + composition benchmark
build_upload.sh                isolated build script
run_auxradix_benchmark.sh      sparse-cleaner baseline benchmark
run_order4_aux256_benchmark.sh order-four evaluator benchmark
summarize_logs.py              compact timing summarizer
ARTIFACT.md                    detailed reproducibility notes
```

## Build

Dependencies are the same as HElib's normal CMake build, including a C++17
compiler, CMake, GMP, and NTL.

```bash
cd /path/to/order-four-cleaner
./build_upload.sh
```

The build script writes only inside this repository:

```text
build/
local/
cache/
```

These paths are ignored by Git.

## Reproduce The Main Benchmark

Run a quick smoke test with one bootstrap:

```bash
REPEAT=1 ./run_order4_aux256_benchmark.sh
```

Run the repeat-three benchmark used in the paper:

```bash
REPEAT=3 ./run_order4_aux256_benchmark.sh
```

Compare against the included sparse aux-256 baseline log:

```bash
python3 scripts/parse_helib_bootstrap_timings.py \
  --baseline results/upload_aux_256_repeat3.log \
  --optimized results/order4_aux256_repeat3.log
```

Compare against the directly measured `/baselines/BGV-Boot-for-Large-p`
baseline log:

```bash
python3 scripts/parse_helib_bootstrap_timings.py \
  --baseline results/baseline_BGV_Boot_for_Large_p_i4_repeat3.log \
  --optimized results/order4_aux256_repeat3.log
```

Expected correctness marker:

```text
### bts finished, everything ok ###
```

Expected order-four trigger:

```text
HELIB_AUX_ORDER4_EVAL enabled: deg(P)=1223, terms(P)=307, deg(Q)=305, terms(Q)=306
```

## Reproduce The New Offline Cleaner/Transform Evidence

### Order-adaptive cleaner sweep

```bash
python3 scripts/order_adaptive_cleaner_sweep.py --p 65537 \
  --json --output results/order_adaptive_cleaner_sweep_p65537.json
```

This produces the conservative threshold story used in the paper:

- order four is the stable global choice at `B=17`;
- order eight becomes support-safe only once `B <= 7`;
- order sixteen is restricted to tiny local supports.

### Character-projected search on the implemented global case

```bash
python3 scripts/character_projected_cleaner_search.py \
  --p 65537 --B 17 --max-order 64 --max-decomposition 16 \
  --aux 35 256 65281 \
  --output results/character_projected_p65537_B17.json --json
```

For the target parameter set, the best valid candidate remains `A=256`. Its
best decomposition is modulo `4`, with estimated evaluator multiplication cost
`41` versus `73` for generic Paterson--Stockmeyer on the degree-1223 cleaner.
The generic radix `A=35` is correct but its best decomposition is only the
odd-polynomial modulo-`2` structure with estimated cost `55`.

### Structured linear-transform planner and exact plaintext check

```bash
python3 scripts/linear_transform_factor_planner.py --all \
  > results/linear_transform_factor_planner_all.txt

python3 scripts/plaintext_rader_evalmap.py --ell 97 --seed 7 \
  > results/plaintext_rader97_seed7.txt
```

The target set E planner output highlights:

```text
D=96 | HElib bsgs auts~18 | radix [2,2,2,2,2,3] auts~7 | primitive-ell Rader auts~9
D=29 | HElib bsgs auts~9  | radix [29] auts~28
```

and the exact plaintext Rader verification reports:

```text
verified=True
operation_ratio=0.7031
conv_length=96 conv_radices=[2, 2, 2, 2, 2, 3]
```

These are design/planner claims, not a demonstrated ciphertext speedup claim.

## Key Files

- Implementation:
  `src/HElib_auxradix_opt/src/extractDigits.cpp`
- Main experiment note:
  `results/order4_aux256_experiment.md`
- Main repeat-three log:
  `results/order4_aux256_repeat3.log`
- Direct Ma baseline repeat-three log:
  `results/baseline_BGV_Boot_for_Large_p_i4_repeat3.log`
- Parsed comparison JSON:
  `results/order4_aux256_vs_aux256_repeat3.json`
  and `results/baselines_BGV_Boot_for_Large_p_vs_order4_aux256_repeat3.json`
- New order-adaptive sweep:
  `results/order_adaptive_cleaner_sweep_p65537.json`
- Character-projected search:
  `scripts/character_projected_cleaner_search.py`
  and `results/character_projected_p65537_B17.json`
- Structured transform planner:
  `results/linear_transform_factor_planner_all.txt`
- Exact plaintext Rader verification:
  `results/plaintext_rader97_seed7.txt`
- Concrete cleaner polynomial forms (regenerable):
  `python3 scripts/export_cleaner_polynomial_forms.py`
  generated by `scripts/export_cleaner_polynomial_forms.py`

## Parameter Boundary

For exact cleaner correctness alone, `p` can be any sufficiently large prime.
For a support bound `B`, choosing

```text
A = 2B + 1
```

is collision-free whenever

```text
p > 4B(B+1).
```

This guarantees a unique degree-`< (2B+1)^2` exact cleaner on
`S_A = {hA+lo : |h|, |lo| <= B}`.

The order-four sparse optimization is stronger and needs more structure:

```text
p prime, p = 1 mod 4, p > 8B^2, choose A with A^2 = -1 mod p.
```

Then the support is automatically collision-free and the canonical cleaner has

```text
P_A(X) = (1/2) X + X^3 Q(X^4).
```

For the broader order-adaptive theory, the paper uses the signed-circulant
sufficient bound

```text
(2B+1) * sum_{j=0}^{d-1} |A|^j < p,  where d = order/2.
```

For the target `p=65537`, this yields the practical thresholds used in the
paper:

```text
order 4  -> B <= 127
order 8  -> B <= 7
order 16 -> B <= 1
```

The positive ciphertext implementation remains the order-four row at `B=17`.

## Claim Boundary

This artifact does not claim a globally optimal FHE bootstrapping method. The
positive implemented claim is specific:

```text
fixed HElib large-prime BGV thin-bootstrap pipeline
+ bounded-support exact cleaner
+ order-four auxiliary radix A=256
+ specialized evaluator for P(X)=c_1 X + X^3 Q(X^4)
```

The broader paper framework additionally claims:

```text
- higher-order cleaner feasibility / threshold results from exact offline sweeps
- Galois-structured linear-transform design evidence from planner output
  and exact plaintext verification
```

Several linear-transform ideas were also tested. They are included as negative
or engineering evidence in `results/linear_transform_attempts_summary.md`.

In particular, the Rader/CRT route should currently be read as:

```text
promising algebraic design route,
not yet a demonstrated ciphertext speedup.
```


## Update: cross-ring validation, moderate primes, and machine-checked proofs

- **`lean/`** — Lean 4 + Mathlib formalization (zero `sorry`): the Order-Four
  Character Filter (coefficient and polynomial forms), the factored form
  `P_A = X/2 + X^3 Q(X^4)`, support characterization, and the Section-4
  framework (isotypic decomposition, monomial support constraint).
  Reproduce: `lake exe cache get && lake build`. See `lean/README.md` for the
  paper-statement ↔ Lean-name interface table.
- **`results/power_of_two/`** — encrypted-domain digit-extraction logs on the
  **power-of-two ring** `m = 2^16 = 65536` for `p = 65537` (32768 slots,
  79-bit security), obtained by porting the order-four evaluator into the
  homomorphic-NTT bootstrapping implementation (`patches/power_of_two_order4_port.patch`):
  generic PS 22.2 s -> order-4 14.0 s (**1.59x**, three-run means); vs. a Ma-style
  baseline radix (`aux=35`, 25.5 s) **1.83x**; thick bootstrapping 34.1 s -> 14.4 s
  (**2.37x**). See `results/power_of_two/suite/` for the definitive sequential
  three-run logs (earlier single-run logs in this directory were taken under
  machine load and overstate absolute times).
- **`results/moderate_primes/`** — encrypted-domain logs for three further
  general-cyclotomic moderate primes (means over three sequential runs):
  `p=1601` (44.2 -> 23.2 s, 1.90x), `p=2917` (28.3 -> 15.4 s, 1.84x),
  `p=8101` (76.3 -> 21.4 s, 3.57x). The order-four evaluator triggers in every
  run (`HELIB_AUX_ORDER4_EVAL enabled: deg(P)=... deg(Q)=...`) and all slots
  decrypt correctly.
- **`encrypted_linear_transform/`** — first **encrypted** instantiation of the
  butterfly (Rader + mixed-radix [6,4,4]) CoeffToSlot block on the full ring
  `m=50731`: correct on all 2784 slots (0 mismatches vs. an independent NTL
  matvec), 41 automorphisms over 17 native rotation classes, 194 bits of noise
  capacity consumed (no overflow even at bits=600); no extended key-switching
  matrices required. At D=96 the monolithic BSGS remains faster (14.9 s vs
  59.1 s), locating the O(log D) crossover at larger D.

## Update: extended cross-parameter evidence (2026-07-23)

- `results/general_rings_new/` — five additional general-ring primes (p = 2521, 4513, 13457, 14401, 15377; sets J–N), OFF/ON with repeat=3: speedups 1.69–1.79x. Note p=2521 runs on m=50731, the same ring as the main p=65537 set.
- `results/power_of_two/expand/` — the power-of-two prime sweep now covers ten primes 1297 <= p <= 65537 on m=2^16 (h=26, degree 727): 1.59–1.74x, same-session sequential OFF/ON pairs.
- `results/degree_scaling/` — key-weight scaling on both ring types (general m=50731: h in {12,18,26}; power-of-two m=2^16: h in {12,18,26,40}); the speedup grows monotonically with the digit-polynomial degree.
- `okada_comparison/` — reproduction of Okada–Player–Pohmann (AC'23, ePrint 2023/1304) from their open-source SEAL implementation (`repro_p257.log`: full p=257 bootstrap, digit extraction 17 key switches, matching their Table 3), plus an encrypted composition benchmark (`compose_bench.cpp`, `t1_FINAL.log`): evaluating the factored Q (deg 55) takes 22 key switches / ~9.0s versus 41 / ~17.2s for the unfactored P_A, bit-identical results.
- `patches/general_ring_sets_J_to_N.patch` — the five new fatboot parameter sets (gens/ords validated via HElib).

## Update: rebuttal-phase experiments (2026-07-24)

- `results/p257_order4/` — boundary demonstration at p=257 inside the power-of-two
  pipeline: lowering the key weight to h=7 makes the support condition hold
  ((2*7+1)^2 = 225 < 257), and the order-four evaluator then reduces encrypted
  digit extraction from 12.6 s to 8.0 s (**1.58x**; deg(P_A)=223 -> deg(Q)=55,
  57 nonzero terms, all 128 slots correct). Reproduce: `./run_p257_boundary.sh`
  (details in `results/p257_order4/SUMMARY.md`).
- `encrypted_linear_transform/crossover/results_trend/` — encrypted
  butterfly-vs-BSGS crossover across block dimensions D = 192, 288, 576, 1152 on
  four purpose-built rings (ring recipe and full table in `TREND_SUMMARY.md`):
  within-row time ratios 1.25 / 1.03 / 0.61 / 0.59, so BSGS wins below and the
  butterfly wins above a crossover just past D=288. Stage chains are
  plaintext-verified (`derive_trend.py`, `stages{192,288,576,1152}.txt`) and both
  encrypted paths decrypt with 0 slot mismatches.
  Reproduce: `./run_crossover_trend.sh`.
- `run_po2_prime_sweep.sh` — self-contained runner for the ten-prime
  power-of-two sweep of `results/power_of_two/expand/` (h=26; expected
  1.59-1.74x), with the prime -> parameter-index -> radix table inline.
