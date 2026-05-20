# Artifact Notes

## What Is Included

This GitHub artifact includes:

- a modified HElib tree under `src/HElib_auxradix_opt`;
- the large-prime BGV fatboot driver under `src/BGV-Boot-auxradix-opt`;
- scripts for polynomial search, degree/term certificates, order-adaptive sweeps,
  transform planning, exact plaintext verification, and log parsing;
- selected benchmark logs and parsed JSON files;
- the LNCS paper source and compiled PDF;
- a short Chinese mathematical theory note for presentations.

The large exploratory external repositories from the working directory are not
included. In particular, the original `upload/external` directory contained
large Fheanor and Lattigo checkouts and was intentionally omitted.

## Main Flags

Sparse cleaner:

```text
HELIB_EXPLICIT_AUX=256
```

Order-four evaluator:

```text
HELIB_AUX_ORDER4_EVAL=1
```

HElib polynomial cache location:

```text
HELIB_ZZX_CACHE_DIR=cache/saved_ZZX
```

The run scripts set `HELIB_ZZX_CACHE_DIR` automatically.

## Main Result Files

```text
results/order4_aux256_experiment.md
results/order4_aux256_repeat3.log
results/order4_aux256_vs_aux256_repeat3.json
results/baseline_BGV_Boot_for_Large_p_i4_repeat3.log
results/baselines_BGV_Boot_for_Large_p_vs_order4_aux256_repeat3.json
results/baselines_BGV_Boot_for_Large_p_vs_aux256_sparse_repeat3.json
results/upload_aux_256_repeat3.log
results/upload_aux_default_repeat3.log
results/order_adaptive_cleaner_sweep_p65537.json
results/linear_transform_factor_planner_all.txt
results/plaintext_rader97_seed7.txt
```

## Rebuilding

```bash
./build_upload.sh
```

This installs the modified HElib only under:

```text
local/helib_auxradix_opt
```

and builds the fatboot executable at:

```text
src/BGV-Boot-auxradix-opt/build/fatboot
```

## Running

Sparse-cleaner baseline:

```bash
REPEAT=3 ./run_auxradix_benchmark.sh
```

Order-four evaluator:

```bash
REPEAT=3 ./run_order4_aux256_benchmark.sh
```

The direct Ma baseline comparison object is the checked-out repository under
`/home/luck/xzy/0424project/baselines/BGV-Boot-for-Large-p`. The repeat-three
baseline log in this artifact was produced with:

```bash
cd /home/luck/xzy/0424project/baselines/BGV-Boot-for-Large-p
timeout 1800s ./build/fatboot i=4 h=12 t=-1 newbts=1 repeat=3 newks=1 \
  > /home/luck/xzy/0424project/github_order4_cleaner/results/baseline_BGV_Boot_for_Large_p_i4_repeat3.log 2>&1
```

## Parsing

```bash
python3 scripts/parse_helib_bootstrap_timings.py \
  --baseline results/upload_aux_256_repeat3.log \
  --optimized results/order4_aux256_repeat3.log
```

Expected output:

```text
metric baseline_mean optimized_mean improvement_percent
linear1 4.380695 4.272161 2.48
linear2 44.512004 43.623283 2.00
extract 147.909934 80.770389 45.39
total 199.397017 131.180978 34.21
```

To compare against `/baselines/BGV-Boot-for-Large-p` directly:

```bash
python3 scripts/parse_helib_bootstrap_timings.py \
  --baseline results/baseline_BGV_Boot_for_Large_p_i4_repeat3.log \
  --optimized results/order4_aux256_repeat3.log
```

Expected output:

```text
metric baseline_mean optimized_mean improvement_percent
linear1 4.428140 4.272161 3.52
linear2 45.268850 43.623283 3.64
extract 160.783267 80.770389 49.76
total 213.102574 131.180978 38.44
```

## Order-Adaptive Cleaner Sweep

```bash
python3 scripts/order_adaptive_cleaner_sweep.py --p 65537 \
  --json --output results/order_adaptive_cleaner_sweep_p65537.json
```

This is offline planner evidence only. It supports the paper's conservative
threshold story:

```text
order 4  -> stable at B=17
order 8  -> first support-safe at B<=7
order 16 -> tiny-local-support only (B<=1)
```

## Structured Linear-Transform Evidence

Planner:

```bash
python3 scripts/linear_transform_factor_planner.py --all \
  > results/linear_transform_factor_planner_all.txt
```

Exact plaintext Rader verification:

```bash
python3 scripts/plaintext_rader_evalmap.py --ell 97 --seed 7 \
  > results/plaintext_rader97_seed7.txt
```

These support the paper's transform-design section, not a positive ciphertext
speedup claim.

## Paper

The paper source is in:

```text
paper/main.tex
```

Build it with:

```bash
cd paper
make
```

The included `paper/main.pdf` was checked with `pdfinfo` and `pdffonts`; rerun
those checks after regeneration.

## Theory Note

The presentation-oriented mathematical explanation is in:

```text
docs/order4_cleaner_theory.tex
docs/order4_cleaner_theory.pdf
```

It explains the fixed-support cleaner viewpoint, the order-four auxiliary radix
identity, the sparse exponent filter, and the resulting evaluator
decomposition.

It also states the general parameter theorem used for claim boundaries:

```text
exact-cleaner correctness:
  p prime, p > 4B(B+1), A = 2B+1

order-four sparse cleaner:
  p prime, p = 1 mod 4, p > 8B^2, A^2 = -1 mod p
```

The complete coefficient lists for the concrete `p=65537, B=17` cleaners are
in:

```text
docs/cleaner_polynomial_forms_p65537_B17.json
```

They can be regenerated with:

```bash
python3 scripts/export_cleaner_polynomial_forms.py \
  --json-output docs/cleaner_polynomial_forms_p65537_B17.json \
  --tex-output docs/cleaner_polynomial_forms_p65537_B17.tex
```

## Claim Boundary

This artifact does not claim a globally optimal FHE bootstrapping method.

The positive implemented claim is specific:

```text
fixed HElib large-prime BGV thin-bootstrap pipeline
+ bounded-support exact cleaner
+ order-four auxiliary radix A=256
+ specialized evaluator for P(X)=c_1 X + X^3 Q(X^4)
```

The broader paper framework additionally uses:

```text
- higher-order cleaner feasibility / threshold results from exact offline sweeps
- Galois-structured linear-transform design evidence from planner output and
  exact plaintext verification
```

Several linear-transform ideas were also tested. They are included as negative
or engineering evidence in `results/linear_transform_attempts_summary.md`.

In particular, the Rader/CRT route should currently be read as:

```text
promising algebraic design route,
not yet a demonstrated ciphertext speedup.
```
