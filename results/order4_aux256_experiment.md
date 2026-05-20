# Order-Four Cleaner Evaluator Experiment

## Implementation

The implementation is isolated under `upload` and does not modify the baseline
trees.  The new flag is:

```text
HELIB_AUX_ORDER4_EVAL=1
```

It is used together with:

```text
HELIB_EXPLICIT_AUX=256
```

The code path is in:

```text
src/HElib_auxradix_opt/src/extractDigits.cpp
```

The helper first certifies that the cleaner has the order-four support pattern:

```text
P(X) = c_1 X + X^3 Q(X^4).
```

If the pattern is not present, the code falls back to the original `polyEvalNew`
path.  In the target experiment the flag is triggered:

```text
HELIB_AUX_ORDER4_EVAL enabled: deg(P)=1223, terms(P)=307, deg(Q)=305, terms(Q)=306
```

## Commands

Build:

```bash
cd .
./build_upload.sh
```

Run the repeat-three benchmark:

```bash
cd .
REPEAT=3 ./run_order4_aux256_benchmark.sh
```

Parse against the previous aux-256 sparse-cleaner baseline:

```bash
python3 scripts/parse_helib_bootstrap_timings.py \
  --baseline results/upload_aux_256_repeat3.log \
  --optimized results/order4_aux256_repeat3.log \
  --output results/order4_aux256_vs_aux256_repeat3.json
```

## Repeat-Three Results

All three bootstraps finished with:

```text
### bts finished, everything ok ###
```

Mean timings in seconds:

| Variant | linear1 | linear2 | extract | total |
|---|---:|---:|---:|---:|
| Ma baseline aux=35 (`/baselines/BGV-Boot-for-Large-p`) | 4.428140 | 45.268850 | 160.783267 | 213.102573 |
| aux=256 sparse cleaner | 4.380695 | 44.512004 | 147.909934 | 199.397017 |
| aux=256 + order-four evaluator | 4.272161 | 43.623283 | 80.770389 | 131.180978 |

Improvement over aux=256 sparse cleaner:

| Metric | Improvement |
|---|---:|
| linear1 | 2.48% |
| linear2 | 2.00% |
| extract | 45.39% |
| total | 34.21% |

Improvement over Ma baseline aux=35:

| Metric | Improvement |
|---|---:|
| linear1 | 3.52% |
| linear2 | 3.64% |
| extract | 49.76% |
| total | 38.44% |

## Interpretation

This turns the previous sparse-cleaner result into a stronger two-layer claim:

1. choose the order-four auxiliary radix to make the exact cleaner sparse;
2. evaluate the resulting cleaner through its algebraic form
   `P(X)=c_1 X + X^3 Q(X^4)`.

The second step is algorithmic, not merely parallel scheduling.  It changes the
homomorphic polynomial-evaluation circuit used inside digit extraction while
preserving the surrounding HElib bootstrapping pipeline.

The linear-transform times change only slightly.  The speedup is concentrated
in `extractDigitsThin`, which drops from `147.909934s` to `80.770389s` against
the aux-256 sparse-cleaner baseline.
