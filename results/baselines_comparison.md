# Baseline Comparison Against `/home/user/experiments/baselines`

Date: 2026-04-29

The main directly comparable baseline is
`/home/user/experiments/baselines/BGV-Boot-for-Large-p`, the Ma large-prime
BGV bootstrapping artifact.  The run uses the same parameter set as the order-four
cleaner experiment:

```text
p = 65537, r = 1, bits = 1500, m = 50731, mvec = [97, 523],
gens = [48117, 5239], ords = [96, 29], h = 12, t = -1, repeat = 3.
```

Command used for the direct baseline:

```bash
cd /home/user/experiments/baselines/BGV-Boot-for-Large-p
timeout 1800s ./build/fatboot i=4 h=12 t=-1 newbts=1 repeat=3 newks=1 \
  > /home/user/order4-artifact/results/baseline_BGV_Boot_for_Large_p_i4_repeat3.log 2>&1
```

All compared logs finish with:

```text
### bts finished, everything ok ###
```

## Timings

| variant | source log | linear1 | linear2 | extract | total |
| --- | --- | ---: | ---: | ---: | ---: |
| Ma baseline, aux=35 | `baseline_BGV_Boot_for_Large_p_i4_repeat3.log` | 4.428140s | 45.268850s | 160.783267s | 213.102573s |
| sparse cleaner, aux=256 | `upload_aux_256_repeat3.log` | 4.380695s | 44.512004s | 147.909934s | 199.397017s |
| order-four evaluator, aux=256 | `order4_aux256_repeat3.log` | 4.272161s | 43.623283s | 80.770389s | 131.180978s |

## Improvement Against `/baselines/BGV-Boot-for-Large-p`

| optimized variant | extract speedup | total speedup |
| --- | ---: | ---: |
| sparse cleaner, aux=256 | 8.01% | 6.43% |
| order-four evaluator, aux=256 | 49.76% | 38.44% |

The polynomial/cost part is independent of timing noise:

| cleaner | degree | nonzero terms | structure |
| --- | ---: | ---: | --- |
| Ma aux=35 | 1223 | 612 | odd powers |
| our aux=256 | 1223 | 307 | `X` and `X^(4j+3)` |

Thus the degree is not lower.  The measured win comes from a sparser exact
cleaner and the order-four decomposition

```text
P_256(X) = 32769 X + X^3 Q(X^4),
```

which changes a generic degree-1223 evaluation into a lower-degree evaluation
of `Q` at `X^4`.

## Other Baselines In The Folder

`HElib` is the canonical upstream HElib code path.  It implements global
digit-extraction polynomials such as `buildDigitPolynomial` and Chen-Han
`compute_magic_poly`, but it is not the same bounded-support large-prime
cleaner benchmark.

`Bootstrapping_Polyfunctions` contains the GIKV/polyfunction and null-lattice
route.  Its cached `poly*.txt` examples are for small `p,e` settings and do not
instantiate the same `p=65537, B=17` support.

`bgv-bootstrapping-with-homomorphic-NTT` targets power-of-two cyclotomics and
has a cached `saved_ZZX/12289_3_ZZX.txt`, which is a different parameter set.

`lattigo`, `openfhe-development`, and `fheanor` are useful modern-library
context, but they are not drop-in comparisons for this HElib large-prime BGV
bounded-support cleaner experiment.
