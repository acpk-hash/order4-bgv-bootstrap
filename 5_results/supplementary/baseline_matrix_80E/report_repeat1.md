# Baseline Parameter Matrix Report

This report is generated from local logs.  Timings are real ciphertext
fatboot runs when the corresponding status is `ok` or `cached` with a
passing correctness marker.

## Ma-large-p

| set | p | h | B | type-B? | generic poly | best/chosen poly | variants |
|---|---:|---:|---:|---|---|---|---|
| A | 17 | 14 | 18 | False | A=37, deg=None, terms=None | A=None, deg=None, terms=None, d=None, est=None | ma_typeA:cached, ma_typeB:skipped, ours_sparse:skipped, ours_order4:skipped |
| B | 127 | 12 | 17 | False | A=35, deg=None, terms=None | A=None, deg=None, terms=None, d=None, est=None | ma_typeA:cached, ma_typeB:skipped, ours_sparse:skipped, ours_order4:skipped |
| C | 257 | 12 | 17 | False | A=35, deg=None, terms=None | A=None, deg=None, terms=None, d=None, est=None | ma_typeA:cached, ma_typeB:skipped, ours_sparse:skipped, ours_order4:skipped |
| D | 8191 | 12 | 17 | True | A=35, deg=1223, terms=612 | A=45, deg=1223, terms=611, d=2, est=55 | ma_typeA:cached, ma_typeB:cached, ours_sparse:cached, ours_order4:skipped |
| E | 65537 | 12 | 17 | True | A=35, deg=1223, terms=612 | A=256, deg=1223, terms=307, d=4, est=41 | ma_typeA:cached, ma_typeB:cached, ours_sparse:cached, ours_order4:cached |

### Ciphertext Timings

| set | variant | status | pass | linear1 | linear2 | extract | total | extract speedup vs Ma-typeB | total speedup vs Ma-typeB |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| A | ma_typeA | cached | True | 2.186 | 11.192 | 28.599 | 43.298 | -- | -- |
| A | ma_typeB | skipped | None | -- | -- | -- | -- | -- | -- |
| A | ours_sparse | skipped | None | -- | -- | -- | -- | -- | -- |
| A | ours_order4 | skipped | None | -- | -- | -- | -- | -- | -- |
| B | ma_typeA | cached | True | 4.180 | 30.041 | 52.408 | 89.127 | -- | -- |
| B | ma_typeB | skipped | None | -- | -- | -- | -- | -- | -- |
| B | ours_sparse | skipped | None | -- | -- | -- | -- | -- | -- |
| B | ours_order4 | skipped | None | -- | -- | -- | -- | -- | -- |
| C | ma_typeA | cached | True | 5.646 | 31.231 | 48.140 | 87.408 | -- | -- |
| C | ma_typeB | skipped | None | -- | -- | -- | -- | -- | -- |
| C | ours_sparse | skipped | None | -- | -- | -- | -- | -- | -- |
| C | ours_order4 | skipped | None | -- | -- | -- | -- | -- | -- |
| D | ma_typeA | cached | True | 3.624 | 19.228 | 28.752 | 53.348 | 77.75 | 68.71 |
| D | ma_typeB | cached | True | 3.525 | 35.819 | 129.200 | 170.470 | 0.00 | 0.00 |
| D | ours_sparse | cached | True | 3.564 | 35.685 | 126.193 | 167.317 | 2.33 | 1.85 |
| D | ours_order4 | skipped | None | -- | -- | -- | -- | -- | -- |
| E | ma_typeA | cached | True | 4.468 | 22.205 | 33.951 | 62.706 | 78.31 | 69.77 |
| E | ma_typeB | cached | True | 4.475 | 44.160 | 156.520 | 207.442 | 0.00 | 0.00 |
| E | ours_sparse | cached | True | 4.528 | 43.702 | 143.937 | 194.446 | 8.04 | 6.26 |
| E | ours_order4 | cached | True | 4.512 | 44.068 | 81.391 | 132.236 | 48.00 | 36.25 |

## Ma-NTT

| set | p | h | B | type-B? | generic poly | best/chosen poly | variants |
|---|---:|---:|---:|---|---|---|---|
| I | 65537 | 26 | 13 | True | A=27, deg=727, terms=364 | A=256, deg=727, terms=183, d=4, est=33 | ntt_article_optimized:ok, ntt_article_baseline:timeout, ours_ciphertext:not_integrated |
| II | 8191 | 24 | 12 | True | A=25, deg=623, terms=312 | A=281, deg=623, terms=311, d=2, est=39 | ntt_article_optimized:ok, ntt_article_baseline:ok, ours_ciphertext:not_integrated |
| III | 131071 | 26 | 13 | True | A=27, deg=727, terms=364 | A=27, deg=727, terms=364, d=2, est=42 | ntt_article_optimized:ok, ntt_article_baseline:ok, ours_ciphertext:not_integrated |

### Ciphertext Timings

| set | variant | status | pass | linear1 | linear2 | extract | total | extract speedup vs Ma-typeB | total speedup vs Ma-typeB |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| I | ntt_article_optimized | ok | True | 3.857 | 10.855 | 5.472 | 20.564 | -- | -- |
| I | ntt_article_baseline | timeout | False | -- | -- | -- | -- | -- | -- |
| I | ours_ciphertext | not_integrated | None | -- | -- | -- | -- | -- | -- |
| II | ntt_article_optimized | ok | True | 2.992 | 12.278 | 5.872 | 21.567 | -- | -- |
| II | ntt_article_baseline | ok | True | 11.734 | 42.759 | 5.670 | 60.543 | -- | -- |
| II | ours_ciphertext | not_integrated | None | -- | -- | -- | -- | -- | -- |
| III | ntt_article_optimized | ok | True | 3.962 | 13.010 | 5.220 | 22.584 | -- | -- |
| III | ntt_article_baseline | ok | True | 34.252 | 137.575 | 5.924 | 178.141 | -- | -- |
| III | ours_ciphertext | not_integrated | None | -- | -- | -- | -- | -- | -- |

## Non-Comparable Local Baselines

- `Bootstrapping_Polyfunctions`: requires Magma; no `magma` binary is available in this environment, so polynomial regeneration cannot be run locally.
- `fheanor/openfhe-development/lattigo`: different scheme/library benchmarks; no drop-in HElib BGV digit-extraction pipeline for same ciphertext timing.

