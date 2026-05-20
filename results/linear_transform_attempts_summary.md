# Linear Transform Combination Attempts

Date: 2026-04-27

Baseline used for comparison:

```text
HELIB_EXPLICIT_AUX=256, i=4, h=12, t=-1, newbts=1, newks=1, thick=0, repeat=1
```

The current stable positive result remains the aux=256 sparse extraction
polynomial path.  After adding the guarded experimental flags, the default
path was re-run and still passes:

| run | correctness | linear1 | linear2 | extract | total |
| --- | --- | ---: | ---: | ---: | ---: |
| `upload_aux_256_after_flags_repeat1.log` | pass | 4.560525s | 44.365585s | 145.713271s | 196.938535s |

Compared with the earlier upload baseline:

| run | correctness | linear1 | linear2 | extract | total |
| --- | --- | ---: | ---: | ---: | ---: |
| `upload_aux_256_repeat1.log` | pass | 4.458518s | 43.957204s | 144.869405s | 195.591277s |

The difference is normal run-to-run variation; the guarded code is off by
default and does not change the working pipeline.

## Attempted Linear-Transform Variants

| idea | flag / option | correctness | linear1 | linear2 | extract | total | conclusion |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| force all dimensions into BSGS | `forcebsgs=1` | pass | 4.875734s | 48.179984s | 145.014466s | 200.366732s | correct but slower |
| stage-batched duplicate coeffToSlot | `HELIB_AUX_BATCHED_COEFF2SLOT=1` | pass | 4.559934s | 45.213966s | 151.478779s | 203.587346s | correct but slower |
| parallel duplicate coeffToSlot | `HELIB_AUX_PARALLEL_COEFF2SLOT=1` | pass | 4.677730s | 24.713672s | 153.675630s | 185.448707s | correct and faster |
| subtract unclean overflow before coeffToSlot | `HELIB_AUX_COMBINED_PRE_SUBTRACT=1` | fail | 4.489666s | 22.327237s | 0.000000s | 29.101530s | leaves uncontrolled residue |
| evaluate sparse extraction before coeffToSlot | `HELIB_AUX_PRE_EXTRACT_SUBTRACT=1` | fail | 4.395869s | 22.097041s | 145.643547s | 174.373047s | extraction does not commute with thin packing map |
| factor D=96 step with Rader/CRT/NTT stages | `HELIB_THINSTEP2_RADER97=1` | fail before timing | n/a | n/a | n/a | n/a | unfused stages break the noise/setup budget |

## Positive Parallel Batch Result

The useful second-path implementation is not an algebraic reordering.  It keeps
the original sparse aux=256 extraction path and runs the two independent
`coeffToSlot` transforms in the newBTS path concurrently:

```text
HELIB_AUX_PARALLEL_COEFF2SLOT=1
```

Repeat-3 command:

```bash
cd .
HELIB_AUX_PARALLEL_COEFF2SLOT=1 HELIB_EXPLICIT_AUX=256 \
  timeout 1500s ./src/BGV-Boot-auxradix-opt/build/fatboot \
  i=4 h=12 t=-1 newbts=1 newks=1 thick=0 repeat=3 \
  > results/parallel_coeff2slot_aux256_repeat3.log 2>&1
```

Repeat-3 result:

| run | correctness | linear1 | linear2 | extract | total |
| --- | --- | ---: | ---: | ---: | ---: |
| aux=35 original repeat-3 | pass | 4.377296s | 45.154748s | 158.456112s | 210.573619s |
| aux=256 sparse repeat-3 | pass | 4.205747s | 43.463309s | 141.778187s | 191.931511s |
| aux=256 + parallel coeffToSlot repeat-3 | pass | 4.603212s | 23.906354s | 151.120953s | 182.353824s |

Relative to aux=256 sparse repeat-3:

```text
linear2: 43.463309 -> 23.906354 = 45.00% faster
total:   191.931511 -> 182.353824 =  4.99% faster
```

Relative to the original aux=35 repeat-3:

```text
total: 210.573619 -> 182.353824 = 13.40% faster
```

## Current Conclusion

The two direct ways of combining the sparse extraction polynomial with the
second coeffToSlot are not valid for full-modulus random BGV plaintexts:

```text
1. Skipping extraction is wrong because Enc(I) is only known modulo aux.
2. Moving extraction before coeffToSlot is wrong because the cleaning
   polynomial is valid after thin packing into slots, not in the original
   coefficient representation.
```

A credible second optimization claim therefore should be stated narrowly:

```text
The current positive linear-transform path is parallel batching of the two
duplicate coeffToSlot transforms in the fixed HElib newBTS pipeline.  It is a
wall-clock executor optimization, not a new algebraic reduction in the number
of automorphisms.  The direct algebraic reorderings should be treated as
negative experimental evidence.
```
