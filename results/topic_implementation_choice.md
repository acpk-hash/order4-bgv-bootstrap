# Topic Implementation Choice

Date: 2026-04-27

This topic currently has one primary algorithmic route and one secondary
implementation-choice route.

## Primary Route

The primary technical route remains sparse exact cleaning polynomial search for
the non-power auxiliary modulus:

```text
HELIB_EXPLICIT_AUX=256
```

For the target parameter set:

```text
p = 65537
B = 17
A = 256
A^2 = -1 mod p
```

this changes the exact cleaning polynomial from the default 612 nonzero terms
to 307 nonzero terms, and gives the stable extraction-time speedup.

## Implementation Choice Route

The second usable optimization is:

```text
HELIB_AUX_PARALLEL_COEFF2SLOT=1
```

It runs the two independent `coeffToSlot` transforms in the fixed HElib newBTS
pipeline concurrently.

This must be framed as:

```text
pipeline/executor scheduling optimization
implementation choice for fixed newBTS pipeline
wall-clock acceleration
```

It must not be framed as:

```text
new homomorphic linear-transform algorithm
new lower asymptotic complexity
reduced number of automorphisms or key switches
```

## Current Best Experimental Line

Repeat-3 benchmark:

| method | correctness | linear2 | total |
| --- | --- | ---: | ---: |
| aux=35 original | pass | 45.154748s | 210.573619s |
| aux=256 sparse | pass | 43.463309s | 191.931511s |
| aux=256 + parallel coeffToSlot | pass | 23.906354s | 182.353824s |

Combined improvement relative to original aux=35:

```text
total: 210.573619 -> 182.353824 = 13.40% faster
```

Incremental improvement relative to aux=256 sparse:

```text
linear2: 43.463309 -> 23.906354 = 45.00% faster
total:   191.931511 -> 182.353824 =  4.99% faster
```

## Writing Boundary

Use the parallel `coeffToSlot` result as an implementation-selection
optimization that complements the sparse polynomial route.  It is useful for
the experiment-driven story, but the main mathematical novelty should stay on
the null-polynomial / sparse exact polynomial side unless a true fused
linear-transform algorithm is later found.
