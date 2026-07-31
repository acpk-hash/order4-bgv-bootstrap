# Measured results

All figures below are measurements. Nothing is extrapolated except where a
line says so. Raw per run output is in the sibling `.jsonl` files and the
full console logs are in the artifact's `measured_raw/` directory.

Tools: the stock Lattice Estimator at commit `3e48ef4`, and the rotated
primal hybrid of ePrint 2026/279 from `github.com/TabOg/mlwe-hybrids`
(`PrimalHybrid`) with its `lattice_estimator` submodule advanced to the same
commit, so both columns rest on one lattice-estimation codebase. Cost model
`RC.MATZOV`, `sigma = 3.2`, sparse ternary secrets, `m = n` samples.


## 1. Digit extraction at the revised power-of-two chain

Ring `m = 2^16`, main key weight 120, encapsulated weight 26, chain built at
1066 bits for six primes and 1071 for four, three trials, single threaded.
`OFF` is the bounded-support evaluator, `ON` is the order-four evaluator;
the pipeline is otherwise identical and the order-four radix `A` is per prime.

| p | A | slots | OFF (s) | ON (s) | speedup | at the shipped chain |
|---:|---:|---:|---:|---:|---:|---:|
| 1297 | 36 | 8 | 23.00 | 11.78 | 1.95x | 1.63x |
| 1601 | 40 | 32 | 22.93 | 12.39 | 1.85x | 1.69x |
| 2521 | 71 | 4 | 24.71 | 12.82 | 1.93x | 1.62x |
| 3137 | 56 | 32 | 21.98 | 11.43 | 1.92x | 1.64x |
| 4513 | 95 | 16 | 20.39 | 11.11 | 1.84x | 1.65x |
| 7057 | 84 | 8 | 21.13 | 11.25 | 1.88x | 1.60x |
| 13457 | 116 | 8 | 21.35 | 10.94 | 1.95x | 1.74x |
| 14401 | 120 | 32 | 24.39 | 11.48 | 2.13x | 1.67x |
| 15377 | 124 | 8 | 21.66 | 11.74 | 1.85x | 1.71x |
| 65537 | 256 | 32768 | 20.06 | 10.81 | 1.85x | 1.59x |

Range 1.84x to 2.13x, mean 1.91x, against 1.59x to 1.74x at the shipped
1314 to 1332 bit chains. Every prime improves. The baseline is essentially
unchanged by the shorter chain while the order-four evaluator gets faster,
because digit-extraction cost is an operation count fixed by `(p, A, B)`
multiplied by a per-operation cost that the number of residue limbs sets.


## 2. Whole bootstrapping, p = 65537

The least favourable case for a total ratio: this prime has the most slots,
so the linear transforms, which the contribution does not touch, dominate.

| stage | OFF (s) | ON (s) | ratio |
|---|---:|---:|---:|
| bootstrapping key switch | 9.50 | 9.58 | 0.99 |
| coefficients to slots | 9.11 | 9.13 | 1.00 |
| digit extraction | 20.06 | 10.81 | 1.85 |
| slots to coefficients | 3.56 | 3.61 | 0.99 |
| total recryption | 42.34 | 33.25 | 1.27 |

The three stages the method does not touch move by 0.99 to 1.00, which is the
measurement confirming it is confined where we claim. Capacity remaining is
344.5 bits for OFF and 319.4 for ON, so on the throughput basis used by the
prior work, remaining capacity divided by total time, the gain is 1.23x
against 1.27x on wall clock alone.


## 3. Digit extraction on general cyclotomic rings

At the recommended chains, same method. These rings are out of scope of
ePrint 2026/279.

| set | p | built log2 Q | OFF (s) | ON (s) | speedup |
|---|---:|---:|---:|---:|---:|
| 128B | 2917 | 1153.0 | 127.23 | 70.68 | 1.80x |
| 80A | 4513 | 1524.2 | 89.31 | 50.12 | 1.78x |
| 80B | 14401 | 1642.0 | 100.05 | 57.40 | 1.74x |

The two families behave differently and we report them separately rather than
assuming one rule. Trimming the power-of-two chain improves the ratio;
lengthening the general-ring chain leaves it flat, 1.79x at the 1055 bit chain
the paper runs against 1.78x at the recommended 1524.


## 4. Where the security numbers come from

| file | contents |
|---|---|
| `rot_2026_279_twelve.jsonl` | the twelve in-scope power-of-two instances, both MitM heuristics |
| `modulus_sweep.jsonl` | exhaustive modulus search at weight 120, step 5 bits |
| `weight_sweep.jsonl` | main-key weight sweeps, at a fixed chain and at chains re-sized for the heavier key |
| `encapsulated_and_limits.jsonl` | minimum encapsulated weight clearing 80 bits, and the probes showing the weight knob does not reach 128 |
| `stock_estimator.jsonl` | the stock estimator, six attacks per instance |
| `cross_verification.jsonl` | independent recomputation and the integer boundary searches |
| `bootstrap_power_of_two.jsonl` | ciphertext bootstrapping, ten primes, both evaluators |
| `bootstrap_general_rings.jsonl` | ciphertext bootstrapping, general rings |

The correct post-attack figure for an instance is the minimum over every
attack run, not the stock figure minus a delta: the stock minimum often comes
from an attack cheaper than the rotated tool's own plain hybrid.

