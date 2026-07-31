# Concrete security of the reported parameter sets

This directory answers the question "what are the concrete bit-security levels
of these parameter sets, and what does ePrint 2026/279 do to them", with the raw
per run output kept alongside every summary.

## Why two tools

The stock Lattice Estimator applies to every ring we use: it consumes only the
LWE parameters and is indifferent to the ring.

The rotated primal hybrid of ePrint 2026/279 does not have that generality. It
is developed for the binomial ring `Z_q[X]/(X^N + 1)`, and the lemma supplying
the coefficient isometries together with the algorithms consuming it are stated
with `N` a power of two. The test for whether one of our rings is in that class
is exact: `Phi_m(X) = X^{phi(m)} + 1` precisely when `m` is a power of two,
since `Phi_{2^k}(X) = X^{2^{k-1}} + 1`, and `X^N - 1` is a cyclotomic polynomial
only for `m = 1`. So the question reduces to whether the cyclotomic index is a
power of two, and it admits no borderline case.

Our sets therefore split in two: general cyclotomic rings, estimated with the
stock tool alone, and power-of-two rings, where both apply and both are reported.

## The model: two LWE instances

BGV bootstrapping with sparse-secret key encapsulation (Ma et al., EUROCRYPT
2024, Sec. 4.3, following Bossuat et al.) defines **two** LWE instances, and
the security of the scheme is the smaller of the two:

| instance | dimension | modulus | secret weight |
|---|---|---|---|
| (i) main key | `phi(m)` | full chain `Q` | `h' = 120` (HElib default) |
| (ii) encapsulated bootstrapping key | `phi(m)` | `q_boot = q0 * R` | `h` (12, 24, 26, ...) |

The `h` reported in the experiment tables is the **encapsulated** weight, which
fixes the noise bound `B` and therefore the digit-extraction degree. It is not
the main secret key. Every run in `results/` prints both: `hwt = 120` is the
main key and `encapHwt = ...` is the encapsulated key. `log2(q_boot)` is the sum
of the printed `final log2(qks)` and `final log2(R)`.

Note that the two baseline papers swap these symbols. ePrint 2024/115 writes the
main key weight as `h'`, which is the convention above; ePrint 2024/164 writes
it as `h`. Any comparison against a published table has to fix one convention
first.

## The post-attack figure is a minimum, not a subtraction

For an instance in scope we run five attacks and take the minimum: the stock
estimator's own minimum, and the plain and rotated hybrids of the 2026/279 tool
each with and without meet in the middle.

Taking the stock figure and subtracting a delta is wrong here, for two reasons
we measured rather than assumed. The stock minimum is often already cheaper
than the rotated tool's own plain hybrid, so subtracting double counts. And on
two of our instances the rotated cost is *higher* than the plain one, so the
delta would be negative: at encapsulated weight 26 the rotated hybrid without
meet in the middle costs 190.88 against 188.28, the enlarged guessing set
outrunning the improved hit probability.

The theoretical ceiling on the gain is `log2 N`, which is 13 to 16 bits on these
rings. Every measured loss is well below it. That ceiling is our derivation from
the hit-probability amplification, not a theorem the paper states.

## Files

| file | contents |
|---|---|
| `estimate_security.py` | stock Lattice Estimator, six attacks per instance |
| `rot_estimate.py` | the rotated primal hybrid of ePrint 2026/279, both meet-in-the-middle heuristics, plus modulus and weight sweeps |
| `admissibility.py` | the order-four admissibility filter, `p = a^2 + b^2` with `max(a,b) > 2B` |
| `sets*.json` | parameter set definitions |
| `make_results.py` | regenerates `results/RESULTS.md` from the JSONL files; rerun after appending new runs |
| `results/RESULTS.md` | summary tables, including the ciphertext benchmarks |
| `results/*.jsonl` | one record per run, structured |
| `raw/` | unedited campaign logs, the source of everything above |
| `estimator_output.txt` | the first stock estimator campaign, kept for reference |

## Reproducing

Stock estimator:

```
git clone https://github.com/malb/lattice-estimator
git -C lattice-estimator checkout 3e48ef4
sage -python estimate_security.py --validate           # reproduce the baseline
sage -python estimate_security.py --sets sets.json
```

Rotated hybrid:

```
git clone https://github.com/TabOg/mlwe-hybrids
cd mlwe-hybrids/PrimalHybrid
git -C lattice_estimator checkout 3e48ef4              # the same commit
cd -
sage -python rot_estimate.py --sets sets.json --out results/rot.jsonl
sage -python rot_estimate.py --sweep-q 32768 1150 1265 5 --h 120
```

Both tools are pinned to the same `lattice_estimator` commit so that the two
columns of every table rest on one lattice-estimation codebase. Cost model
`RC.MATZOV`, GSA shape, `sigma = 3.2`, sparse ternary secrets, `m = n` samples
throughout.

## Independent check on our use of the attack tool

The author of ePrint 2026/279 computed the rotated square-root column for the
same instances. All eleven comparable values agree with ours to within 0.01
bits, and the optimal guessing dimensions agree exactly: 4893, 3890, 3208,
2627, 2130 and 1950 across the weight sweep. Two baseline values differ, by
0.05 and 0.21 bits, which we attribute to two code paths settling in different
local optima of the same optimiser; neither is the binding column for its row.

## Two results worth stating separately

**The gain is exactly zero when a non-hybrid attack is already cheapest.** This
is visible in the released code, where the rotated cost function short-circuits
to the plain estimate at guessing dimension zero, and it is visible in our data:
on the two evaluator microbenchmarks the optimal guessing dimension is 0 and 491
and the gain is 0.00 and 0.02 bits, while on the long-chain main-key rows the
dimension is in the thousands and the gain is large. Those two microbenchmark
rows sit far below 80 bits, but not because of this attack; they carry a chain
sized for a much larger ring and were never parameter proposals. We withdraw
them as such.

**Raising the main key weight does not restore the target.** Holding the chain
fixed while raising the weight is not self-consistent, because the noise bound
grows as `sqrt(h)` and the chain must lengthen to preserve the same
bootstrapping depth. Charging even the smallest chain growth consistent with
that law, weight 270 reaches only 79.07 bits and the curve is flat from 330; the
guessing dimension never approaches zero, so the attack is never switched off,
only made gradually less profitable. The recommendations therefore leave the
weight at the library default and carry the security target on the modulus,
which is also what the prior work does.
