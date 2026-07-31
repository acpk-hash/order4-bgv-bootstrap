# Benchmark

Reproduces the ciphertext bootstrapping measurements in the paper: the
order-four evaluator switched on and off in an otherwise identical pipeline, so
that the ratio between two runs isolates the contribution.

## Why this is a wrapper and not a standalone driver

The order-four evaluator needs an auxiliary modulus that is **not** a power of
the plaintext prime. Stock HElib chooses that modulus in `RecryptData::setAE`,
which only considers powers of `p`; on a power-of-two ring at `p = 65537` it
fails outright with

```
helib::RuntimeError: setAE: cannot find suitable e
```

before a context can even be constructed. The fork under `../baselines/` extends
that choice, and the `t = -1` flag together with `HELIB_EXPLICIT_AUX` is how the
extension is selected. A driver written against stock HElib therefore does not
merely produce different numbers, it does not run at all. We tried; that is why
this directory ships a wrapper around the fork's own executable rather than a
fresh driver.

## What is here

| file | contents |
|---|---|
| `parameter.h` | every recommended set with its security under both estimators, the plaintext primes and their order-four radices, and the sets we withdraw |
| `run_bench.sh` | expands a set name into the exact `fatboot` invocation used for the reported numbers |
| `README.md` | this file |

## Building the executable it drives

```
cd ../baselines/bgv-bootstrapping-with-homomorphic-NTT
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -Dhelib_DIR=<prefix>/share/cmake/helib
cmake --build build -j
```

`<prefix>` is the HElib install carrying the order-four evaluator.

## Running

```
./run_bench.sh list
./run_bench.sh po2-80 off
./run_bench.sh po2-80 on
./run_bench.sh po2-80 on 13457        # a different plaintext prime
FATBOOT=/path/to/fatboot ./run_bench.sh po2-80 on
TRIALS=1 ./run_bench.sh po2-80 on     # default is 3
```

Compare two runs differing only in `off` and `on`. A verification run of
`po2-80 on` reproduces the recorded numbers: chain built at 1066 bits, 32768
slots, digit extraction 9.94 s against the 9.93 s recorded.

## Three failure modes that do not announce themselves

**The order-four radix is per prime.** `A` must satisfy `A^2 = -1 mod p`.
Passing one prime's radix to another does not error: the evaluator falls back to
the generic path and the run reports a speedup near 1.0, which reads as a
negative result rather than a misconfiguration. We lost a batch of nine primes to
exactly this. `run_bench.sh` refuses primes it has no radix for.

**The requested `bits` is not the chain you get.** HElib builds a longer chain
than requested, by roughly 300 bits at `m = 2^16` and 210 at `m = 2^17`, and the
offset is not a property of the ring alone: the same request at a different `t`
gives a different chain. Always quote the built value, which the run prints.

**A power-of-two ring with `p = 1 mod 4` needs a two dimensional hypercube.**
Generators `-1` of order 2 together with `5`. A one dimensional form compiles,
builds a context, and then fails inside the linear transform with `p = 1 mod 4
should have multiple blocks`. For `p = 3 mod 4` the requirement is the opposite.

## Sets not covered by this script

**`po2-100` and `po2-128`** live on `m = 2^17`, for which the baseline's
parameter table has no entry: it carries `2^14`, `2^15`, `2^16` and, in an array
that is not selected, `2^18`. We added and validated one. Appended to
`po2_real_params`, so that the existing indices keep their meaning, it is

```cpp
Parameters( // d = 2, m = 131072
        65537, 1, 3, 980, 0,
            1, 0, 1, 100,
            0, 0, force_chen_han, 0, 0,
        {-1, 5}, {2, 16384}, {131072}, 1, 5),
```

Validated by a run reporting `m = 131072`, `ordP = 2`, 32768 slots and a
bootstrap completing without error.

**General cyclotomic sets** need their own `mvec`, `gens` and `ords`, which
depend on `m = q1*q2` and on `ord_m(p)`. They live in the baseline's own
parameter table and are selected by index there.

## Resource cost

Measured peak resident set for one run, single threaded:

| ring | slots | peak RSS | wall per run |
|---|---:|---:|---:|
| `m = 2^16`, p = 65537 | 32768 | 31 GB | about 10 min |
| `m = 2^16`, p = 14401 | 32 | 16 GB | about 6 min |
| `m = 2^16`, p = 2521 | 4 | 119 GB | about 10 min |
| `m = 2^17`, p = 65537 | 32768 | 75 GB | about 25 min |

Memory scales inversely with the slot count, not with the ring. Size any
parallel scheduling from the low slot cases; taking `p = 65537` as
representative understates the worst case by a factor of four.

## Reproducing the reported numbers

Security figures in `parameter.h` come from `../security`, which carries the
estimator drivers and the raw per run output. Timing figures come from this
wrapper; the measured values are tabulated in `../security/results/RESULTS.md`.
