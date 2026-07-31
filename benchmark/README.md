# Benchmark

Ciphertext bootstrapping on the revised parameter sets, with the order-four
evaluator switched on and off in the identical pipeline so that the ratio
isolates the contribution.

## What is here

| file | contents |
|---|---|
| `parameter.h` | every recommended set, with its security under both estimators, the plaintext primes and their order-four radices, and the sets we withdraw |
| `bench_boot.cpp` | the driver: builds the context, bootstraps, verifies decryption, reports per stage timings |
| `CMakeLists.txt` | build |

The parameter table and the evaluator selection live here rather than in the
caller's shell, because a mismatch between the two is easy to miss and produces
a run that silently measures nothing. See the note on the radix below.

## Building

Bootstrapping cannot be vendored, so HElib is the one external dependency. It
must be the build carrying the order-four evaluator.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -Dhelib_DIR=<prefix>/share/cmake/helib
cmake --build build -j
```

## Running

```
./build/bench_boot list
./build/bench_boot set=po2-80 mode=off trials=3
./build/bench_boot set=po2-80 mode=on  trials=3
./build/bench_boot set=po2-80 mode=on  p=13457 trials=3
```

`mode=off` selects the bounded-support digit-extraction evaluator, `mode=on`
the order-four evaluator. Compare two runs that differ only in `mode`.

## Three things that will bite

**The order-four radix is per prime.** `A` must satisfy `A^2 = -1 mod p`. Passing
one prime's radix to another does not fail loudly: the evaluator falls back to
the generic path and the run reports a speedup of about 1.0, which looks like a
negative result rather than a configuration error. The driver refuses primes it
has no radix for, and `parameter.h` carries the table.

**The requested `bits` is not the chain you get.** HElib builds a longer chain
than requested, by roughly 300 bits at `m = 2^16` and 210 at `m = 2^17`, and the
offset is not a property of the ring alone. Always quote the built value, which
the driver prints. The `bits_request` field in `parameter.h` is the value that
produced the recorded `log2Q` on our machine, not a portable constant.

**Power-of-two rings need a two dimensional hypercube when `p = 1 mod 4`.** The
generators are `-1` of order 2 together with `5`. A one dimensional form
compiles, builds a context and then fails inside the linear transform with
`p = 1 mod 4 should have multiple blocks`. For `p = 3 mod 4` the requirement is
the opposite. The driver sets this correctly for the power-of-two sets.

## Resource cost

Measured peak resident set for one run, single threaded:

| ring | slots | peak RSS | wall per run |
|---|---:|---:|---:|
| `m = 2^16`, p = 65537 | 32768 | 31 GB | about 10 min |
| `m = 2^16`, low slot primes | 4 to 32 | 16 to 119 GB | 6 to 10 min |
| `m = 2^17` | 32768 | 75 GB | about 25 min |

Memory scales inversely with the slot count, not with the ring: `p = 2521` with
4 slots peaks at 119 GB while `p = 14401` with 32 slots peaks at 16 GB. Size any
parallel scheduling from the low slot cases, not from `p = 65537`.

## Reproducing the reported numbers

The security figures in `parameter.h` come from `../security`, which carries the
estimator drivers and the raw per run output. The timing figures come from this
benchmark; the measured values are tabulated in `../security/results/RESULTS.md`.
