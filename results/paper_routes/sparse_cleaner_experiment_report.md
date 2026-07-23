# Sparse Exact Cleaner Experiment Report

Date: 2026-04-27

This report records the stable experiment route for the main polynomial/null
polynomial contribution.  All scripts and outputs are under
`.`, so the baseline code collection is not
modified.

## Claim Boundary

For fixed `p`, bounded support parameter `B`, and auxiliary radix `A`, the
exact cleaner is

```text
C_A(hi*A + lo) = lo mod p,       hi, lo in [-B, B].
```

The target support size is `(2B+1)^2`.  When the support points are distinct,
the degree-`< |S|` interpolant is unique.  Therefore the current optimization
is not a lower-degree exact polynomial for the same support.  The stable claim
is:

```text
choose algebraically structured aux so the unique exact cleaner is much
sparser and cheaper to evaluate in the existing BGV bootstrapping pipeline.
```

## Scripts Added

| script | purpose |
| --- | --- |
| `scripts/cleaner_degree_certificate.py` | certificate that the canonical exact cleaner degree is minimal for the fixed support |
| `scripts/structured_aux_cleaner_experiment.py` | structured search over low-order/algebraic aux candidates |
| `scripts/cleaner_term_structure.py` | inspect the monomial support pattern of each cleaner |
| `scripts/parse_helib_bootstrap_timings.py` | parse HElib fatboot timing logs and compute repeat statistics |

## Degree Certificates

### Target HElib Parameter Set

`p=65537`, `B=17`, support size `1225`.

| aux | tag | degree | nonzero terms | degree-minimal? | output |
| ---: | --- | ---: | ---: | --- | --- |
| 35 | generic baseline | 1223 | 612 | yes | `cleaner_degree_certificate_aux35.json` |
| 256 | `sqrt(-1)`, order 4 | 1223 | 307 | yes | `cleaner_degree_certificate_aux256.json` |

The degree stays `1223` in both cases.  The useful change is the monomial
support: `612 -> 307` nonzero coefficients, a `49.84%` sparsity reduction.

### Small Complete-Search Sanity Check

`p=257`, `B=5`, support size `121`.

| aux | tag | degree | nonzero terms | degree-minimal? | output |
| ---: | --- | ---: | ---: | --- | --- |
| 35 | generic | 119 | 60 | yes | `cleaner_degree_certificate_p257_B5_aux35.json` |
| 16 | `sqrt(-1)`, order 4 | 119 | 31 | yes | `cleaner_degree_certificate_p257_B5_aux16.json` |

This smaller field reproduces the same phenomenon: the best `sqrt(-1)` radix
does not lower degree, but it roughly halves the number of terms.

## Search Results

### Complete Search At Small Size

Command:

```bash
python3 scripts/aux_cleaning_global_sweep.py \
  --p 257 --B 5 --mode full --top 20 --json \
  --output results/paper_routes/full_aux_search_p257_B5.json
```

Result:

| p | B | candidates | valid | best aux | best nonzero terms | best generic terms |
| ---: | ---: | ---: | ---: | --- | ---: | ---: |
| 257 | 5 | 255 | 130 | 16, 241 | 31 | 59 |

The full search selects exactly the order-4 roots `A^2=-1 mod 257` as the best
structured radices.

### Structured Search At Target Size

Command:

```bash
timeout 240s python3 scripts/structured_aux_cleaner_experiment.py \
  --p 65537 --B 17 --max-order 64 \
  --aux 35 256 65281 65502 --top 30 --json \
  --output results/paper_routes/structured_aux_search_p65537_B17.json
```

Result:

| p | B | candidates | valid | best aux | best nonzero terms | baseline terms |
| ---: | ---: | ---: | ---: | --- | ---: | ---: |
| 65537 | 17 | 65 | 44 | 256, 65281 | 307 | 612 |

No scanned structured candidate beats `aux=256`.  The residue-class equivalent
candidate `65281 = -256 mod 65537` has the same polynomial over `F_p`, but it
is not the practical HElib choice because HElib uses the integer aux in the
mod-up/overflow scaling.  The practical radix is the small positive root
`aux=256`.

## Monomial Support Structure

Output files:

```text
results/paper_routes/cleaner_term_structure_p257_B5.json
results/paper_routes/cleaner_term_structure_p65537_B17.json
```

For the target parameter set:

| aux | degree | terms | exponents `mod 4` |
| ---: | ---: | ---: | --- |
| 35 | 1223 | 612 | `1 mod 4`: 306, `3 mod 4`: 306 |
| 256 | 1223 | 307 | `1 mod 4`: 1, `3 mod 4`: 306 |

So the order-4 radix removes almost all `1 mod 4` monomials while preserving
exactness on the bounded support.  This gives a concrete mathematical
explanation for the sparse exact cleaner.

## Real HElib Ciphertext Benchmark

Existing target ciphertext logs:

```text
results/baseline_BGV_Boot_for_Large_p_i4_repeat3.log
results/upload_aux_default_repeat1.log
results/upload_aux_256_repeat1.log
results/upload_aux_default_repeat3.log
results/upload_aux_256_repeat3.log
results/baselines_BGV_Boot_for_Large_p_vs_aux256_sparse_repeat3.json
results/paper_routes/helib_aux35_vs_aux256_repeat3_summary.json
```

All listed ciphertext runs finish with `### bts finished, everything ok ###`.
For the paper-facing comparison, the direct baseline is
`/home/user/experiments/baselines/BGV-Boot-for-Large-p`, whose repeat-3
log is `results/baseline_BGV_Boot_for_Large_p_i4_repeat3.log`.  The
`upload_aux_default_repeat3.log` file is retained as an internal same-source
reproduction, but it is not the primary comparison object.

### Single-Run Baseline

| run | aux | linear1 | linear2 | extract | total |
| --- | ---: | ---: | ---: | ---: | ---: |
| baseline repeat1 | 35 | 4.488775s | 44.084850s | 156.793361s | 207.677903s |
| sparse cleaner repeat1 | 256 | 4.458518s | 43.957204s | 144.869405s | 195.591277s |

Observed improvement:

| metric | improvement |
| --- | ---: |
| nonzero monomial terms | 49.84% fewer |
| extraction time | 7.60% faster |
| total bootstrapping time | 5.82% faster |

### Repeat-3 Stability Run

Commands:

```bash
cd /home/user/experiments/baselines/BGV-Boot-for-Large-p
timeout 1800s ./build/fatboot i=4 h=12 t=-1 newbts=1 repeat=3 newks=1 \
  > /home/user/order4-artifact/results/baseline_BGV_Boot_for_Large_p_i4_repeat3.log 2>&1

HELIB_ZZX_CACHE_DIR=./cache/saved_ZZX \
HELIB_EXPLICIT_AUX=256 \
timeout 2400s \
./src/BGV-Boot-auxradix-opt/build/fatboot \
  i=4 h=12 t=-1 newbts=1 newks=1 thick=0 repeat=3 \
  > results/upload_aux_256_repeat3.log 2>&1
```

Per-run timings:

| run | aux | linear1 | linear2 | extract | total |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 35 direct Ma baseline | 4.614521s | 45.368782s | 160.520337s | 212.870691s |
| 2 | 35 direct Ma baseline | 4.321101s | 45.314570s | 161.361991s | 213.740190s |
| 3 | 35 direct Ma baseline | 4.348797s | 45.123198s | 160.467474s | 212.696840s |
| 1 | 256 | 4.519253s | 44.159891s | 146.529575s | 197.502668s |
| 2 | 256 | 4.262789s | 44.465031s | 148.003047s | 199.430020s |
| 3 | 256 | 4.360044s | 44.911091s | 149.197179s | 201.258364s |

Repeat-3 mean and sample standard deviation:

| aux | linear1 mean | linear2 mean | extract mean | total mean | total std |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 35 direct Ma baseline | 4.428140s | 45.268850s | 160.783267s | 213.102573s | 0.558992s |
| 256 | 4.380695s | 44.512004s | 147.909934s | 199.397017s | 1.878065s |

Repeat-3 observed improvement:

| metric | improvement |
| --- | ---: |
| extraction time | 8.01% faster |
| total bootstrapping time | 6.43% faster |

The extraction speedup is smaller than the term-count reduction because the
encrypted polynomial evaluation is still dominated by ciphertext multiplications
and relinearization/mod-down work, not only by the number of plaintext
coefficients.

## Paper Position

Use this as the main experimental contribution:

```text
global/structured search over aux-induced bounded supports finds an order-4
radix whose exact cleaner is degree-minimal but nearly half as sparse, and this
translates into a real HElib BGV bootstrapping speedup on the target parameter
set.
```

Do not claim:

```text
lower minimal degree for the same fixed support;
global optimality over all bootstrapping pipelines;
linear-transform algebraic speedup from the negative LC-AKS attempts.
```

The parallel `coeffToSlot` executor can remain an implementation supplement,
but it should be separated from the main mathematical claim.
