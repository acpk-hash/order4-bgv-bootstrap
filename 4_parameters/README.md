# 4_parameters — parameter selection and concrete security

This folder contains everything behind the parameter sets of the paper
(Sec. 6.1, Tables 3 and 4; full version also Tables 9 and 10):
the pinned Lattice Estimator, our estimation scripts, every transcript, and
the recommended parameters in machine-readable form.

```
4_parameters/
  recommended_parameters.csv     the recommended sets (table below), one row per set
  lattice-estimator-3e48ef4/     Lattice Estimator, unmodified, commit 3e48ef4 (LGPLv3+)
  scripts/
    built_chain.py               per-attack security at the chain HElib built  (Table 4, λ*)
    estimate_security.py         stock estimator, six attacks per instance      (Tables 3-4)
    rot_estimate.py              rotated primal hybrid of ePrint 2026/279       (Table 3)
    admissibility.py             order-four admissibility per prime             (Thm 3, full Table 9)
    make_results.py              regenerates transcripts/summaries/RESULTS.md
    sets*.json                   instance definitions
  transcripts/
    built_chain/<set>.json       output of built_chain.py, nine sets
    summaries/*.jsonl            one record per estimator run; RESULTS.md summarises them
    raw/                         unedited campaign logs
  README_security_scripts.md     details of the security model and of the two tools
```

## 1. How a parameter set is chosen

1. **Plaintext prime and radix.** `p ≡ 1 (mod 4)` and `A` with `A² ≡ −1 (mod p)`
   (Theorem 3 (i)). `scripts/admissibility.py` lists `A` and the largest
   admissible `B` for every prime used.
2. **Noise bound.** The encapsulated key weight `h` fixes the bound `B`
   (`h = 12 → B = 17`, `h = 26 → B = 24` on general rings; `h = 26 → B = 13` on
   `m = 2^16`), and `(2B+1)² < p` must hold (Theorem 3 (ii)); `deg P_A = (2B+1)² − 2`.
3. **Ring.** `m = 2^16` or `2^17`, or a general `m = q₁q₂` with five distinct
   primes per level and enough slots.
4. **Security.** Two LWE instances per set (Ma et al., EUROCRYPT 2024, Sec. 4.3):
   the main key (`n = φ(m)`, full chain `Q`, weight `h' = 120`) and the
   encapsulated bootstrapping key (`n = φ(m)`, short modulus `q_boot`, weight `h`).
   The security λ* is the smaller of the two, and for each instance the minimum
   over the estimator's attacks at their own optimal parameters (uSVP, BDD, dual,
   primal hybrid with and without meet-in-the-middle and Babai). On power-of-two
   rings the rotated primal hybrid of ePrint 2026/279 is applied as well.
   Cost model `RC.MATZOV`, `σ = 3.2`, sparse ternary secrets, `m = n` samples.
5. **Built chain.** HElib does not build the chain it is asked for, so λ* of every
   set with a ciphertext run is recomputed at the chain HElib printed
   (`INFO: total bits = …`); this is the last column of Table 4.

## 2. Recommended parameters

From `recommended_parameters.csv`. `log2 Q` recommended = the chain length
recommended by the paper; built = the chain HElib actually built in the
ciphertext runs. "stock" is the stock Lattice Estimator on the main key,
"rotated" the same instance under ePrint 2026/279 (power-of-two rings only).

**Power-of-two rings** (Table 3): `p = 65537`, `h' = 120`, `h = 26`.

| set | m | slots | target | log2 Q recommended | log2 Q built | main, stock | main, rotated | encapsulated | **λ*** |
|---|---|---|---|---|---|---|---|---|---|
| po2-80 | 2^16 | 32768 | 80 | 1078 | 1066 | 95.6 | 86.7 | 138.7 | **86.7** |
| po2-100 | 2^17 | 32768 | 100 | 1676 | – | 117.8 | 106.0 | 149.6 | **106.0** |
| po2-128 | 2^17 | 32768 | 128 | 1277 | – | 143.6 | 131.4 | 149.6 | **131.4** |

**General cyclotomic rings** (Table 4): `h' = 120`,
five distinct primes per level. λ* is at the built chain.

| set | p | m | φ(m) | slots | h | log2 Q recommended | log2 Q built | main | encapsulated | **λ*** |
|---|---|---|---|---|---|---|---|---|---|---|
| 80-A | 4513 | 37637 | 36960 | 9240 | 12 | 1510 | 1524.2 | 81.1 | 87.6 | **80.5** |
| 80-B | 14401 | 43033 | 41140 | 3740 | 12 | 1680 | 1642.0 | 81.3 | 88.7 | **82.75** |
| 80-C | 13457 | 45193 | 44100 | 6300 | 12 | 1800 | 1786.8 | 81.4 | 89.7 | **81.94** |
| 80-D | 2521 | 50731 | 50112 | 2784 | 12 | 2040 | 2026.2 | 81.9 | 91.7 | **82.44** |
| 80-E | 65537 | 50731 | 50112 | 2784 | 12 | 2017 | – | 82.6 | 90.9 | – |
| 128-A | 13457 | 65047 | 62776 | 31388 | 26 | 1427 | 1406.9 | 128.5 | 158.5 | **129.81** |
| 128-B | 2917 | 53983 | 52488 | 13122 | 26 | 1184 | 1153.0 | 128.6 | 153.4 | **131.0** |
| 128-C | 15377 | 55459 | 53640 | 13410 | 26 | 1212 | 1193.1 | 128.4 | 154.1 | **129.90** |
| 128-D | 8101 | 57521 | 56832 | 7104 | 26 | 1287 | 1267.5 | 128.5 | 155.7 | **129.84** |
| 128-E | 14401 | 54359 | 51480 | 5148 | 26 | 1161 | 1132.6 | 128.5 | 153.0 | **130.77** |

80-E has no ciphertext run at its own chain (paper, Sec. 6.3). On every row with
a built chain the binding attack is the primal hybrid with meet-in-the-middle
and without Babai (TF); `transcripts/built_chain/<set>.json` has all six attacks.
These are the sets measured in Tables 5 and 6.

## 3. Reproducing

All scripts need SageMath (9 or later) and use the pinned estimator copy in this
folder by default (override with `LATTICE_ESTIMATOR=/path`).

```sh
cd 4_parameters/scripts

# λ* at the built chain, one set (~15-40 min each, single-threaded)
sage -python built_chain.py 80B > ../transcripts/built_chain/80B.json
python3 ../../3_experiments/tests/test_security.py      # compare all nine with the paper

# order-four radix and admissible B per prime (instant, no Sage)
python3 admissibility.py

# stock estimator on the instance definitions (hours)
sage -python estimate_security.py --sets sets_recommend.json --out ../transcripts/summaries/new.jsonl

# summary tables from the transcripts (instant, no Sage)
python3 make_results.py
```

`rot_estimate.py` additionally needs the ePrint 2026/279 tool
(`git clone https://github.com/TabOg/mlwe-hybrids`, see
`README_security_scripts.md`); it is not vendored because it carries no licence.

## 4. Status of the transcripts

| quantity | transcript | reproduced by us |
|---|---|---|
| Table 4 λ* (built chain) and full-version Table 10 | `transcripts/built_chain/` | all nine sets, every attack within 0.01 bit (`3_experiments/tests/test_security.py`) |
| Table 3 boundary chains, sq-root / estim. columns, rotated losses | `summaries/cross_verification.jsonl`, `rot_2026_279_twelve.jsonl`, `modulus_sweep.jsonl` | stored transcripts |
| Table 3 `λ_est` (95.6 / 117.8 / 143.6) | `raw/estimator_a1.jsonl` | stored transcript |
| Table 4 encapsulated column, sets 80-A to 80-E and 128-A | `estimator_output.txt` | stored transcript |
| Table 4 encapsulated column, sets 128-B to 128-E; Table 4 "main" column at the recommended chain; Table 3 caption values 111.2 and 139.9 | not in the transcripts | not reproduced here. In the paper these values exceed λ* of the same row by 6 to 29 bits, so they do not determine λ*; λ* itself (main key at the built chain) is reproduced for every set with a ciphertext run |
