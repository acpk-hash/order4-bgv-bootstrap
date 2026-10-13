# 1_implementation — our code

| path | what it is | upstream it starts from |
|---|---|---|
| `helib_order4/` | HElib with the order-four digit-extraction evaluator, used for the general rings and `m = 50731` (Tables 6 and 8) | HElib `3e337a6` + Ma et al. `bgvboot-final.patch` (both in `2_baselines/`) |
| `driver_order4/` | bootstrapping driver for `helib_order4`: parameter sets of Tables 6 and 8, the option `bits=` | Ma et al. Large-p driver `fatboot.cpp` (MIT) |
| `patches/helib-ntt-order4.patch` | the same evaluator ported to Ma et al.'s power-of-two HElib (2 files, +122 lines), Table 5 | HElib `3e337a6` + Ma et al. `github_20241003.patch` |
| `patches/ma-ntt-fatboot-order4.patch` | parameter sets of Table 5 and the options `bits=`, `mainhwt=` for Ma et al.'s NTT driver | Ma et al. NTT driver `fatboot.cpp` |
| `patches/helib-auxradix-order4.REFERENCE.patch` | for reading only: the complete diff of `helib_order4/` against its upstream | – |
| `butterfly/` | Section 5: the Step2 (CoeffToSlot) block as a mixed-radix/Rader butterfly on ciphertexts (`enc_butterfly.cpp`), the plaintext derivation of its stages (`derive_butterfly.py` → `stages.txt`), and the HElib matrix dump it is checked against | builds against `helib_order4/` |

## Where the method lives

* `helib_order4/src/extractDigits.cpp`: the evaluator. With
  `HELIB_AUX_ORDER4_EVAL=1` and an order-four radix, `P_A` is evaluated as
  `X/2 + X^3·Q(X^4)`: one squaring chain to `X^4`, Paterson–Stockmeyer on `Q`
  (degree `(deg P_A − 3)/4`), one multiplication by `X^3`. Otherwise, or if
  `A² ≢ −1 (mod p)`, it prints a message and evaluates `P_A` exactly as the
  baseline does. Search for `HELIB_AUX_ORDER4_EVAL`.
* `helib_order4/src/recryption.cpp`: `HELIB_EXPLICIT_AUX=A` sets the auxiliary
  radix instead of the driver default `2B+1` (search for `HELIB_EXPLICIT_AUX`).

All other changes in `helib_order4/` (see the reference diff) are experimental
options that are off unless their environment variable is set, and no
experiment of the paper sets them; the one exception is
`HELIB_AUX_PARALLEL_COEFF2SLOT`, used only by the "combined" row of Table 8.

## OFF and ON come from the same binary

Each executable contains both evaluators; the environment chooses. So OFF
(control) and ON (treatment) share the binary, the keys' distribution, the
modulus chain and every other stage of bootstrapping. The baselines' original
code is built separately from `2_baselines/` (group G0 of
`3_experiments/compare_baselines.sh`) to show that our executables with the
evaluator OFF behave as the original code (group G1).
