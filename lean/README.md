# Machine-Checked Lean 4 Formalization
### *Order-Four Symmetry in BGV Bootstrapping — Faster Digit Extraction for Large Primes*

This directory contains a **Lean 4 + Mathlib** formalization of the algebraic core
of the paper. Every result below is proved with **zero `sorry`, zero `admit`,
zero `axiom`**; `lake build` completes with no errors.

## Build / reproduce

```bash
# Lean toolchain is pinned in lean-toolchain (leanprover/lean4:v4.31.0-rc2)
lake exe cache get          # fetch prebuilt Mathlib oleans
lake build                  # builds BootstrapAlgebra.OrderFour + .Character
```

A green build is the certificate: it re-checks every proof from the Mathlib
kernel up.

## Files

| file | contents |
|---|---|
| `BootstrapAlgebra/OrderFour.lean` | order-four element, support stability, the character filter, the **Order-Four Character Filter** theorem, factored form |
| `BootstrapAlgebra/Character.lean` | the Section-4 character-decomposition framework: **Isotypic decomposition** and the **Monomial Support Constraint** (covariance) |

## Interface: paper statement ↔ Lean name

Everything is formalized at the **coefficient level**: the substitution operator
`T_A : X^k ↦ A^k X^k` is modeled by its action on a coefficient sequence
`c : ℕ → K`, and the digit-extraction polynomial by its coefficients.

### `OrderFour.lean`  (namespace `OrderFour`)

| paper statement | Lean name |
|---|---|
| `A² = -1 ⟹ A⁴ = 1` (order-four element) | `sq_neg_one_pow_four` |
| `Aᵏ = A^(k mod 4)` (4-periodic action) | `pow_mod_four` |
| `A³ = -A` | `pow_three` |
| **Lemma (Support stability)** `A·(hA+ℓ) = ℓA + (-h)` — the quarter-turn `(h,ℓ)↦(ℓ,-h)`, giving `A·S_A = S_A` | `stability_pointwise` |
| the filter *passes* `k ≡ 3 (mod 4)`: `Aᵏ + A = 0` | `filter_factor_eq_zero` |
| the filter *kills* `k ≢ 3 (mod 4)`: `Aᵏ + A ≠ 0` | `filter_factor_ne_zero` |
| character filter as an iff: `Aᵏ + A = 0 ⇔ k ≡ 3 (mod 4)` | `filter_factor_zero_iff` |
| `(P.comp (C A · X)).coeff k = Aᵏ · P.coeff k` (bridge) | `coeff_comp_C_mul_X` |
| **Theorem (Order-Four Character Filter)** — coefficient form: from `(Aᵏ+A)·cₖ = A·[k=1]`, support ⊆ `{1}∪{k≡3 mod 4}` and `c₁ = 1/2` | `order_four_filter_coeff` |
| **Theorem (Order-Four Character Filter)** — polynomial form: from `P(AX)+A·P(X)=AX` | `order_four_filter` |
| support corollary: `supp(P) ⊆ {1}∪{k≡3 mod 4}` | `support_subset` |
| **Theorem (Factored evaluation form)** `P = ½X + X³·Q(X⁴)` (coefficient statement) | `factored_form` |

### `Character.lean`  (namespace `OrderFour.Character`)

| paper statement | Lean name |
|---|---|
| `A^r = 1 ⟹ Aᵏ = A^(k mod r)` (order-`r` action) | `pow_mod_order` |
| **Theorem (Isotypic decomposition)** — eigenvalue on `Xᵏ` depends only on `k mod r` | `isotypic_eigenvalue` |
| diagonal action of `T_A` on coefficients | `isotypic_action` |
| covariance factoring `(Aᵏ - ωʲ)·cₖ = hₖ` | `covariance_factor` |
| **Prop (Monomial Support Constraint) (i)** — determined components `cₖ = hₖ/(Aᵏ-ωʲ)` for `Aᵏ ≠ ωʲ` | `covariance_determined` |
| **Prop (Monomial Support Constraint) (ii)** — consistency `Aᵏ = ωʲ ⟹ hₖ = 0` | `covariance_consistency` |
| the full dichotomy over every index `k` | `covariance_dichotomy` |

## Notes

- The formalization also *confirms a reviewer observation*: the functional
  equation `P(AX)+A·P(X)=AX` **alone** forces the support into `{1}∪{k≡3 mod 4}`
  and `c₁ = 1/2` (see `order_four_filter`), so the separate oddness argument is
  not logically required for the support characterization.
- The coefficient-level modeling is faithful to the paper's own proof, whose
  key step (Theorem, Step 2) is exactly the coefficient comparison
  `(Aᵏ + A)·cₖ = A·[k=1]`.
