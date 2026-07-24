/-
  Character.lean

  Formal verification (Lean 4 + Mathlib) of the "Section 4 character-decomposition
  framework" for BGV bootstrapping (Asiacrypt paper), generalizing the order-four
  special case of `OrderFour.lean` to arbitrary order `r`.

  Everything is modelled at the *coefficient level*, exactly as in
  `OrderFour.lean`'s `order_four_filter_coeff`: a formal element `f = Σ c_k X^k`
  is represented by its coefficient sequence `c : ℕ → K`, and the substitution
  operator `T_A : X^k ↦ A^k X^k` acts on coefficients by `c_k ↦ A^k · c_k`.

  This file formalizes, with *no* `sorry`:

  * `pow_mod_order`
        For an order-`r` element `A` (`A ^ r = 1`, `r ≠ 0`), the power `A ^ k`
        depends only on `k mod r`.        (generalizes `OrderFour.pow_mod_four`)

  * `isotypic_eigenvalue`, `isotypic_action`
        (Theorem "Isotypic decomposition", coefficient form)
        `V = ⊕_{j=0}^{r-1} V_j` with `V_j = {f : T_A f = ω^j f}`.  Since
        `T_A (X^k) = A^k X^k = A^{k mod r} X^k`, the eigenvalue attached to the
        monomial `X^k` depends only on `k mod r`, so `X^k ∈ V_{k mod r}`.

  * `covariance_determined`, `covariance_consistency`, `covariance_dichotomy`
        (Proposition "Monomial support constraint")
        If `T_A f = ω^j f + h` (coefficient-wise `A^k c_k = ω^j c_k + h_k`), then
        for every `k` with `A^k ≠ ω^j` the coefficient is uniquely determined,
        `c_k = h_k / (A^k − ω^j)`; and where `A^k = ω^j` the relation forces the
        consistency condition `h_k = 0`, leaving `c_k` free.
-/
import Mathlib

namespace OrderFour.Character

open Polynomial

/-! ## 1. Order-`r` element:  `A ^ r = 1`  ⟹  `A ^ k` is `r`-periodic in `k` -/
section Ring
variable {R : Type*} [CommRing R] {A : R}

/-- The order-`r` action is `r`-periodic in the exponent: `A ^ k = A ^ (k % r)`.
    This generalizes `OrderFour.pow_mod_four` (the `r = 4` case) to an arbitrary
    order-`r` element `A` (with `A ^ r = 1` and `r ≠ 0`). -/
lemma pow_mod_order (hA : A ^ r = 1) (hr : r ≠ 0) (k : ℕ) : A ^ k = A ^ (k % r) := by
  conv_lhs => rw [← Nat.div_add_mod k r]
  rw [pow_add, pow_mul, hA, one_pow, one_mul]

end Ring

/-! ## 2. Theorem "Isotypic decomposition" (coefficient / eigenvalue form)

The character-decomposition framework writes the ambient space as
`V = ⊕_{j=0}^{r-1} V_j`, where `V_j = {f : T_A f = ω^j f}` is the `ω^j`-eigenspace
of the substitution operator `T_A : X^k ↦ A^k X^k`.  Because `A` is a primitive
`r`-th root of unity (`A ^ r = 1`), we have `T_A (X^k) = A^k X^k = A^{k mod r} X^k`,
so the monomial `X^k` is an eigenvector whose eigenvalue depends only on `k mod r`;
hence `X^k ∈ V_{k mod r}`.  The two lemmas below make this precise at the level of
eigenvalues. -/
section Isotypic
variable {R : Type*} [CommRing R] {A : R}

/-- The eigenvalue `A ^ k` attached to the monomial `X^k` under `T_A` depends only
    on the residue `k mod r`.  Concretely, whenever `k ≡ l (mod r)` the two
    monomials `X^k` and `X^l` carry the *same* eigenvalue, so they lie in the same
    isotypic component `V_{k mod r}`.                (Theorem "Isotypic decomposition") -/
lemma isotypic_eigenvalue (hA : A ^ r = 1) (hr : r ≠ 0) {k l : ℕ}
    (hkl : k % r = l % r) : A ^ k = A ^ l := by
  rw [pow_mod_order hA hr k, pow_mod_order hA hr l, hkl]

/-- Coefficient form of the isotypic action of `T_A : X^k ↦ A^k X^k`.  The operator
    scales the `k`-th coefficient by the eigenvalue `A ^ k`, which by `r`-periodicity
    equals `A ^ (k % r)`.  Thus `T_A` acts on the coefficient sequence `c` diagonally,
    with the eigenvalue at index `k` depending only on `k mod r`. -/
lemma isotypic_action (hA : A ^ r = 1) (hr : r ≠ 0) (c : ℕ → R) (k : ℕ) :
    A ^ k * c k = A ^ (k % r) * c k := by
  rw [pow_mod_order hA hr k]

end Isotypic

/-! ## 3. Proposition "Monomial support constraint"

Suppose the covariance relation `T_A f = ω^j f + h` holds, i.e. coefficient-wise
`A^k · c_k = ω^j · c_k + h_k` for every `k`.  Reading off the `k`-th coefficient,
this is `(A^k − ω^j)·c_k = h_k`.  There are two regimes:

  (i)  If `A^k ≠ ω^j` (`X^k` sits in a different isotypic component than the target
       `ω^j`), the coefficient is *uniquely determined*: `c_k = h_k / (A^k − ω^j)`.

  (ii) If `A^k = ω^j` (`X^k` is in the target component `V_j`), the relation forces
       the *consistency condition* `h_k = 0`, and `c_k` is left completely free.
-/
section Field
variable {K : Type*} [Field K] {A ω : K} {j : ℕ}

/-- The covariance relation, read off coefficient by coefficient, is equivalent to the
    factored linear equation `(A^k − ω^j)·c_k = h_k`. -/
lemma covariance_factor (c h : ℕ → K)
    (hcov : ∀ k, A ^ k * c k = ω ^ j * c k + h k) (k : ℕ) :
    (A ^ k - ω ^ j) * c k = h k := by
  have := hcov k
  linear_combination this

/-- Proposition "Monomial support constraint", part (i).  If the covariance relation
    `A^k c_k = ω^j c_k + h_k` holds and `A^k ≠ ω^j`, then the coefficient is uniquely
    determined by `c_k = h_k / (A^k − ω^j)`. -/
theorem covariance_determined (c h : ℕ → K)
    (hcov : ∀ k, A ^ k * c k = ω ^ j * c k + h k) {k : ℕ} (hne : A ^ k ≠ ω ^ j) :
    c k = h k / (A ^ k - ω ^ j) := by
  have hfac : (A ^ k - ω ^ j) * c k = h k := covariance_factor c h hcov k
  have hden : A ^ k - ω ^ j ≠ 0 := sub_ne_zero.mpr hne
  rw [eq_div_iff hden]
  linear_combination hfac

/-- Proposition "Monomial support constraint", part (ii), consistency direction.
    Where `A^k = ω^j` (the monomial `X^k` lies in the target isotypic component
    `V_j`), the covariance relation forces the consistency condition `h_k = 0`.
    The coefficient `c_k` itself is unconstrained (it is the free parameter of the
    eigenspace `V_j`). -/
theorem covariance_consistency (c h : ℕ → K)
    (hcov : ∀ k, A ^ k * c k = ω ^ j * c k + h k) {k : ℕ} (heq : A ^ k = ω ^ j) :
    h k = 0 := by
  have hfac : (A ^ k - ω ^ j) * c k = h k := covariance_factor c h hcov k
  rw [heq, sub_self, zero_mul] at hfac
  exact hfac.symm

/-- The full dichotomy of the "Monomial support constraint": for every index `k`,
    either `A^k ≠ ω^j` and the coefficient is pinned down as `h_k / (A^k − ω^j)`,
    or `A^k = ω^j` and the inhomogeneous term satisfies the consistency `h_k = 0`
    (with `c_k` free). -/
theorem covariance_dichotomy (c h : ℕ → K)
    (hcov : ∀ k, A ^ k * c k = ω ^ j * c k + h k) (k : ℕ) :
    (A ^ k ≠ ω ^ j ∧ c k = h k / (A ^ k - ω ^ j)) ∨ (A ^ k = ω ^ j ∧ h k = 0) := by
  by_cases hk : A ^ k = ω ^ j
  · exact Or.inr ⟨hk, covariance_consistency c h hcov hk⟩
  · exact Or.inl ⟨hk, covariance_determined c h hcov hk⟩

end Field

end OrderFour.Character
