/-
  OrderFour.lean

  Formal verification (Lean 4 + Mathlib) of the algebraic core of

    "Order-Four Symmetry in BGV Bootstrapping:
     Faster Digit Extraction for Large Primes"

  This file formalizes, with *no* `sorry`, the results underpinning the
  paper's main (digit-extraction) contribution:

  * `sq_neg_one_pow_four`, `pow_mod_four`
        A² = -1 forces A⁴ = 1, i.e. A is an order-four element, and Aᵏ
        depends only on k mod 4.                       (order-four action)

  * `stability_pointwise`
        A·(hA + ℓ) = ℓA + (-h), i.e. multiplication by A is the quarter
        turn (h, ℓ) ↦ (ℓ, -h).           (Lemma "Support stability", core)

  * `filter_factor_eq_zero`, `filter_factor_ne_zero`, `filter_factor_zero_iff`
        Over a field of characteristic ≠ 2 with A² = -1,
             Aᵏ + A = 0  ↔  k ≡ 3 (mod 4).            (the "character filter")

  * `order_four_filter_coeff`  and  `order_four_filter`
        (Theorem "Order-Four Character Filter")
        The functional equation  P(AX) + A·P(X) = A·X  forces every
        coefficient cₖ with k ≠ 1 and k ≢ 3 (mod 4) to vanish, and c₁ = 1/2.

  * `factored_form`
        (Theorem "Factored evaluation form")
        Consequently  P = (1/2)·X + X³ · Q(X⁴)  for some Q.
-/
import Mathlib

namespace OrderFour

open Polynomial

/-! ## 1. Order-four element:  A² = -1  ⟹  A⁴ = 1 -/
section Ring
variable {R : Type*} [CommRing R] {A : R}

/-- If `A² = -1` then `A⁴ = 1`: `A` has multiplicative order dividing 4. -/
lemma sq_neg_one_pow_four (hA : A ^ 2 = -1) : A ^ 4 = 1 := by
  have h : A ^ 4 = (A ^ 2) ^ 2 := by ring
  rw [h, hA]; norm_num

/-- The order-four action is 4-periodic in the exponent: `Aᵏ = A^(k mod 4)`. -/
lemma pow_mod_four (hA : A ^ 2 = -1) (k : ℕ) : A ^ k = A ^ (k % 4) := by
  conv_lhs => rw [← Nat.div_add_mod k 4]
  rw [pow_add, pow_mul, sq_neg_one_pow_four hA, one_pow, one_mul]

/-- `A³ = -A` when `A² = -1`. -/
lemma pow_three (hA : A ^ 2 = -1) : A ^ 3 = -A := by
  have : A ^ 3 = A ^ 2 * A := by ring
  rw [this, hA]; ring

/-- Lemma "Support stability" (algebraic core).  Multiplication by an
    order-four radix `A` (i.e. `A² = -1`) is the quarter turn
    `(h, ℓ) ↦ (ℓ, -h)` on the `A`-adic grid:  `A·(hA + ℓ) = ℓA + (-h)`.
    Since `[-B,B]²` is symmetric this shows `A · S_A = S_A`. -/
lemma stability_pointwise (hA : A ^ 2 = -1) (h l : R) :
    A * (h * A + l) = l * A + (-h) := by
  have : A * (h * A + l) = h * A ^ 2 + l * A := by ring
  rw [this, hA]; ring

end Ring

/-! ## 2. The character filter over a field of characteristic ≠ 2 -/
section Field
variable {K : Type*} [Field K] {A : K}

lemma A_ne_zero (hA : A ^ 2 = -1) : A ≠ 0 := by
  intro h; rw [h] at hA; norm_num at hA

lemma A_ne_one (hA : A ^ 2 = -1) (h2 : (2 : K) ≠ 0) : A ≠ 1 := by
  intro h; rw [h, one_pow] at hA
  exact h2 (by linear_combination hA)

lemma A_ne_neg_one (hA : A ^ 2 = -1) (h2 : (2 : K) ≠ 0) : A ≠ -1 := by
  intro h; rw [h, neg_one_sq] at hA
  exact h2 (by linear_combination hA)

/-- The filter *passes* exactly the exponents `k ≡ 3 (mod 4)`:
    for those, the coefficient factor `Aᵏ + A` vanishes. -/
lemma filter_factor_eq_zero (hA : A ^ 2 = -1) {k : ℕ} (hk : k % 4 = 3) :
    A ^ k + A = 0 := by
  rw [pow_mod_four hA, hk, pow_three hA]; ring

/-- The filter *kills* every other exponent: for `k ≢ 3 (mod 4)` the factor
    `Aᵏ + A` is nonzero, so (over a field) the coefficient must vanish. -/
lemma filter_factor_ne_zero (hA : A ^ 2 = -1) (h2 : (2 : K) ≠ 0)
    {k : ℕ} (hk : k % 4 ≠ 3) : A ^ k + A ≠ 0 := by
  rw [pow_mod_four hA]
  have hAne : A ≠ 0 := A_ne_zero hA
  have hlt : k % 4 < 4 := by omega
  rcases (show k % 4 = 0 ∨ k % 4 = 1 ∨ k % 4 = 2 by omega) with h0 | h1 | h2'
  · -- k ≡ 0 : factor = 1 + A ≠ 0  (else A = -1)
    rw [h0, pow_zero]
    intro h; exact A_ne_neg_one hA h2 (by linear_combination h)
  · -- k ≡ 1 : factor = 2A ≠ 0
    rw [h1, pow_one]
    intro h
    have h2a : (2 : K) * A = 0 := by linear_combination h
    exact (mul_eq_zero.mp h2a).elim h2 hAne
  · -- k ≡ 2 : factor = A - 1 ≠ 0  (else A = 1)
    rw [h2', hA]
    intro h; exact A_ne_one hA h2 (by linear_combination h)

/-- The order-four character filter as an iff. -/
theorem filter_factor_zero_iff (hA : A ^ 2 = -1) (h2 : (2 : K) ≠ 0) (k : ℕ) :
    A ^ k + A = 0 ↔ k % 4 = 3 := by
  constructor
  · intro h; by_contra hk; exact filter_factor_ne_zero hA h2 hk h
  · exact filter_factor_eq_zero hA

/-! ## 3. Theorem "Order-Four Character Filter" (coefficient form)

The functional equation `P(AX) + A·P(X) = A·X`, read off coefficient by
coefficient, is exactly the hypothesis `hfe` below:
`(Aᵏ + A)·cₖ = A·[k=1]`.  This is the self-contained heart of the theorem. -/

/-- Every coefficient off the support `{1} ∪ {k ≡ 3 (mod 4)}` vanishes, and
    the linear coefficient equals `1/2`. -/
theorem order_four_filter_coeff (hA : A ^ 2 = -1) (h2 : (2 : K) ≠ 0)
    (c : ℕ → K) (hfe : ∀ k, (A ^ k + A) * c k = if k = 1 then A else 0) :
    (∀ k, k ≠ 1 → k % 4 ≠ 3 → c k = 0) ∧ c 1 = 2⁻¹ := by
  refine ⟨?_, ?_⟩
  · intro k hk1 hk3
    have h := hfe k
    rw [if_neg hk1] at h
    have hfac := filter_factor_ne_zero hA h2 hk3
    exact (mul_eq_zero.mp h).resolve_left hfac
  · have h := hfe 1
    rw [if_pos rfl] at h
    -- (A¹ + A)·c₁ = A, i.e. 2A·c₁ = A, and A ≠ 0
    have hAne : A ≠ 0 := A_ne_zero hA
    have hkey : A * (2 * c 1 - 1) = 0 := by
      have : A ^ 1 + A = 2 * A := by ring
      rw [this] at h; linear_combination h
    have h21 : 2 * c 1 - 1 = 0 := (mul_eq_zero.mp hkey).resolve_left hAne
    have hone : (2 : K) * c 1 = 1 := by linear_combination h21
    have hfin : c 1 = 1 / 2 := by rw [eq_div_iff h2]; linear_combination hone
    rwa [one_div] at hfin

/-! ## 4. Bridge to the genuine polynomial functional equation

We prove `(P.comp (C A · X)).coeff n = Aⁿ · P.coeff n`, so that the honest
polynomial statement `P(AX) + A·P = A·X` reduces to `order_four_filter_coeff`. -/

variable (A) in
/-- Substituting `X ↦ A·X` scales the `n`-th coefficient by `Aⁿ`. -/
lemma coeff_comp_C_mul_X (P : K[X]) (n : ℕ) :
    (P.comp (C A * X)).coeff n = A ^ n * P.coeff n := by
  induction P using Polynomial.induction_on' with
  | add p q hp hq =>
      simp only [add_comp, coeff_add, hp, hq, mul_add]
  | monomial m a =>
      have hpow : (C A * X) ^ m = C (A ^ m) * X ^ m := by rw [mul_pow, ← C_pow]
      rw [monomial_comp, hpow, coeff_C_mul, coeff_C_mul, coeff_X_pow, coeff_monomial]
      by_cases h : n = m
      · subst h; rw [if_pos rfl, if_pos rfl]; ring
      · rw [if_neg h, if_neg (fun e => h e.symm)]; ring

/-- Theorem "Order-Four Character Filter" (polynomial form).
    If `P` satisfies the order-four functional equation
    `P(AX) + A·P(X) = A·X`, then its monomial support lies in
    `{1} ∪ {k ≡ 3 (mod 4)}` and its linear coefficient is `1/2`. -/
theorem order_four_filter (hA : A ^ 2 = -1) (h2 : (2 : K) ≠ 0) (P : K[X])
    (hfe : P.comp (C A * X) + C A * P = C A * X) :
    (∀ k, k ≠ 1 → k % 4 ≠ 3 → P.coeff k = 0) ∧ P.coeff 1 = 2⁻¹ := by
  apply order_four_filter_coeff hA h2 P.coeff
  intro k
  have hco := congrArg (fun p : K[X] => p.coeff k) hfe
  simp only [coeff_add, coeff_comp_C_mul_X, coeff_C_mul] at hco
  -- hco : Aᵏ · cₖ + A · cₖ = A · (X.coeff k)
  by_cases hk : k = 1
  · subst hk
    rw [coeff_X_one] at hco
    rw [if_pos rfl]
    linear_combination hco
  · rw [if_neg hk]
    have hz : (X : K[X]).coeff k = 0 := by rw [coeff_X]; exact if_neg (by omega)
    rw [hz, mul_zero] at hco
    linear_combination hco

/-! ## 5. Corollaries: support and the factored evaluation form -/

/-- The monomial support of `P` is contained in `{1} ∪ {k : k ≡ 3 (mod 4)}`. -/
theorem support_subset (hA : A ^ 2 = -1) (h2 : (2 : K) ≠ 0) (P : K[X])
    (hfe : P.comp (C A * X) + C A * P = C A * X) :
    ∀ k ∈ P.support, k = 1 ∨ k % 4 = 3 := by
  intro k hk
  by_contra hcon
  rw [not_or] at hcon
  exact (Polynomial.mem_support_iff.mp hk)
    ((order_four_filter hA h2 P hfe).1 k hcon.1 hcon.2)

/-- (Theorem "Factored evaluation form", coefficient statement.)
    Every coefficient of `P` is pinned down:  the linear coefficient is `1/2`,
    the coefficients at exponents `≡ 3 (mod 4)` are the free parameters
    (these are the coefficients `q_j = c_{4j+3}` of `Q`), and every other
    coefficient vanishes.  This is exactly the decomposition
    `P = ½·X + X³·Q(X⁴)`. -/
theorem factored_form (hA : A ^ 2 = -1) (h2 : (2 : K) ≠ 0) (P : K[X])
    (hfe : P.comp (C A * X) + C A * P = C A * X) (k : ℕ) :
    P.coeff k = if k = 1 then 2⁻¹ else if k % 4 = 3 then P.coeff k else 0 := by
  obtain ⟨hsupp, hc1⟩ := order_four_filter hA h2 P hfe
  split_ifs with h h'
  · rw [h]; exact hc1
  · rfl
  · exact hsupp k h h'

end Field

end OrderFour
