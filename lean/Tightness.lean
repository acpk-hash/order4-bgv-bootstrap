/-
  Tightness.lean

  Formal verification (Lean 4 + Mathlib) of the quantitative "tightness"
  results of

    "Order-Four Symmetry in BGV Bootstrapping:
     Faster Digit Extraction for Large Primes"

  complementing `OrderFour.lean` (the Order-Four Character Filter) and
  `Character.lean` (the general order-`r` framework).  This file
  formalizes, with *no* `sorry`:

  * `card_v3`, `card_v3_mod`
        (Theorem "Tight term count", counting part)
        Among the exponents `k < |S_A| = 4q + 1` exactly
        `q = (|S_A| - 1)/4` satisfy `k ≡ 3 (mod 4)` — namely
        `3, 7, …, 4q - 1`.  (`|S_A| = (2B+1)²` is an odd square, hence
        `≡ 1 (mod 4)`, so this is exactly the paper's case.)

  * `support_card_le`, `term_count_bound`
        (Theorem "Tight term count", bound `|T| ≤ (|S_A| - 1)/4 + 1`)
        A coefficient sequence supported in `{1} ∪ {k ≡ 3 (mod 4)}` — the
        conclusion of the Order-Four Character Filter — has at most
        `q + 1` nonzero coefficients below `4q + 1`.  Chained with
        `OrderFour.order_four_filter_coeff`, this machine-checks the
        paper's term-count bound directly from the functional equation.

  * `encode_abs_le`, `encode_injective_int`
        (Parameter Admissibility Theorem, injectivity core)
        The digit encoding `(η, λ) ↦ ηA + λ` with `|η|, |λ| ≤ B` has
        absolute value at most `|A|·B + B`; hence if
        `2(|A|·B + B) < p`, two encodings that agree modulo `p` are
        already equal over `ℤ`.  This is the sufficiency direction of
        the admissibility condition (Definition 1).

  * `ratio_eq`, `ratio_lt_two`, `ratio_tendsto`
        (Complexity theorem, asymptotic part)
        The evaluation-cost ratio `2√d / (2√(d/4) + 4)` equals
        `2√d / (√d + 4)`, is `< 2` for every `d ≥ 0`, and tends to `2`
        as `d → ∞`: the order-four method asymptotically halves the
        number of non-constant multiplications, approaching the
        speedup factor `2` from below.
-/
import Mathlib
import BootstrapAlgebra.OrderFour

namespace OrderFour.Tightness

/-! ## 1. Theorem "Tight term count" — the counting part

The paper's set of exponents is `{0, 1, …, |S_A| - 1}` with
`|S_A| = (2B+1)² = 4q + 1` (an odd square is `≡ 1 (mod 4)`).  The filter
passes the residue class `k ≡ 3 (mod 4)`, i.e. `3, 7, …, 4q - 1`:
exactly `q` exponents. -/

/-- Among the exponents `k < 4q + 1` exactly `q` satisfy `k ≡ 3 (mod 4)`
    (namely `3, 7, …, 4q - 1`).  Counting core of the tight-term-count
    theorem, with `|S_A| = (2B+1)² = 4q + 1`. -/
theorem card_v3 (q : ℕ) :
    ((Finset.range (4 * q + 1)).filter (fun k => k % 4 = 3)).card = q := by
  induction q with
  | zero => decide
  | succ q ih =>
      -- Going from `4q + 1` to `4(q+1) + 1` adds the exponents
      -- `4q+1, 4q+2, 4q+3, 4q+4`, of which only `4q+3 ≡ 3 (mod 4)`.
      have hstep :
          (Finset.range (4 * (q + 1) + 1)).filter (fun k => k % 4 = 3)
            = insert (4 * q + 3)
                ((Finset.range (4 * q + 1)).filter (fun k => k % 4 = 3)) := by
        ext k
        simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_insert]
        omega
      have hnot : 4 * q + 3 ∉
          (Finset.range (4 * q + 1)).filter (fun k => k % 4 = 3) := by
        simp only [Finset.mem_filter, Finset.mem_range]
        omega
      rw [hstep, Finset.card_insert_of_notMem hnot, ih]

/-- General-`n` form: for `n ≡ 1 (mod 4)` (in the paper
    `n = |S_A| = (2B+1)²`, an odd square), the number of exponents
    `k < n` with `k ≡ 3 (mod 4)` is `(n - 1)/4`. -/
theorem card_v3_mod {n : ℕ} (hn : n % 4 = 1) :
    ((Finset.range n).filter (fun k => k % 4 = 3)).card = (n - 1) / 4 := by
  obtain ⟨q, rfl⟩ : ∃ q, n = 4 * q + 1 := ⟨n / 4, by omega⟩
  rw [card_v3]
  omega

/-- (Theorem "Tight term count", upper bound.)  If a coefficient
    sequence `c` vanishes off `{1} ∪ {k ≡ 3 (mod 4)}` — precisely the
    conclusion of `OrderFour.order_four_filter_coeff` — then it has at
    most `q + 1 = (|S_A| - 1)/4 + 1` nonzero entries below
    `|S_A| = 4q + 1`. -/
theorem support_card_le {K : Type*} [Field K] [DecidableEq K]
    (c : ℕ → K) (q : ℕ) (hc : ∀ k, k ≠ 1 → k % 4 ≠ 3 → c k = 0) :
    ((Finset.range (4 * q + 1)).filter (fun k => c k ≠ 0)).card ≤ q + 1 := by
  have hsub : (Finset.range (4 * q + 1)).filter (fun k => c k ≠ 0) ⊆
      insert 1 ((Finset.range (4 * q + 1)).filter (fun k => k % 4 = 3)) := by
    intro k hk
    rw [Finset.mem_filter, Finset.mem_range] at hk
    rw [Finset.mem_insert, Finset.mem_filter, Finset.mem_range]
    by_cases h1 : k = 1
    · exact Or.inl h1
    · by_cases h3 : k % 4 = 3
      · exact Or.inr ⟨hk.1, h3⟩
      · exact absurd (hc k h1 h3) hk.2
  calc ((Finset.range (4 * q + 1)).filter (fun k => c k ≠ 0)).card
      ≤ (insert 1 ((Finset.range (4 * q + 1)).filter
            (fun k => k % 4 = 3))).card := Finset.card_le_card hsub
    _ ≤ ((Finset.range (4 * q + 1)).filter
            (fun k => k % 4 = 3)).card + 1 := Finset.card_insert_le _ _
    _ = q + 1 := by rw [card_v3]

/-- (Theorem "Tight term count", end to end.)  Directly from the
    order-four functional equation, read coefficient-wise as
    `(Aᵏ + A)·cₖ = A·[k = 1]`: the digit-extraction polynomial has at
    most `(|S_A| - 1)/4 + 1` nonzero coefficients below
    `|S_A| = 4q + 1`.  This chains the Order-Four Character Filter with
    the counting bound above. -/
theorem term_count_bound {K : Type*} [Field K] [DecidableEq K] {A : K}
    (hA : A ^ 2 = -1) (h2 : (2 : K) ≠ 0) (c : ℕ → K)
    (hfe : ∀ k, (A ^ k + A) * c k = if k = 1 then A else 0) (q : ℕ) :
    ((Finset.range (4 * q + 1)).filter (fun k => c k ≠ 0)).card ≤ q + 1 :=
  support_card_le c q (order_four_filter_coeff hA h2 c hfe).1

/-! ## 2. Parameter Admissibility — the injectivity core

Definition 1 requires the digit encoding `(η, λ) ↦ ηA + λ`,
`|η|, |λ| ≤ B`, to be injective modulo `p`.  The checkable arithmetic
core: every encoding lies in `[-(|A|B + B), |A|B + B]`, so as soon as
`2(|A|B + B) < p` two encodings that agree mod `p` differ by a multiple
of `p` of absolute value `< p`, hence are equal over `ℤ`. -/

/-- Size bound on the digit encoding: if `|η|, |λ| ≤ B` then
    `|ηA + λ| ≤ |A|·B + B`. -/
theorem encode_abs_le (A B η l : ℤ) (hη : |η| ≤ B) (hl : |l| ≤ B) :
    |η * A + l| ≤ |A| * B + B := by
  calc |η * A + l| ≤ |η * A| + |l| := abs_add_le _ _
    _ = |η| * |A| + |l| := by rw [abs_mul]
    _ ≤ B * |A| + B :=
        add_le_add (mul_le_mul_of_nonneg_right hη (abs_nonneg A)) hl
    _ = |A| * B + B := by ring

/-- An integer multiple of `p` of absolute value `< p` is zero. -/
private lemma eq_zero_of_dvd_of_abs_lt {p x : ℤ} (hdvd : p ∣ x)
    (habs : |x| < p) : x = 0 := by
  rcases hdvd with ⟨t, rfl⟩
  have hp0 : 0 < p := lt_of_le_of_lt (abs_nonneg _) habs
  by_contra hne
  have ht0 : t ≠ 0 := by
    rintro rfl
    exact hne (mul_zero p)
  have h1 : (1 : ℤ) ≤ |t| := Int.one_le_abs ht0
  have hle : p ≤ |p * t| := by
    rw [abs_mul, abs_of_pos hp0]
    calc p = p * 1 := (mul_one p).symm
      _ ≤ p * |t| := mul_le_mul_of_nonneg_left h1 hp0.le
  linarith

/-- (Parameter Admissibility, sufficiency of `2(|A|B + B) < p`.)
    Two digit encodings `η₁A + λ₁` and `η₂A + λ₂` with digits in
    `[-B, B]` that agree modulo `p` (their difference is a multiple of
    `p`) are equal over `ℤ`.  In particular the encoding
    `(η, λ) ↦ ηA + λ mod p` is injective on `[-B,B]²`. -/
theorem encode_injective_int (A B p : ℤ) (hp : 2 * (|A| * B + B) < p)
    {η₁ l₁ η₂ l₂ : ℤ}
    (hη₁ : |η₁| ≤ B) (hl₁ : |l₁| ≤ B) (hη₂ : |η₂| ≤ B) (hl₂ : |l₂| ≤ B)
    (hdvd : p ∣ (η₁ * A + l₁) - (η₂ * A + l₂)) :
    η₁ * A + l₁ = η₂ * A + l₂ := by
  have h1 := encode_abs_le A B η₁ l₁ hη₁ hl₁
  have h2 := encode_abs_le A B η₂ l₂ hη₂ hl₂
  have habs : |(η₁ * A + l₁) - (η₂ * A + l₂)| < p := by
    calc |(η₁ * A + l₁) - (η₂ * A + l₂)|
        = |(η₁ * A + l₁) + -(η₂ * A + l₂)| := by rw [sub_eq_add_neg]
      _ ≤ |η₁ * A + l₁| + |-(η₂ * A + l₂)| := abs_add_le _ _
      _ = |η₁ * A + l₁| + |η₂ * A + l₂| := by rw [abs_neg]
      _ ≤ (|A| * B + B) + (|A| * B + B) := add_le_add h1 h2
      _ = 2 * (|A| * B + B) := by ring
      _ < p := hp
  have hzero := eq_zero_of_dvd_of_abs_lt hdvd habs
  linarith

/-! ## 3. Complexity theorem — the asymptotic evaluation-cost ratio

The baseline digit extraction costs `2√d` non-constant multiplications;
the order-four method costs `2√(d/4) + 4`.  Since
`2√(d/4) + 4 = √d + 4`, the cost ratio is `2√d/(√d + 4) = 2/(1 + 4/√d)`:
it is always below `2` and tends to `2` as the degree `d → ∞`. -/

/-- The evaluation-cost ratio of the paper's complexity theorem:
    baseline `2√d` non-constant multiplications versus `2√(d/4) + 4`
    for the order-four method. -/
noncomputable def ratio (d : ℝ) : ℝ :=
  (2 * Real.sqrt d) / (2 * Real.sqrt (d / 4) + 4)

private lemma sqrt_div_four {d : ℝ} (hd : 0 ≤ d) :
    Real.sqrt (d / 4) = Real.sqrt d / 2 := by
  have h : d / 4 = (Real.sqrt d / 2) ^ 2 := by
    rw [div_pow, Real.sq_sqrt hd]
    norm_num
  rw [h]
  exact Real.sqrt_sq (by positivity)

/-- Algebraic form of the cost ratio: for `d ≥ 0`,
    `2√d / (2√(d/4) + 4) = 2√d / (√d + 4)`. -/
theorem ratio_eq {d : ℝ} (hd : 0 ≤ d) :
    ratio d = 2 * Real.sqrt d / (Real.sqrt d + 4) := by
  unfold ratio
  rw [sqrt_div_four hd]
  ring

/-- The speedup never exceeds the asymptotic factor: `ratio d < 2` for
    every degree `d ≥ 0` (the limit `2` is approached from below). -/
theorem ratio_lt_two {d : ℝ} (hd : 0 ≤ d) : ratio d < 2 := by
  rw [ratio_eq hd]
  have h1 : (0 : ℝ) ≤ Real.sqrt d := Real.sqrt_nonneg d
  have h2 : (0 : ℝ) < Real.sqrt d + 4 := by linarith
  have key : 2 * Real.sqrt d / (Real.sqrt d + 4) - 2
      = -8 / (Real.sqrt d + 4) := by
    rw [eq_div_iff h2.ne', sub_mul, div_mul_cancel₀ _ h2.ne']
    ring
  have hneg : -8 / (Real.sqrt d + 4) < 0 :=
    div_neg_of_neg_of_pos (by norm_num) h2
  linarith

private lemma tendsto_sqrt_atTop' :
    Filter.Tendsto Real.sqrt Filter.atTop Filter.atTop := by
  apply Filter.tendsto_atTop.mpr
  intro b
  filter_upwards [Filter.eventually_ge_atTop ((max 0 b) ^ 2)] with x hx
  have hb : (0 : ℝ) ≤ max 0 b := le_max_left 0 b
  have h1 : Real.sqrt ((max 0 b) ^ 2) ≤ Real.sqrt x := Real.sqrt_le_sqrt hx
  rw [Real.sqrt_sq hb] at h1
  exact le_trans (le_max_right 0 b) h1

private lemma two_mul_div_add_four_tendsto :
    Filter.Tendsto (fun x : ℝ => 2 * x / (x + 4))
      Filter.atTop (nhds 2) := by
  have h1 : Filter.Tendsto (fun x : ℝ => x + 4)
      Filter.atTop Filter.atTop := by
    apply Filter.tendsto_atTop.mpr
    intro b
    filter_upwards [Filter.eventually_ge_atTop (b - 4)] with x hx
    linarith
  have h2 : Filter.Tendsto (fun x : ℝ => 8 / (x + 4))
      Filter.atTop (nhds 0) := Filter.Tendsto.div_atTop tendsto_const_nhds h1
  have h3 : Filter.Tendsto (fun x : ℝ => 2 - 8 / (x + 4))
      Filter.atTop (nhds (2 - 0)) := Filter.Tendsto.sub tendsto_const_nhds h2
  rw [sub_zero] at h3
  apply h3.congr'
  filter_upwards [Filter.eventually_ge_atTop (1 : ℝ)] with x hx
  have hx4 : x + 4 ≠ 0 := by linarith
  rw [eq_div_iff hx4, sub_mul, div_mul_cancel₀ _ hx4]
  ring

/-- (Complexity theorem, asymptotic part.)  The evaluation-cost ratio
    `2√d / (2√(d/4) + 4)` tends to `2` as the degree `d → ∞`: the
    order-four method asymptotically halves the number of non-constant
    multiplications. -/
theorem ratio_tendsto :
    Filter.Tendsto ratio Filter.atTop (nhds 2) := by
  have hcomp : Filter.Tendsto ((fun x : ℝ => 2 * x / (x + 4)) ∘ Real.sqrt)
      Filter.atTop (nhds 2) :=
    two_mul_div_add_four_tendsto.comp tendsto_sqrt_atTop'
  apply hcomp.congr'
  filter_upwards [Filter.eventually_ge_atTop (0 : ℝ)] with d hd
  simp only [Function.comp_apply]
  exact (ratio_eq hd).symm

end OrderFour.Tightness
