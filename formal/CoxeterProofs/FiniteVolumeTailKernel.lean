import Std
import Init.Grind.Ordered

/-!
Scalar higher-tail logic in finite_volume57.tex.
Compilation and axiom-audit execution are recorded separately from this source.

The finite table below defines exact eightfold finite caps and shape
coefficients for every dimension 59..163, and checks their comparisons
by kernel reduction. Actual rules, cusp envelopes, concentration and table
coverage remain application inputs. The ordered-ring part closes the
narrow-wall scalar contradiction from the cleared Lorentz identity and
the near/far deficit estimates. No dimension ceiling is an input field.
-/

namespace CoxeterProofs.FiniteVolumeTailKernel

/-- Eight times the stated finite cap, meaningful here on 59..163. -/
def u8 (n : Nat) : Nat :=
  match n with
  | 59 => 3132 | 60 => 3232 | 61 => 3360 | 62 => 3472
  | 63 => 3592 | 64 => 3712 | 65 => 3832 | 66 => 3960
  | 67 => 4032 | 68 => 4096 | 69 => 4352
  | n =>
    if n ≤ 77 then 64*n else
    if n ≤ 116 then 72*n else
    if n ≤ 126 then 80*n else
    if n ≤ 136 then 88*n else
    if n ≤ 146 then 96*n else
    if n ≤ 156 then 104*n else 112*n

/-- Eight times beta_n; natural subtraction is exact in the admitted range. -/
def beta8 (n : Nat) : Nat :=
  if n % 2 = 0 then n*(n-6) else (n-1)*(n-5)

/-- All 105 comparisons are computed, rather than supplied as a premise.
`decide` uses kernel reduction, not native_decide or an external checker. -/
theorem bounded_cap_table_check :
    ∀ i : Fin 105, u8 (59+i.val) ≤ beta8 (59+i.val) := by
  decide

theorem bounded_cap_comparison (n : Nat) (lower : 59 ≤ n) (upper : n ≤ 163) :
    u8 n ≤ beta8 n := by
  let i : Fin 105 := ⟨n-59, by omega⟩
  have index : 59+i.val = n := by dsimp [i]; omega
  have checked := bounded_cap_table_check i
  simpa only [index] using checked

/-- A small sum lemma keeps finite-vertex multiplicity explicit. -/
theorem sum_le_cap (xs : List Nat) (cap : Nat)
    (bound : ∀ x ∈ xs, x ≤ cap) : xs.sum ≤ cap*xs.length := by
  induction xs with
  | nil => simp
  | cons x xs ih =>
    have head := bound x (by simp)
    have tail := ih (by intro y hy; exact bound y (by simp [hy]))
    simp only [List.sum_cons, List.length_cons, Nat.mul_add, Nat.mul_one]
    omega

/-- Each pair stores (actual cusp cost, curvature), with a common scale. -/
def cuspCosts (xs : List (Nat × Nat)) : Nat := (xs.map Prod.fst).sum
def cuspCurvatures (xs : List (Nat × Nat)) : Nat := (xs.map Prod.snd).sum

theorem cusp_sum_comparison (xs : List (Nat × Nat))
    (bound : ∀ x ∈ xs, x.1 ≤ x.2) : cuspCosts xs ≤ cuspCurvatures xs := by
  induction xs with
  | nil => simp [cuspCosts, cuspCurvatures]
  | cons x xs ih =>
    have head := bound x (by simp)
    have tail := ih (by intro y hy; exact bound y (by simp [hy]))
    simp only [cuspCosts, cuspCurvatures, List.map_cons, List.sum_cons] at *
    omega

/-- Common payment subtraction, including an empty cusp list. All values
are nonnegative integers after one common denominator is cleared. -/
theorem payment_contradiction (finite : List Nat) (cusps : List (Nat × Nat))
    (deficit correction U beta : Nat)
    (finite_cap : ∀ x ∈ finite, x ≤ U)
    (cusp_cap : ∀ x ∈ cusps, x.1 ≤ x.2)
    (shape : beta*finite.length+cuspCurvatures cusps+correction ≤ deficit)
    (payment : deficit ≤ finite.sum+cuspCosts cusps)
    (positive : 0 < correction) (comparison : U ≤ beta) : False := by
  have finite_sum := sum_le_cap finite beta
    (by intro x hx; exact Nat.le_trans (finite_cap x hx) comparison)
  have cusp_sum := cusp_sum_comparison cusps cusp_cap
  omega

/-- In the bounded tail even U<=beta is now derived from the literal table. -/
theorem bounded_tail_payment_contradiction (n : Nat)
    (lower : 59 ≤ n) (upper : n ≤ 163)
    (scale : Nat)
    (finite : List Nat) (cusps : List (Nat × Nat)) (deficit correction : Nat)
    (finite_cap : ∀ x ∈ finite, x ≤ scale*u8 n)
    (cusp_cap : ∀ x ∈ cusps, x.1 ≤ x.2)
    (shape : (scale*beta8 n)*finite.length+cuspCurvatures cusps+correction ≤ deficit)
    (payment : deficit ≤ finite.sum+cuspCosts cusps)
    (positive : 0 < correction) : False := by
  exact payment_contradiction finite cusps deficit correction (scale*u8 n) (scale*beta8 n)
    finite_cap cusp_cap shape payment positive
    (Nat.mul_le_mul_left scale (bounded_cap_comparison n lower upper))

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- Scalar Lorentz closure in the near range. T is the sum of omission
targets, W the sum of squared displacements, H=q+q^{-1}, p=-(x,y).
The nontrivial all-factor loss<=4 is an explicit geometric/table input. -/
theorem near_pairing_bound (N M T W H p : R)
    (nonnegative : 0 ≤ M) (target : 9*T = 2*N-2)
    (loss : T-W ≤ 4) (lorentz : 2*p = M*W-2*H) (height : H ≤ 2*M) :
    M*(N-37) ≤ 9*p := by
  have onePositive : (0 : R) < 1 := OrderedRing.zero_lt_one
  have deficit : 2*N ≤ 9*W+38 := by grind
  have scaled := OrderedRing.mul_le_mul_of_nonneg_left deficit nonnegative
  have heightScaled := OrderedRing.mul_le_mul_of_nonneg_left height
    (show (0 : R) ≤ 18 by grind)
  grind

/-- The far loss is derived using both endpoint budgets separately;
the all-factor coefficient-four-thirds estimate remains an input. -/
theorem far_pairing_bound (N M T W H p C0 C1 : R)
    (nonnegative : 0 ≤ M) (target : 9*T = 2*N-2)
    (budget0 : C0 ≤ 2) (budget1 : C1 ≤ 2)
    (loss : 3*(T-W) ≤ 4*(C0+C1))
    (lorentz : 2*p = M*W-2*H) (height : H ≤ 2*M) :
    M*(N-43) ≤ 9*p := by
  have onePositive : (0 : R) < 1 := OrderedRing.zero_lt_one
  have deficit : 2*N ≤ 9*W+50 := by grind
  have scaled := OrderedRing.mul_le_mul_of_nonneg_left deficit nonnegative
  have heightScaled := OrderedRing.mul_le_mul_of_nonneg_left height
    (show (0 : R) ≤ 18 by grind)
  grind

/-- N>=163 and M>=1 force p>=14 from the near estimate. -/
theorem near_narrow_contradiction (N M T W H p : R)
    (rank : 163 ≤ N) (ratio : 1 ≤ M) (target : 9*T = 2*N-2)
    (loss : T-W ≤ 4) (lorentz : 2*p = M*W-2*H)
    (height : H ≤ 2*M) (narrow : p < 14) : False := by
  have onePositive : (0 : R) < 1 := OrderedRing.zero_lt_one
  have nonnegative : 0 ≤ M := by grind
  have pairing := near_pairing_bound N M T W H p
    nonnegative target loss lorentz height
  have product := OrderedRing.mul_nonneg
    (show 0 ≤ M-1 by grind) (show 0 ≤ N-37 by grind)
  grind

/-- The farther range needs 20M>=21: at N=163 it reaches exactly14. -/
theorem far_narrow_contradiction (N M T W H p C0 C1 : R)
    (rank : 163 ≤ N) (ratio : 21 ≤ 20*M) (target : 9*T = 2*N-2)
    (budget0 : C0 ≤ 2) (budget1 : C1 ≤ 2)
    (loss : 3*(T-W) ≤ 4*(C0+C1))
    (lorentz : 2*p = M*W-2*H) (height : H ≤ 2*M)
    (narrow : p < 14) : False := by
  have onePositive : (0 : R) < 1 := OrderedRing.zero_lt_one
  have nonnegative : 0 ≤ M := by grind
  have pairing := far_pairing_bound N M T W H p C0 C1
    nonnegative target budget0 budget1 loss lorentz height
  have product := OrderedRing.mul_nonneg nonnegative (show 0 ≤ N-163 by grind)
  have scaled : 126 ≤ 120*M := by grind
  grind

/-- The near/far case split is exhaustive for the common height ratio.
Framing-wall existence, N=n-1, and the narrowness theorem stay external. -/
theorem narrow_cusp_contradiction (N M T W H p C0 C1 : R)
    (rank : 163 ≤ N) (ratio : 1 ≤ M) (target : 9*T = 2*N-2)
    (budget0 : C0 ≤ 2) (budget1 : C1 ≤ 2)
    (near_loss : 20*M ≤ 21 → T-W ≤ 4)
    (far_loss : 21 ≤ 20*M → 3*(T-W) ≤ 4*(C0+C1))
    (lorentz : 2*p = M*W-2*H) (height : H ≤ 2*M)
    (narrow : p < 14) : False := by
  by_cases near : 20*M ≤ 21
  · exact near_narrow_contradiction N M T W H p rank ratio target
      (near_loss near) lorentz height narrow
  · have far : 21 ≤ 20*M := by grind
    exact far_narrow_contradiction N M T W H p C0 C1 rank far target
      budget0 budget1 (far_loss far) lorentz height narrow

end CoxeterProofs.FiniteVolumeTailKernel

#print axioms CoxeterProofs.FiniteVolumeTailKernel.bounded_cap_table_check
#print axioms CoxeterProofs.FiniteVolumeTailKernel.bounded_cap_comparison
#print axioms CoxeterProofs.FiniteVolumeTailKernel.payment_contradiction
#print axioms CoxeterProofs.FiniteVolumeTailKernel.bounded_tail_payment_contradiction
#print axioms CoxeterProofs.FiniteVolumeTailKernel.near_pairing_bound
#print axioms CoxeterProofs.FiniteVolumeTailKernel.far_pairing_bound
#print axioms CoxeterProofs.FiniteVolumeTailKernel.near_narrow_contradiction
#print axioms CoxeterProofs.FiniteVolumeTailKernel.far_narrow_contradiction
#print axioms CoxeterProofs.FiniteVolumeTailKernel.narrow_cusp_contradiction
