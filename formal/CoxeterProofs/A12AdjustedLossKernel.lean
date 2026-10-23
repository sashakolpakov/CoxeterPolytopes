import Std

/-!
PREPARED SOURCE ONLY. No compilation or actual axiom audit is credited.

Continuous adjusted-loss and simultaneous-foot logic for the accepted
A12/actual-typed-seven component estimate. Primitive radical equations,
rather than a precomputed maximum, supply the full interval inequalities.

External: actual geometry and field/helper inventory, reference coordinates,
the individual coefficient/loss interfaces, all-wall feasible points, the
one-vector pressure/norm estimates, and defect-specific spectator floors.
No real square-root function, number field or polytope instance is supplied.
-/

namespace CoxeterProofs.A12AdjustedLossKernel

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- The circle identity proves the linear bound from a strict rational gap. -/
theorem circle_upper (a c r u b : R)
    (root_square : b ^ 2 = u * (1 - u))
    (r_nonnegative : 0 ≤ r) (r_gap : a ^ 2 + c ^ 2 < r ^ 2) :
    2 * a * u + 2 * c * b < a + r := by
  have complement : 0 ≤ (c * (2 * u - 1) - 2 * a * b) ^ 2 :=
    OrderedRing.sq_nonneg
  have identity :
      (a * (2 * u - 1) + 2 * c * b) ^ 2 +
        (c * (2 * u - 1) - 2 * a * b) ^ 2 = a ^ 2 + c ^ 2 := by
    grind only
  by_cases target : 2 * a * u + 2 * c * b < a + r
  · exact target
  · have first : 0 ≤ a * (2 * u - 1) + 2 * c * b - r := by grind
    have second : 0 ≤ a * (2 * u - 1) + 2 * c * b + r := by grind
    have product := OrderedRing.mul_nonneg first second
    grind

/-- The radical never exceeds one half; no sampled root bound is assumed. -/
theorem radical_half (u b : R)
    (b_nonnegative : 0 ≤ b) (root_square : b ^ 2 = u * (1 - u)) :
    2 * b ≤ 1 := by
  have square : 0 ≤ (2 * u - 1) ^ 2 := OrderedRing.sq_nonneg
  by_cases target : 2 * b ≤ 1
  · exact target
  · have product := OrderedRing.mul_pos
      (show 0 < 2 * b - 1 by grind) (show 0 < 2 * b + 1 by grind)
    grind

/-- The outer coefficient uses beta >= .08 on the entire closed interval. -/
theorem outer_adjusted_polynomial (u b : R)
    (u_lower : 0 ≤ u) (u_upper : u ≤ 1)
    (b_nonnegative : 0 ≤ b) (root_square : b ^ 2 = u * (1 - u)) :
    264 * u + 80 * u ^ 2 + 2000 * b < 1201 := by
  have half := radical_half u b b_nonnegative root_square
  have unit_product := OrderedRing.mul_nonneg u_lower (show 0 ≤ 1 - u by grind)
  by_cases left : 5 * u ≤ 2
  · grind
  · have linear := circle_upper (R := R) 212 1000 1023 u b
      root_square (by grind) (by grind)
    have product := OrderedRing.mul_nonneg
      (show 0 ≤ 5 * u - 2 by grind) (show 0 ≤ 8 - 5 * u by grind)
    have correction : 16 ≤ 25 * (2 * u - u ^ 2) := by grind
    grind

/-- Three analytic ranges prove the sharper inner coefficient bound. -/
theorem inner_adjusted_polynomial (u b : R)
    (u_lower : 0 ≤ u) (u_upper : u ≤ 1)
    (b_nonnegative : 0 ≤ b) (root_square : b ^ 2 = u * (1 - u)) :
    698 * u + 61 * u ^ 2 + 4000 * b < 2402 := by
  have half := radical_half u b b_nonnegative root_square
  have unit_product := OrderedRing.mul_nonneg u_lower (show 0 ≤ 1 - u by grind)
  by_cases left : 2 * u ≤ 1
  · have product := OrderedRing.mul_nonneg
      (show 0 ≤ 1 - 2 * u by grind) (show 0 ≤ 1 + 2 * u by grind)
    grind
  · by_cases middle : 10 * u ≤ 7
    · have chord_product := OrderedRing.mul_nonneg
        (show 0 ≤ 2 * u - 1 by grind) (show 0 ≤ 7 - 10 * u by grind)
      have chord : 20 * u ^ 2 ≤ 24 * u - 7 := by grind
      have linear := circle_upper (R := R) 1928 10000 10185 u b
        root_square (by grind) (by grind)
      grind
    · have positive_product := OrderedRing.mul_pos
        (show 0 < 2 * u - 1 by grind) (show 0 < 13 * u - 9 by grind)
      have chord : 5 * b < 3 - u := by
        by_cases target : 5 * b < 3 - u
        · exact target
        · have first : 0 ≤ 5 * b - (3 - u) := by grind
          have second : 0 ≤ 5 * b + (3 - u) := by grind
          have product := OrderedRing.mul_nonneg first second
          grind
      grind

/-- Actual individual losses inherit the bound with their own coefficient. -/
theorem individual_adjusted_loss (u b beta loss : R)
    (u_lower : 0 ≤ u) (u_upper : u ≤ 1)
    (b_nonnegative : 0 ≤ b) (root_square : b ^ 2 = u * (1 - u))
    (individual_role :
      (160 ≤ 2000 * beta ∧ 2000 * loss < 424 * u + 2000 * b) ∨
      (2000 * beta = 61 ∧ 4000 * loss < 820 * u + 4000 * b)) :
    2000 * loss - 1000 * beta * (2 * u - u ^ 2) < 1201 := by
  cases individual_role with
  | inl outer =>
    have scalar := outer_adjusted_polynomial u b u_lower u_upper
      b_nonnegative root_square
    have q_nonnegative : 0 ≤ 2 * u - u ^ 2 := by
      have product := OrderedRing.mul_nonneg u_lower (show 0 ≤ 2 - u by grind)
      grind
    have price := OrderedRing.mul_nonneg
      (show 0 ≤ 1000 * beta - 80 by grind) q_nonnegative
    grind
  | inr inner =>
    have scalar := inner_adjusted_polynomial u b u_lower u_upper
      b_nonnegative root_square
    have scaled_beta :
        2000 * beta * (2 * u - u ^ 2) = 61 * (2 * u - u ^ 2) := by
      rw [inner.1]
    grind

/-- Each foot carries primitive data, not an assumed final gain. -/
structure Foot (R : Type u) [CommRing R] [LE R] [LT R]
    [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R] where
  u : R
  b : R
  beta : R
  loss : R
  u_lower : 0 ≤ u
  u_upper : u ≤ 1
  b_nonnegative : 0 ≤ b
  root_square : b ^ 2 = u * (1 - u)
  beta_positive : 0 < beta
  individual_role :
    (160 ≤ 2000 * beta ∧ 2000 * loss < 424 * u + 2000 * b) ∨
    (2000 * beta = 61 ∧ 4000 * loss < 820 * u + 4000 * b)

def credit (f : Foot R) : R := 1000 * f.beta ^ 2 * (2 * f.u - f.u ^ 2)
def debit (f : Foot R) : R := f.beta * (2000 * f.loss - 1201)

def sumCredit : List (Foot R) → R
  | [] => 0
  | f :: fs => credit f + sumCredit fs

def sumDebit : List (Foot R) → R
  | [] => 0
  | f :: fs => debit f + sumDebit fs

/-- Finite lists include adjacent and arbitrarily many weak feet. -/
theorem finite_gain (fs : List (Foot R)) :
    0 ≤ sumCredit fs - sumDebit fs ∧
      (fs ≠ [] → 0 < sumCredit fs - sumDebit fs) := by
  induction fs with
  | nil =>
    simp only [sumCredit, sumDebit]
    constructor <;> grind
  | cons f fs ih =>
    have scalar := individual_adjusted_loss f.u f.b f.beta f.loss
      f.u_lower f.u_upper f.b_nonnegative f.root_square f.individual_role
    have gain := OrderedRing.mul_pos f.beta_positive
      (show 0 < 1201 - 2000 * f.loss + 1000 * f.beta * (2 * f.u - f.u ^ 2) by
        grind)
    have difference : 0 < credit f - debit f := by
      unfold credit debit
      grind
    simp only [sumCredit, sumDebit]
    constructor
    · grind
    · intro _
      grind

/-- One pressure, one norm and one completed square absorb the entire list. -/
theorem simultaneous_foot_completion
    (fs : List (Foot R)) (base pressure normZ energy : R)
    (one_pressure : 2000 * base - sumDebit fs ≤ 2000 * pressure)
    (one_norm : 1000 * normZ ≤ 1000 * base - sumCredit fs)
    (one_square : 0 ≤ energy - 2 * pressure + normZ)
    (empty_strict : fs = [] → base < pressure) :
    base < energy := by
  have gain := finite_gain fs
  by_cases empty : fs = []
  · have strict := empty_strict empty
    subst fs
    simp only [sumCredit, sumDebit] at one_norm one_pressure
    grind
  · have strict := gain.2 empty
    grind

/-- The rational reference energy is calculated, rather than assumed as a floor. -/
theorem reference_energy (h5 h6 base : R)
    (outer_height : 2000 * h5 = 3759)
    (inner_height : 2000 * h6 = 3881)
    (reference_norm : 5 * base = h5 ^ 2 + 5 * (h6 - h5) ^ 2) :
    20000000 * base = 14204501 := by
  grind

/-- The literal rational reference and the two defect-specific floors suffice. -/
theorem typed_seven_budget_contradiction (a defect spectator total : R)
    (component_bound : 14204501 < 20000000 * a)
    (actual_defect_floor : 9 < 100 * defect)
    (conditional_spectator_floor : 1 < 5 * spectator)
    (one_orthogonal_sum : a + defect + spectator ≤ total)
    (original_norm : total < 1) : False := by
  grind

/-- Complete scalar chain from the reference and individual feet to contradiction. -/
theorem complete_typed_seven_scalar_contradiction
    (fs : List (Foot R)) (h5 h6 base pressure normZ a defect spectator total : R)
    (outer_height : 2000 * h5 = 3759)
    (inner_height : 2000 * h6 = 3881)
    (reference_norm : 5 * base = h5 ^ 2 + 5 * (h6 - h5) ^ 2)
    (one_pressure : 2000 * base - sumDebit fs ≤ 2000 * pressure)
    (one_norm : 1000 * normZ ≤ 1000 * base - sumCredit fs)
    (one_square : 0 ≤ a - 2 * pressure + normZ)
    (empty_strict : fs = [] → base < pressure)
    (actual_defect_floor : 9 < 100 * defect)
    (conditional_spectator_floor : 1 < 5 * spectator)
    (one_orthogonal_sum : a + defect + spectator ≤ total)
    (original_norm : total < 1) : False := by
  have component := simultaneous_foot_completion fs base pressure normZ a
    one_pressure one_norm one_square empty_strict
  have reference := reference_energy h5 h6 base outer_height inner_height reference_norm
  have component_bound : 14204501 < 20000000 * a := by grind
  exact typed_seven_budget_contradiction a defect spectator total component_bound
    actual_defect_floor conditional_spectator_floor one_orthogonal_sum original_norm

end CoxeterProofs.A12AdjustedLossKernel

#print axioms CoxeterProofs.A12AdjustedLossKernel.circle_upper
#print axioms CoxeterProofs.A12AdjustedLossKernel.radical_half
#print axioms CoxeterProofs.A12AdjustedLossKernel.outer_adjusted_polynomial
#print axioms CoxeterProofs.A12AdjustedLossKernel.inner_adjusted_polynomial
#print axioms CoxeterProofs.A12AdjustedLossKernel.individual_adjusted_loss
#print axioms CoxeterProofs.A12AdjustedLossKernel.finite_gain
#print axioms CoxeterProofs.A12AdjustedLossKernel.simultaneous_foot_completion
#print axioms CoxeterProofs.A12AdjustedLossKernel.reference_energy
#print axioms CoxeterProofs.A12AdjustedLossKernel.typed_seven_budget_contradiction
#print axioms CoxeterProofs.A12AdjustedLossKernel.complete_typed_seven_scalar_contradiction
