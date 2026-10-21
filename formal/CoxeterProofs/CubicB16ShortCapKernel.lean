import Std

/-!
Prepared source only: uncompiled, with zero actual axiom-audit outputs.

Internal integer deductions and ordered-ring energy algebra for the written
cubic B16+actual-I obstruction. Geometry, field embeddings/integral basis,
root brackets, the unit/J identities and the paid-reference inequality are
explicit application inputs; no polytope or field instance is constructed.
The Int variables below are integral coefficients, not substitutes for real
geometric coordinates. Only Std is imported.
-/

namespace CoxeterProofs.CubicB16ShortCapKernel

/-- The complete finite trace problem, after deriving the two field inequalities. -/
theorem cubic_coefficients_unique (S T : Int)
    (trace_lower : 17 ≤ S) (trace_upper : S ≤ 28)
    (coefficient_nonnegative : 0 ≤ T)
    (quarter_test : 4 * T ≤ S + 15) (root_test : 4 * S ≤ T * T) :
    S = 25 ∧ T = 10 := by
  have coefficient_cases :
      T = 0 ∨ T = 1 ∨ T = 2 ∨ T = 3 ∨ T = 4 ∨ T = 5 ∨
      T = 6 ∨ T = 7 ∨ T = 8 ∨ T = 9 ∨ T = 10 := by omega
  rcases coefficient_cases with rfl | rfl | rfl | rfl | rfl | rfl |
    rfl | rfl | rfl | rfl | rfl
  all_goals omega

/--
The integral coefficient strip derived from the two nonphysical embeddings.
The bounds are cleared inequalities, not a supplied coefficient conclusion.
-/
theorem i_strip_coefficients_unique (m n : Int)
    (coefficient_lower : -6 ≤ n) (coefficient_upper : n ≤ 6)
    (nonnegative_strip : 0 ≤ n →
      2205 * n < 10000 * m ∧ 10000 * m < 1845 * n + 2470)
    (negative_strip : n < 0 →
      1845 * n < 10000 * m ∧ 10000 * m < 2205 * n + 2470) :
    m = 0 ∧ n = -1 := by
  by_cases sign : 0 ≤ n
  · obtain ⟨lower, upper⟩ := nonnegative_strip sign
    have index_cases : n = 0 ∨ n = 1 ∨ n = 2 ∨ n = 3 ∨ n = 4 ∨ n = 5 ∨ n = 6 := by
      omega
    rcases index_cases with rfl | rfl | rfl | rfl | rfl | rfl | rfl
    all_goals omega
  · have negative : n < 0 := by omega
    obtain ⟨lower, upper⟩ := negative_strip negative
    have index_cases : n = -6 ∨ n = -5 ∨ n = -4 ∨ n = -3 ∨ n = -2 ∨ n = -1 := by
      omega
    rcases index_cases with rfl | rfl | rfl | rfl | rfl | rfl
    all_goals omega

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- The actual J=6-1/B identity forces the strict I-height lower bound. -/
theorem i_height_lower_from_unit (B eta h : R)
    (physical_root : 24 < B) (inverse_positive : 0 < eta)
    (inverse_identity : B * eta = 1) (height_identity : 4 * h = 2 - eta) :
    47 < 96 * h := by
  have positive := OrderedRing.mul_pos
    (show 0 < B - 24 by grind) inverse_positive
  have inverse_bound : 24 * eta < 1 := by grind
  grind

/-- The half-tangent relation yields r>.099, with no monotonicity axiom. -/
theorem i_energy_lower_from_height (h r : R)
    (height_lower : 47 < 96 * h)
    (energy_nonnegative : 0 ≤ r) (energy_below_one : r < 1)
    (height_relation : h * (1 - r) ^ 2 = 4 * r) :
    99 < 1000 * r := by
  have positive_square : 0 < (1 - r) ^ 2 := by
    have square : (1 - r) ^ 2 = (1 - r) * (1 - r) := by grind only
    rw [square]
    exact OrderedRing.mul_pos (by grind) (by grind)
  have positive := OrderedRing.mul_pos
    (show 0 < 96 * h - 47 by grind) positive_square
  have expansion : (96 * h - 47) * (1 - r) ^ 2 =
      96 * (h * (1 - r) ^ 2) - 47 * (1 - r) ^ 2 := by grind only
  rw [expansion, height_relation] at positive
  have polynomial_negative : 47 * r ^ 2 - 478 * r + 47 < 0 := by grind
  by_cases target : 99 < 1000 * r
  · exact target
  · have factor_one : 0 ≤ 99 - 1000 * r := by grind
    have factor_two : 0 ≤ 478000 - 47 * (1000 * r + 99) := by grind
    have product_nonnegative := OrderedRing.mul_nonneg factor_one factor_two
    have completion :
        (99 - 1000 * r) * (478000 - 47 * (1000 * r + 99)) =
          1000000 * (47 * r ^ 2 - 478 * r + 47) - 138647 := by grind only
    rw [completion] at product_nonnegative
    grind

/-- One paid B reference plus one I energy gives the exact total margin. -/
theorem paid_reference_norm_margin (base f S eB eI total : R)
    (reference_floor : 89026 ≤ 100000 * base)
    (endpoint_upper : 1000 * f ≤ 3311) (short_pressure : 7 < 2 * S)
    (paid_reference : 5000 * base + 1300 * (S - f) - 169 ≤ 5000 * eB)
    (i_energy : 99 < 1000 * eI) (original_total : total = eB + eI) :
    100460 < 100000 * total ∧ 1 < total := by
  have exact_margin : 100460 < 100000 * total := by grind
  grind

/--
Conditional actual-cubic contradiction after the explicit field/geometry bridge.
Neither an I energy margin nor the desired total bound is a premise.
-/
theorem cubic_short_cap_contradiction (B eta h r base f S eB eI total : R)
    (physical_root : 24 < B) (inverse_positive : 0 < eta)
    (inverse_identity : B * eta = 1) (height_identity : 4 * h = 2 - eta)
    (energy_nonnegative : 0 ≤ r) (energy_below_one : r < 1)
    (height_relation : h * (1 - r) ^ 2 = 4 * r) (actual_i_lower : r ≤ eI)
    (reference_floor : 89026 ≤ 100000 * base)
    (endpoint_upper : 1000 * f ≤ 3311) (short_pressure : 7 < 2 * S)
    (paid_reference : 5000 * base + 1300 * (S - f) - 169 ≤ 5000 * eB)
    (original_total : total = eB + eI) (norm_below_one : total < 1) : False := by
  have height := i_height_lower_from_unit B eta h physical_root inverse_positive
    inverse_identity height_identity
  have r_lower := i_energy_lower_from_height h r height energy_nonnegative
    energy_below_one height_relation
  have i_energy : 99 < 1000 * eI := by grind
  have margin := paid_reference_norm_margin base f S eB eI total reference_floor
    endpoint_upper short_pressure paid_reference i_energy original_total
  grind

end CoxeterProofs.CubicB16ShortCapKernel

#print axioms CoxeterProofs.CubicB16ShortCapKernel.cubic_coefficients_unique
#print axioms CoxeterProofs.CubicB16ShortCapKernel.i_strip_coefficients_unique
#print axioms CoxeterProofs.CubicB16ShortCapKernel.i_height_lower_from_unit
#print axioms CoxeterProofs.CubicB16ShortCapKernel.i_energy_lower_from_height
#print axioms CoxeterProofs.CubicB16ShortCapKernel.paid_reference_norm_margin
#print axioms CoxeterProofs.CubicB16ShortCapKernel.cubic_short_cap_contradiction
