import CoxeterProofs.RightAngledProjectionHeightKernel

/-!
Coordinate proof of the ideal 24-cell equality step. Source preparation
only; compilation and actual axiom outputs require a separate bounded run.

The two normals below are an explicit opposite spatial pair in the written
24-cell model. Their unit norms, height identities, future time bound and
center conclusion are derived. Metric identification with an actual facet,
published rigidity, and projection geometry remain application inputs.
-/

namespace CoxeterProofs.RegularIdeal24CellHeightKernel

open Lean.Grind
open CoxeterProofs.RightAngledCubicScalarKernel
open CoxeterProofs.RightAngledProjectionHeightKernel

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

def positiveNormal : Vector R := (1, fun i => if i < 2 then 1 else 0)
def negativeNormal : Vector R := (1, fun i => if i < 2 then -1 else 0)

/-- Unit normals and their two actual coordinate heights, including the sum. -/
theorem opposite_pair_identities (q : Vector R) :
    lorentz 4 (positiveNormal (R := R)) positiveNormal = 1 ∧
    lorentz 4 (negativeNormal (R := R)) negativeNormal = 1 ∧
    height 4 q positiveNormal = q.1 - q.2 0 - q.2 1 ∧
    height 4 q negativeNormal = q.1 + q.2 0 + q.2 1 ∧
    height 4 q positiveNormal + height 4 q negativeNormal = 2 * q.1 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  all_goals
    simp [height, lorentz, positiveNormal, negativeNormal, coordinateDot, sumN] <;> grind

/-- Future unit timelikeness derives the time bound, rather than assuming it. -/
theorem future_unit_time_lower_bound (q : Vector R)
    (future : 0 < q.1) (unit : lorentz 4 q q = -1) : 1 ≤ q.1 := by
  by_cases bound : 1 ≤ q.1
  · exact bound
  · have small : q.1 < 1 := by grind
    have square_small := OrderedRing.mul_lt_mul_of_pos_right small future
    have spatial_nonnegative := coordinateDot_self_nonneg 4 q.2
    unfold lorentz at unit
    grind

/-- The two primitive height caps force time one and every spatial coordinate zero. -/
theorem opposite_height_caps_force_center (q : Vector R)
    (future : 0 < q.1) (unit : lorentz 4 q q = -1)
    (positive_cap : height 4 q positiveNormal ≤ 1)
    (negative_cap : height 4 q negativeNormal ≤ 1) :
    q.1 = 1 ∧ ∀ i, i < 4 → q.2 i = 0 := by
  have lower := future_unit_time_lower_bound q future unit
  have pair := (opposite_pair_identities q).2.2.2.2
  have time_one : q.1 = 1 := by grind
  have nonnegative (i : Nat) : 0 ≤ (q.2 i)^2 := OrderedRing.sq_nonneg
  have nonnegative0 := nonnegative 0
  have nonnegative1 := nonnegative 1
  have nonnegative2 := nonnegative 2
  have nonnegative3 := nonnegative 3
  have square_sum : (q.2 0)^2 + (q.2 1)^2 +
      (q.2 2)^2 + (q.2 3)^2 = 0 := by
    unfold lorentz coordinateDot at unit
    simp only [sumN] at unit
    rw [time_one] at unit
    grind
  refine ⟨time_one, ?_⟩
  intro i hi
  have cases : i = 0 ∨ i = 1 ∨ i = 2 ∨ i = 3 := by omega
  have square_zero : (q.2 i)^2 = 0 := by
    rcases cases with h | h | h | h <;> subst i <;> grind
  by_cases zero : q.2 i = 0
  · exact zero
  · have positive := OrderedRing.sq_pos zero
    grind

/-- At that coordinate center, every normal with time coordinate one has height one. -/
theorem center_has_unit_heights (q m : Vector R)
    (time_one : q.1 = 1) (spatial_zero : ∀ i, i < 4 → q.2 i = 0)
    (normal_time : m.1 = 1) : height 4 q m = 1 := by
  have zero0 := spatial_zero 0 (by omega)
  have zero1 := spatial_zero 1 (by omega)
  have zero2 := spatial_zero 2 (by omega)
  have zero3 := spatial_zero 3 (by omega)
  unfold height lorentz coordinateDot
  simp only [sumN]
  rw [time_one, normal_time, zero0, zero1, zero2, zero3]
  grind

/-- Strict caps on just this opposite pair contradict future unit timelikeness. -/
theorem no_strict_opposite_height_caps (q : Vector R)
    (future : 0 < q.1) (unit : lorentz 4 q q = -1)
    (positive_cap : height 4 q positiveNormal < 1)
    (negative_cap : height 4 q negativeNormal < 1) : False := by
  have lower := future_unit_time_lower_bound q future unit
  have pair := (opposite_pair_identities q).2.2.2.2
  grind

end CoxeterProofs.RegularIdeal24CellHeightKernel

#print axioms CoxeterProofs.RegularIdeal24CellHeightKernel.opposite_pair_identities
#print axioms CoxeterProofs.RegularIdeal24CellHeightKernel.future_unit_time_lower_bound
#print axioms CoxeterProofs.RegularIdeal24CellHeightKernel.opposite_height_caps_force_center
#print axioms CoxeterProofs.RegularIdeal24CellHeightKernel.center_has_unit_heights
#print axioms CoxeterProofs.RegularIdeal24CellHeightKernel.no_strict_opposite_height_caps
