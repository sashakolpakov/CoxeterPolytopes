import Std

/-!
Prepared scalar logic for the classical dimension-seventeen census closure.
No compilation or actual axiom audit is credited by this source header.

The full actual geometry, arithmetic field/cycle reductions, exhaustive
two-case census, radical interpretation, legal P5/P14 payments and the
one-vector paid-reference inequality remain explicit application inputs.
This module proves the strict scalar implications from those inputs.
It constructs no hyperbolic polytope or instance of the written theorem.
-/

namespace CoxeterProofs.ClassicalSeventeenCensusKernel

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- The exact .043056 head margin, from the three primitive scalar seals. -/
theorem survivor_head_gain (height baseCoordinate norm gain : R)
    (height_lower : 482 < 1000 * height)
    (base_upper : 1000 * baseCoordinate ≤ 361)
    (norm_upper : 100 * norm ≤ 54)
    (gain_identity :
      3125 * gain = 2250 * (height - baseCoordinate) - 405 * norm + 81) :
    43056 < 1000000 * gain ∧ 43 < 1000 * gain := by
  have exact_margin : 43056 < 1000000 * gain := by grind
  constructor
  · exact exact_margin
  · grind

/-- A strict positive-square comparison gives the fork pressure above3.7. -/
theorem fork_pressure_seal (pressure : R)
    (nonnegative : 0 ≤ pressure)
    (square_lower : 357 < 26 * pressure ^ 2) :
    37 < 10 * pressure := by
  by_cases target : 37 < 10 * pressure
  · exact target
  · have first : 0 ≤ 37 - 10 * pressure := by grind
    have second : 0 ≤ 37 + 10 * pressure := by grind
    have product := OrderedRing.mul_nonneg first second
    have expansion :
        (37 - 10 * pressure) * (37 + 10 * pressure) =
          1369 - 100 * pressure ^ 2 := by grind only
    rw [expansion] at product
    grind

/-- The full exact conservative1.00184 budget, with each payment already legal. -/
theorem coupled_energy_margin (base gain pressure oldHeight energy : R)
    (base_lower : 9084 ≤ 10000 * base)
    (gain_lower : 43 < 1000 * gain)
    (pressure_lower : 37 < 10 * pressure)
    (old_height_upper : 1000 * oldHeight ≤ 3311)
    (one_paid_reference :
      10000 * base + 10000 * gain + 2600 * (pressure - oldHeight) - 507 ≤
        10000 * energy) :
    100184 < 100000 * energy ∧ 1 < energy := by
  have exact_margin : 100184 < 100000 * energy := by grind
  constructor
  · exact exact_margin
  · grind

/-- Either census survivor is impossible once its actual budget is supplied. -/
theorem no_bounded_survivor (base gain pressure oldHeight energy : R)
    (base_lower : 9084 ≤ 10000 * base)
    (gain_lower : 43 < 1000 * gain)
    (pressure_lower : 37 < 10 * pressure)
    (old_height_upper : 1000 * oldHeight ≤ 3311)
    (one_paid_reference :
      10000 * base + 10000 * gain + 2600 * (pressure - oldHeight) - 507 ≤
        10000 * energy)
    (actual_norm_budget : energy < 1) : False := by
  have margin := coupled_energy_margin base gain pressure oldHeight energy
    base_lower gain_lower pressure_lower old_height_upper one_paid_reference
  grind

end CoxeterProofs.ClassicalSeventeenCensusKernel

#print axioms CoxeterProofs.ClassicalSeventeenCensusKernel.survivor_head_gain
#print axioms CoxeterProofs.ClassicalSeventeenCensusKernel.fork_pressure_seal
#print axioms CoxeterProofs.ClassicalSeventeenCensusKernel.coupled_energy_margin
#print axioms CoxeterProofs.ClassicalSeventeenCensusKernel.no_bounded_survivor
