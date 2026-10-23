import Std

/-!
PREPARED SOURCE ONLY: no compilation or actual axiom audit is credited.

The full continuous scalar budget for the arithmetic A13+actual-I+ordinary-T2
consumer. Literal gain identities, two polynomial squares and interval factors
derive all ranges. Geometry, inventories, actual pressures, individual prices,
weak-foot absorption and the interpretation of ONE completed square stay external.
No Real/polytope/number-field instance or geometric existence is supplied.
-/

namespace CoxeterProofs.A13IEndpointGraftScalarKernel

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- A single exact square proves the nonhelper middle gain on ALL parameters. -/
theorem nonhelper_middle_quadratic (s : R) :
    217000 < 2500000 * s ^ 2 + (1163 - 3615 * s) ^ 2 := by
  have square : 0 ≤ (15568225 * s - 4204245) ^ 2 := OrderedRing.sq_nonneg
  have identity :
      15568225 * (2500000 * s ^ 2 + (1163 - 3615 * s) ^ 2 - 217000) =
        (15568225 * s - 4204245) ^ 2 + 3117675000 := by
    grind only
  have positive :
      0 < 15568225 * (2500000 * s ^ 2 + (1163 - 3615 * s) ^ 2 - 217000) := by
    rw [identity]
    grind
  grind

/-- The helper middle gain is derived from its own exact square and gap. -/
theorem helper_middle_quadratic (s : R) :
    385000 < 5000000 * s ^ 2 + (1181 - 3615 * s) ^ 2 := by
  have square : 0 ≤ (18068225 * s - 4269315) ^ 2 := OrderedRing.sq_nonneg
  have identity :
      18068225 * (5000000 * s ^ 2 + (1181 - 3615 * s) ^ 2 - 385000) =
        (18068225 * s - 4269315) ^ 2 + 17538375000 := by
    grind only
  have positive :
      0 < 18068225 * (5000000 * s ^ 2 + (1181 - 3615 * s) ^ 2 - 385000) := by
    rw [identity]
    grind
  grind

/-- The fixed two-end graft covers the ENTIRE low interval, including equality. -/
theorem nonhelper_fixed_gain (s gain : R)
    (upper : 5 * s ≤ 1) (gain_identity : 3125 * gain = 963 - 3615 * s) :
    73 ≤ 625 * (s ^ 2 + gain) := by
  have product := OrderedRing.mul_nonneg
    (show 0 ≤ 1 - 5 * s by grind) (show 0 ≤ 598 - 625 * s by grind)
  have identity :
      3125 * (s ^ 2 + gain) - 365 = (1 - 5 * s) * (598 - 625 * s) := by
    grind only
  rw [← identity] at product
  grind

/-- The opposite-block helper graft has a separate exact low-interval factor. -/
theorem helper_fixed_gain (s gain : R)
    (upper : 5 * s ≤ 1) (gain_identity : 6250 * gain = 981 - 3615 * s) :
    254 ≤ 3125 * (s ^ 2 + gain) := by
  have product := OrderedRing.mul_nonneg
    (show 0 ≤ 1 - 5 * s by grind) (show 0 ≤ 473 - 1250 * s by grind)
  have identity :
      6250 * (s ^ 2 + gain) - 508 = (1 - 5 * s) * (473 - 1250 * s) := by
    grind only
  rw [← identity] at product
  grind

/-- High actual-I parameter values pay at least .09 directly. -/
theorem high_parameter_square (s : R) (nonnegative : 0 ≤ s)
    (lower : 3 ≤ 10 * s) : 9 ≤ 100 * s ^ 2 := by
  have product := OrderedRing.mul_nonneg
    (show 0 ≤ 10 * s - 3 by grind) (show 0 ≤ 10 * s + 3 by grind)
  grind

/-- The energy is deduced from ONE paid pressure/norm and ONE square. -/
theorem paid_energy_completion (base gain pressure normW energy : R)
    (paid_pressure_norm : base + gain ≤ 2 * pressure - normW)
    (one_square : 0 ≤ energy - 2 * pressure + normW) :
    base + gain ≤ energy := by
  grind

/-- Ordinary data, not an instance: actual joint-vector interfaces stay visible. -/
structure PaidGraftData (R : Type u) [CommRing R] [LE R] [LT R]
    [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]
    (s energy : R) where
  base : R
  lowGain : R
  middleGain : R
  lowPressure : R
  lowNorm : R
  middlePressure : R
  middleNorm : R
  component_floor : base ≤ energy
  low_paid : 5 * s ≤ 1 → base + lowGain ≤ 2 * lowPressure - lowNorm
  low_square : 5 * s ≤ 1 → 0 ≤ energy - 2 * lowPressure + lowNorm
  middle_paid : 1 ≤ 5 * s → 10 * s ≤ 3 →
    base + middleGain ≤ 2 * middlePressure - middleNorm
  middle_square : 1 ≤ 5 * s → 10 * s ≤ 3 →
    0 ≤ energy - 2 * middlePressure + middleNorm

/-- Trichotomy derives a uniform nonhelper gain; no final gain bound is assumed. -/
theorem nonhelper_all_ranges (s energy isolated : R)
    (data : PaidGraftData R s energy) (parameter_positive : 0 < s)
    (isolated_floor : s ^ 2 ≤ isolated)
    (low_identity : 3125 * data.lowGain = 963 - 3615 * s)
    (middle_identity : 2500000 * data.middleGain = (1163 - 3615 * s) ^ 2) :
    217 < 2500 * (energy + isolated - data.base) := by
  by_cases low : 5 * s ≤ 1
  · have gain := nonhelper_fixed_gain s data.lowGain low low_identity
    have paid := paid_energy_completion data.base data.lowGain
      data.lowPressure data.lowNorm energy (data.low_paid low) (data.low_square low)
    grind
  · have middle_lower : 1 ≤ 5 * s := by grind
    by_cases middle : 10 * s ≤ 3
    · have scalar := nonhelper_middle_quadratic (R := R) s
      rw [← middle_identity] at scalar
      have gain : 217 < 2500 * (s ^ 2 + data.middleGain) := by grind
      have paid := paid_energy_completion data.base data.middleGain
        data.middlePressure data.middleNorm energy
        (data.middle_paid middle_lower middle) (data.middle_square middle_lower middle)
      grind
    · have high := high_parameter_square s (by grind) (show 3 ≤ 10 * s by grind)
      have floor := data.component_floor
      grind

/-- The helper branch uses its one-end gain, with the same full trichotomy. -/
theorem helper_all_ranges (s energy isolated : R)
    (data : PaidGraftData R s energy) (parameter_positive : 0 < s)
    (isolated_floor : s ^ 2 ≤ isolated)
    (low_identity : 6250 * data.lowGain = 981 - 3615 * s)
    (middle_identity : 5000000 * data.middleGain = (1181 - 3615 * s) ^ 2) :
    77 < 1000 * (energy + isolated - data.base) := by
  by_cases low : 5 * s ≤ 1
  · have gain := helper_fixed_gain s data.lowGain low low_identity
    have paid := paid_energy_completion data.base data.lowGain
      data.lowPressure data.lowNorm energy (data.low_paid low) (data.low_square low)
    grind
  · have middle_lower : 1 ≤ 5 * s := by grind
    by_cases middle : 10 * s ≤ 3
    · have scalar := helper_middle_quadratic (R := R) s
      rw [← middle_identity] at scalar
      have gain : 77 < 1000 * (s ^ 2 + data.middleGain) := by grind
      have paid := paid_energy_completion data.base data.middleGain
        data.middlePressure data.middleNorm energy
        (data.middle_paid middle_lower middle) (data.middle_square middle_lower middle)
      grind
    · have high := high_parameter_square s (by grind) (show 3 ≤ 10 * s by grind)
      have floor := data.component_floor
      grind

/-- The literal nonhelper base plus the DERIVED continuous gain contradicts one norm. -/
theorem nonhelper_total_contradiction (base energy isolated dihedral total : R)
    (base_identity : 100000 * base = 75719)
    (derived_gain : 217 < 2500 * (energy + isolated - base))
    (ordinary_dihedral : 3 < 19 * dihedral)
    (one_orthogonal_sum : energy + isolated + dihedral ≤ total)
    (original_norm : total < 1) : False := by
  grind

/-- The helper base and its separately derived .077 gain have their own strict margin. -/
theorem helper_total_contradiction (base energy isolated dihedral total : R)
    (base_identity : 200000 * base = 153801)
    (derived_gain : 77 < 1000 * (energy + isolated - base))
    (ordinary_dihedral : 3 < 19 * dihedral)
    (one_orthogonal_sum : energy + isolated + dihedral ≤ total)
    (original_norm : total < 1) : False := by
  grind

/-- Complete scalar chain from primitive gain identities in EITHER actual case. -/
theorem complete_group_nine_scalar_contradiction (s energy isolated dihedral total : R)
    (data : PaidGraftData R s energy) (parameter_positive : 0 < s)
    (isolated_floor : s ^ 2 ≤ isolated)
    (ordinary_dihedral : 3 < 19 * dihedral)
    (one_orthogonal_sum : energy + isolated + dihedral ≤ total)
    (original_norm : total < 1)
    (actual_case :
      (100000 * data.base = 75719 ∧
        3125 * data.lowGain = 963 - 3615 * s ∧
        2500000 * data.middleGain = (1163 - 3615 * s) ^ 2) ∨
      (200000 * data.base = 153801 ∧
        6250 * data.lowGain = 981 - 3615 * s ∧
        5000000 * data.middleGain = (1181 - 3615 * s) ^ 2)) : False := by
  cases actual_case with
  | inl nonhelper =>
    have gain := nonhelper_all_ranges s energy isolated data parameter_positive
      isolated_floor nonhelper.2.1 nonhelper.2.2
    exact nonhelper_total_contradiction data.base energy isolated dihedral total
      nonhelper.1 gain ordinary_dihedral one_orthogonal_sum original_norm
  | inr helper =>
    have gain := helper_all_ranges s energy isolated data parameter_positive
      isolated_floor helper.2.1 helper.2.2
    exact helper_total_contradiction data.base energy isolated dihedral total
      helper.1 gain ordinary_dihedral one_orthogonal_sum original_norm

end CoxeterProofs.A13IEndpointGraftScalarKernel

#print axioms CoxeterProofs.A13IEndpointGraftScalarKernel.nonhelper_middle_quadratic
#print axioms CoxeterProofs.A13IEndpointGraftScalarKernel.helper_middle_quadratic
#print axioms CoxeterProofs.A13IEndpointGraftScalarKernel.nonhelper_fixed_gain
#print axioms CoxeterProofs.A13IEndpointGraftScalarKernel.helper_fixed_gain
#print axioms CoxeterProofs.A13IEndpointGraftScalarKernel.high_parameter_square
#print axioms CoxeterProofs.A13IEndpointGraftScalarKernel.paid_energy_completion
#print axioms CoxeterProofs.A13IEndpointGraftScalarKernel.nonhelper_all_ranges
#print axioms CoxeterProofs.A13IEndpointGraftScalarKernel.helper_all_ranges
#print axioms CoxeterProofs.A13IEndpointGraftScalarKernel.nonhelper_total_contradiction
#print axioms CoxeterProofs.A13IEndpointGraftScalarKernel.helper_total_contradiction
#print axioms CoxeterProofs.A13IEndpointGraftScalarKernel.complete_group_nine_scalar_contradiction
