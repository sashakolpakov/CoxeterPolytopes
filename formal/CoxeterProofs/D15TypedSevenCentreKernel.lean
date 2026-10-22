import Std

/-!
PREPARED source: no compilation or actual axiom audit is credited here.

Ordered-ring logic for the tight D15+actual-typed-seven centre argument.
The radical chord and tangent are derived from their square equations;
the continuous margin is proved for the entire interval, not sampled.

External application inputs: actual field/white classifications, the
original parameter and radical interpretation, legal simultaneous P5/P11
payments and orthogonality, defect energy and positive reference floors.
No geometric point, quadratic form or polytope instance is constructed.
-/

namespace CoxeterProofs.D15TypedSevenCentreKernel

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- Rational bounds for the nonnegative root of six, proved by squares. -/
theorem sqrt_six_seals (s : R) (nonnegative : 0 ≤ s) (square : s ^ 2 = 6) :
    17 < 7 * s ∧ 2000 * s < 4899 := by
  have lower : 17 < 7 * s := by
    by_cases target : 17 < 7 * s
    · exact target
    · have first : 0 ≤ 17 - 7 * s := by grind
      have second : 0 ≤ 17 + 7 * s := by grind
      have product := OrderedRing.mul_nonneg first second
      have identity : (17 - 7 * s) * (17 + 7 * s) = -5 := by
        grind only
      rw [identity] at product
      grind
  have upper : 2000 * s < 4899 := by
    by_cases target : 2000 * s < 4899
    · exact target
    · have first : 0 ≤ 2000 * s - 4899 := by grind
      have second : 0 ≤ 2000 * s + 4899 := by grind
      have product := OrderedRing.mul_nonneg first second
      have identity : (2000 * s - 4899) * (2000 * s + 4899) = -201 := by
        grind only
      rw [identity] at product
      grind
  exact ⟨lower, upper⟩

/-- A stronger interval seal proves the actual radical chord internally. -/
theorem first_radical_chord (x a : R)
    (interval_lower : 1 ≤ x) (interval_upper : 4 * x ≤ 5)
    (root_nonnegative : 0 ≤ a) (root_square : a ^ 2 = 15 - 2 * x) :
    777 - 56 * x < 200 * a := by
  have product := OrderedRing.mul_nonneg
    (show 0 ≤ x - 1 by grind) (show 0 ≤ 5 - 4 * x by grind)
  have identity :
      40000 * a ^ 2 - (777 - 56 * x) ^ 2 =
        784 * (x - 1) * (5 - 4 * x) + 191 - 32 * x := by
    grind only
  have positive : 0 < 40000 * a ^ 2 - (777 - 56 * x) ^ 2 := by
    rw [identity]
    grind
  by_cases target : 777 - 56 * x < 200 * a
  · exact target
  · have first : 0 ≤ (777 - 56 * x) - 200 * a := by grind
    have second : 0 ≤ (777 - 56 * x) + 200 * a := by grind
    have contradiction := OrderedRing.mul_nonneg first second
    have expansion :
        ((777 - 56 * x) - 200 * a) * ((777 - 56 * x) + 200 * a) =
          - (40000 * a ^ 2 - (777 - 56 * x) ^ 2) := by grind only
    rw [expansion] at contradiction
    grind

/-- The square (s*b-1)^2 supplies the full radical tangent inequality. -/
theorem second_radical_tangent (s z b : R)
    (s_nonnegative : 0 ≤ s) (s_square : s ^ 2 = 6)
    (root_square : 6 * b ^ 2 = 1 - 36 * z ^ 2) :
    6 * b + 18 * s * z ^ 2 ≤ s := by
  have nonnegative_square : 0 ≤ (s * b - 1) ^ 2 := OrderedRing.sq_nonneg
  have expansion : (s * b - 1) ^ 2 = 2 - 36 * z ^ 2 - 2 * s * b := by
    grind only
  rw [expansion] at nonnegative_square
  have tangent : s * b + 18 * z ^ 2 ≤ 1 := by grind
  have scaled := OrderedRing.mul_le_mul_of_nonneg_left tangent s_nonnegative
  have reduction : s * (s * b + 18 * z ^ 2) = 6 * b + 18 * s * z ^ 2 := by
    grind only
  rw [reduction] at scaled
  grind

/-- Uniform 14G-R > 31.59 from primitive square equations and interval. -/
theorem centre_radical_gap (x z a b s norm : R)
    (interval_lower : 1 ≤ x) (interval_upper : 4 * x ≤ 5)
    (substitution : 6 * z = 6 * x - 7)
    (first_nonnegative : 0 ≤ a) (first_square : a ^ 2 = 15 - 2 * x)
    (s_nonnegative : 0 ≤ s) (s_square : s ^ 2 = 6)
    (second_square : 6 * b ^ 2 = 1 - 36 * z ^ 2)
    (norm_identity : 2 * norm = 25 - 4 * z + 12 * z ^ 2) :
    3159 < 100 * (14 * (a - b) - norm) := by
  obtain ⟨s_lower, s_upper⟩ := sqrt_six_seals s s_nonnegative s_square
  have chord := first_radical_chord x a interval_lower interval_upper
    first_nonnegative first_square
  have tangent := second_radical_tangent s z b s_nonnegative s_square second_square
  have polynomial_lower :
      223900 - 14000 * s - 11520 * z + (252000 * s - 36000) * z ^ 2 <
        6000 * (14 * (a - b) - norm) := by grind
  have z_nonnegative : 0 ≤ z ^ 2 := OrderedRing.sq_nonneg
  have coefficient_nonnegative : 0 ≤ 252000 * s - 612000 := by grind
  have coefficient_product := OrderedRing.mul_nonneg coefficient_nonnegative z_nonnegative
  have simpler :
      189607 - 11520 * z + 576000 * z ^ 2 <
        6000 * (14 * (a - b) - norm) := by grind
  have square_nonnegative : 0 ≤ (100 * z - 1) ^ 2 := OrderedRing.sq_nonneg
  have completion :
      5 * (189607 - 11520 * z + 576000 * z ^ 2) =
        288 * (100 * z - 1) ^ 2 + 947747 := by grind only
  have sharp : 947747 < 30000 * (14 * (a - b) - norm) := by grind
  grind

/-- One legal joint head/centre reference and one actual defect charge. -/
theorem funded_centre_margin (base gain gap height eD eZ total : R)
    (base_lower : 84816 ≤ 100000 * base)
    (head_gain : 2080125 ≤ 100000000 * gain)
    (scalar_gap : 3159 < 100 * gap)
    (height_upper : 1000 * height ≤ 2886)
    (paid_reference :
      49 * base + 49 * gain + gap - 14 * height + 11 ≤ 49 * eD)
    (defect_lower : 9 < 100 * eZ) (original_total : total = eD + eZ) :
    100357 < 100000 * total ∧ 1 < total := by
  have exact_margin : 100357 < 100000 * total := by grind
  constructor
  · exact exact_margin
  · grind

/-- After the complete centre classification, the final positive reference. -/
theorem strong_centre_margin (base gain eD eZ total : R)
    (base_lower : 88989 ≤ 100000 * base)
    (head_gain : 2080125 ≤ 100000000 * gain)
    (paid_reference : base + gain ≤ eD)
    (defect_lower : 9 < 100 * eZ) (original_total : total = eD + eZ) :
    100069125 < 100000000 * total ∧ 1 < total := by
  have exact_margin : 100069125 < 100000000 * total := by grind
  constructor
  · exact exact_margin
  · grind

end CoxeterProofs.D15TypedSevenCentreKernel

#print axioms CoxeterProofs.D15TypedSevenCentreKernel.sqrt_six_seals
#print axioms CoxeterProofs.D15TypedSevenCentreKernel.first_radical_chord
#print axioms CoxeterProofs.D15TypedSevenCentreKernel.second_radical_tangent
#print axioms CoxeterProofs.D15TypedSevenCentreKernel.centre_radical_gap
#print axioms CoxeterProofs.D15TypedSevenCentreKernel.funded_centre_margin
#print axioms CoxeterProofs.D15TypedSevenCentreKernel.strong_centre_margin
