import Std

/-!
PREPARED SOURCE ONLY: no compilation or actual axiom audit is credited.

Ordered-ring scalar logic for the settled whole-A18 two-isolated-dotted
endpoint theorem. Three radical chords, a radical tangent, all quadratic
seals, both endpoint gains and the complete three-budget contradiction
are derived from primitive square equations and rational data.

External: actual endpoint/core geometry, the complete fresh N19 inventory,
all-wall feasibility, legal individual boundary prices, one original
reference/vector, residual orthogonality and reference norm floors.
No real square-root function, polytope or dimension instance is supplied.
-/

namespace CoxeterProofs.A18TwoDottedEndpointScalarKernel

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

private theorem chord_of_positive_square_gap (a l : R)
    (root_nonnegative : 0 ≤ a)
    (square_gap : 0 < 100000000 * a ^ 2 - l ^ 2) :
    l < 10000 * a := by
  by_cases target : l < 10000 * a
  · exact target
  · have first : 0 ≤ l - 10000 * a := by grind
    have second : 0 ≤ l + 10000 * a := by grind
    have product := OrderedRing.mul_nonneg first second
    have expansion :
        (l - 10000 * a) * (l + 10000 * a) =
          -(100000000 * a ^ 2 - l ^ 2) := by grind only
    rw [expansion] at product
    grind

/-- All three strict chords follow from the primitive root equation. -/
theorem first_radical_chords (u a : R)
    (interval_lower : 0 ≤ u) (interval_upper : u ≤ 1)
    (root_nonnegative : 0 ≤ a)
    (root_square : 19 * a ^ 2 = 36 - 19 * u) :
    (2 * u ≤ 1 → 13760 - 3902 * u < 10000 * a) ∧
    (1 ≤ 2 * u → 4 * u ≤ 3 → 14029 - 4440 * u < 10000 * a) ∧
    (3 ≤ 4 * u → 14419 - 4960 * u < 10000 * a) := by
  have left : 2 * u ≤ 1 → 13760 - 3902 * u < 10000 * a := by
    intro half_upper
    have product := OrderedRing.mul_nonneg interval_lower
      (show 0 ≤ 1 - 2 * u by grind)
    have identity :
        38 * (100000000 * a ^ 2 - (13760 - 3902 * u) ^ 2) =
          289286476 * u * (1 - 2 * u) +
            4365478 * (1 - 2 * u) + 805722 := by grind only
    have scaled_positive :
        0 < 38 * (100000000 * a ^ 2 - (13760 - 3902 * u) ^ 2) := by
      rw [identity]
      grind
    have gap : 0 < 100000000 * a ^ 2 - (13760 - 3902 * u) ^ 2 := by
      grind
    exact chord_of_positive_square_gap a (13760 - 3902 * u)
      root_nonnegative gap
  have middle :
      1 ≤ 2 * u → 4 * u ≤ 3 → 14029 - 4440 * u < 10000 * a := by
    intro half_lower three_quarters_upper
    have product := OrderedRing.mul_nonneg
      (show 0 ≤ 2 * u - 1 by grind) (show 0 ≤ 3 - 4 * u by grind)
    have identity :
        19 * (100000000 * a ^ 2 - (14029 - 4440 * u) ^ 2) =
          46819800 * (2 * u - 1) * (3 - 4 * u) +
            306280 * (3 - 4 * u) + 96581 := by grind only
    have scaled_positive :
        0 < 19 * (100000000 * a ^ 2 - (14029 - 4440 * u) ^ 2) := by
      rw [identity]
      grind
    have gap : 0 < 100000000 * a ^ 2 - (14029 - 4440 * u) ^ 2 := by
      grind
    exact chord_of_positive_square_gap a (14029 - 4440 * u)
      root_nonnegative gap
  have right : 3 ≤ 4 * u → 14419 - 4960 * u < 10000 * a := by
    intro three_quarters_lower
    have product := OrderedRing.mul_nonneg
      (show 0 ≤ 4 * u - 3 by grind) (show 0 ≤ 1 - u by grind)
    have identity :
        19 * (100000000 * a ^ 2 - (14419 - 4960 * u) ^ 2) =
          116857600 * (4 * u - 3) * (1 - u) +
            310080 * (1 - u) + 19061 := by grind only
    have scaled_positive :
        0 < 19 * (100000000 * a ^ 2 - (14419 - 4960 * u) ^ 2) := by
      rw [identity]
      grind
    have gap : 0 < 100000000 * a ^ 2 - (14419 - 4960 * u) ^ 2 := by
      grind
    exact chord_of_positive_square_gap a (14419 - 4960 * u)
      root_nonnegative gap
  exact ⟨left, middle, right⟩

/-- The full tangent is the square (2*b-1)^2, not an assumed estimate. -/
theorem second_radical_tangent (u b : R)
    (root_square : b ^ 2 = u * (1 - u)) :
    4 * b ≤ 1 + 4 * u - 4 * u ^ 2 := by
  have square_nonnegative : 0 ≤ (2 * b - 1) ^ 2 := OrderedRing.sq_nonneg
  have expansion : (2 * b - 1) ^ 2 = 1 + 4 * u - 4 * u ^ 2 - 4 * b := by
    grind only
  rw [expansion] at square_nonnegative
  grind

/-- Interval factors and one global completion prove all three seals. -/
theorem three_quadratic_seals (u : R) :
    (2 * u ≤ 1 → 0 < 7704 - 22804 * u + 17500 * u ^ 2) ∧
    (0 < 721175 - 2089500 * u + 1531250 * u ^ 2) ∧
    (3 ≤ 4 * u → 0 < 36088 - 99680 * u + 70000 * u ^ 2) := by
  have left : 2 * u ≤ 1 → 0 < 7704 - 22804 * u + 17500 * u ^ 2 := by
    intro upper
    have product := OrderedRing.mul_nonneg
      (show 0 ≤ 1 - 2 * u by grind) (show 0 ≤ 7027 - 8750 * u by grind)
    have identity :
        7704 - 22804 * u + 17500 * u ^ 2 =
          677 + (1 - 2 * u) * (7027 - 8750 * u) := by grind only
    rw [identity]
    grind
  have middle : 0 < 721175 - 2089500 * u + 1531250 * u ^ 2 := by
    have square_nonnegative : 0 ≤ (875 * u - 597) ^ 2 := OrderedRing.sq_nonneg
    have identity :
        721175 - 2089500 * u + 1531250 * u ^ 2 =
          2 * (875 * u - 597) ^ 2 + 8357 := by grind only
    rw [identity]
    grind
  have right : 3 ≤ 4 * u → 0 < 36088 - 99680 * u + 70000 * u ^ 2 := by
    intro lower
    have product := OrderedRing.mul_nonneg
      (show 0 ≤ 4 * u - 3 by grind) (show 0 ≤ 3500 * u - 2359 by grind)
    have identity :
        36088 - 99680 * u + 70000 * u ^ 2 =
          703 + 5 * (4 * u - 3) * (3500 * u - 2359) := by grind only
    rw [identity]
    grind
  exact ⟨left, middle, right⟩

/-- Entire closed-interval radical gap from primitive nonnegative roots. -/
theorem continuous_core_gap (u a b : R)
    (interval_lower : 0 ≤ u) (interval_upper : u ≤ 1)
    (first_nonnegative : 0 ≤ a) (second_nonnegative : 0 ≤ b)
    (first_square : 19 * a ^ 2 = 36 - 19 * u)
    (second_square : b ^ 2 = u * (1 - u)) :
    7408 - 2500 * u + 1250 * u ^ 2 < 10000 * (a - b) := by
  obtain ⟨left_chord, middle_chord, right_chord⟩ :=
    first_radical_chords u a interval_lower interval_upper first_nonnegative first_square
  have tangent := second_radical_tangent u b second_square
  obtain ⟨left_seal, middle_seal, right_seal⟩ := three_quadratic_seals (R := R) u
  by_cases left : 2 * u ≤ 1
  · have chord := left_chord left
    have positive := left_seal left
    grind
  · have half_lower : 1 ≤ 2 * u := by grind
    by_cases middle : 4 * u ≤ 3
    · have chord := middle_chord half_lower middle
      grind
    · have three_quarters_lower : 3 ≤ 4 * u := by grind
      have chord := right_chord three_quarters_lower
      have positive := right_seal three_quarters_lower
      grind

/-- Both genuine block gains, including the exact 1099/60000 difference. -/
theorem endpoint_gains (u a b g5 g6 : R)
    (interval_lower : 0 ≤ u) (interval_upper : u ≤ 1)
    (first_nonnegative : 0 ≤ a) (second_nonnegative : 0 ≤ b)
    (first_square : 19 * a ^ 2 = 36 - 19 * u)
    (second_square : b ^ 2 = u * (1 - u))
    (five_gain_identity :
      10000 * g5 = 5000 * (a - b) - 3104 + 1250 * u - 625 * u ^ 2)
    (six_gain_identity :
      60000 * g6 = 30000 * (a - b) - 17525 + 7500 * u - 3750 * u ^ 2) :
    600 < 10000 * g5 ∧ 4699 < 60000 * g6 ∧
      60000 * g6 = 60000 * g5 + 1099 := by
  have core := continuous_core_gap u a b interval_lower interval_upper
    first_nonnegative second_nonnegative first_square second_square
  constructor
  · grind
  · constructor
    · grind
    · grind

/-- ONE pressure, ONE squared norm and ONE nonnegative completed square. -/
theorem simultaneous_paid_completion
    (base kL kR hL hR sL sR pressure normW gainL gainR total : R)
    (one_pressure : 4 * base + kL + kR - hL - hR ≤ 4 * pressure)
    (one_norm : 16 * normW = 16 * base + sL + sR)
    (left_gain_identity : 16 * gainL = 8 * (kL - hL) - sL)
    (right_gain_identity : 16 * gainR = 8 * (kR - hR) - sR)
    (one_square : 0 ≤ total - 2 * pressure + normW) :
    base + gainL + gainR ≤ total := by
  grind

/-- k=0: one paid original reference with two five-block contributions. -/
theorem zero_weak_end_budget (base left right total : R)
    (base_lower : 9240352 ≤ 10000000 * base)
    (left_gain : 600 < 10000 * left) (right_gain : 600 < 10000 * right)
    (one_paid_budget : base + left + right ≤ total) :
    10440352 < 10000000 * total ∧ 1 < total := by
  have exact_margin : 10440352 < 10000000 * total := by grind
  constructor
  · exact exact_margin
  · grind

/-- k=1: either orientation has the same exact one-reference margin. -/
theorem one_weak_end_budget (base five six total : R)
    (base_lower : 8937 ≤ 10000 * base)
    (five_gain : 600 < 10000 * five) (six_gain : 4699 < 60000 * six)
    (one_paid_budget : base + five + six ≤ total) :
    61921 < 60000 * total ∧ 1 < total := by
  have exact_margin : 61921 < 60000 * total := by grind
  constructor
  · exact exact_margin
  · grind

/-- k=2: both six-block gains are paid on the same original vector. -/
theorem two_weak_end_budget (base left right total : R)
    (base_lower : 86414 ≤ 100000 * base)
    (left_gain : 4699 < 60000 * left) (right_gain : 4699 < 60000 * right)
    (one_paid_budget : base + left + right ≤ total) :
    153116 < 150000 * total ∧ 1 < total := by
  have exact_margin : 153116 < 150000 * total := by grind
  constructor
  · exact exact_margin
  · grind

/-- Complete scalar contradiction; no final gain or norm margin is a premise. -/
theorem paired_endpoint_budget_contradiction
    (u v a b c d g5L g6L g5R g6R base total : R)
    (left_lower : 0 ≤ u) (left_upper : u ≤ 1)
    (left_first_nonnegative : 0 ≤ a) (left_second_nonnegative : 0 ≤ b)
    (left_first_square : 19 * a ^ 2 = 36 - 19 * u)
    (left_second_square : b ^ 2 = u * (1 - u))
    (left_five_identity :
      10000 * g5L = 5000 * (a - b) - 3104 + 1250 * u - 625 * u ^ 2)
    (left_six_identity :
      60000 * g6L = 30000 * (a - b) - 17525 + 7500 * u - 3750 * u ^ 2)
    (right_lower : 0 ≤ v) (right_upper : v ≤ 1)
    (right_first_nonnegative : 0 ≤ c) (right_second_nonnegative : 0 ≤ d)
    (right_first_square : 19 * c ^ 2 = 36 - 19 * v)
    (right_second_square : d ^ 2 = v * (1 - v))
    (right_five_identity :
      10000 * g5R = 5000 * (c - d) - 3104 + 1250 * v - 625 * v ^ 2)
    (right_six_identity :
      60000 * g6R = 30000 * (c - d) - 17525 + 7500 * v - 3750 * v ^ 2)
    (original_norm_bound : total < 1)
    (legal_pattern_budget :
      (9240352 ≤ 10000000 * base ∧ base + g5L + g5R ≤ total) ∨
      (8937 ≤ 10000 * base ∧
        (base + g5L + g6R ≤ total ∨ base + g6L + g5R ≤ total)) ∨
      (86414 ≤ 100000 * base ∧ base + g6L + g6R ≤ total)) :
    False := by
  obtain ⟨left_five, left_six, left_difference⟩ := endpoint_gains u a b g5L g6L
    left_lower left_upper left_first_nonnegative left_second_nonnegative
    left_first_square left_second_square left_five_identity left_six_identity
  obtain ⟨right_five, right_six, right_difference⟩ := endpoint_gains v c d g5R g6R
    right_lower right_upper right_first_nonnegative right_second_nonnegative
    right_first_square right_second_square right_five_identity right_six_identity
  cases legal_pattern_budget with
  | inl zero =>
    have margin := zero_weak_end_budget base g5L g5R total
      zero.1 left_five right_five zero.2
    grind
  | inr rest =>
    cases rest with
    | inl one =>
      cases one.2 with
      | inl left_strong =>
        have margin := one_weak_end_budget base g5L g6R total
          one.1 left_five right_six left_strong
        grind
      | inr right_strong =>
        have paid : base + g5R + g6L ≤ total := by grind
        have margin := one_weak_end_budget base g5R g6L total
          one.1 right_five left_six paid
        grind
    | inr two =>
      have margin := two_weak_end_budget base g6L g6R total
        two.1 left_six right_six two.2
      grind

end CoxeterProofs.A18TwoDottedEndpointScalarKernel

#print axioms CoxeterProofs.A18TwoDottedEndpointScalarKernel.first_radical_chords
#print axioms CoxeterProofs.A18TwoDottedEndpointScalarKernel.second_radical_tangent
#print axioms CoxeterProofs.A18TwoDottedEndpointScalarKernel.three_quadratic_seals
#print axioms CoxeterProofs.A18TwoDottedEndpointScalarKernel.continuous_core_gap
#print axioms CoxeterProofs.A18TwoDottedEndpointScalarKernel.endpoint_gains
#print axioms CoxeterProofs.A18TwoDottedEndpointScalarKernel.simultaneous_paid_completion
#print axioms CoxeterProofs.A18TwoDottedEndpointScalarKernel.zero_weak_end_budget
#print axioms CoxeterProofs.A18TwoDottedEndpointScalarKernel.one_weak_end_budget
#print axioms CoxeterProofs.A18TwoDottedEndpointScalarKernel.two_weak_end_budget
#print axioms CoxeterProofs.A18TwoDottedEndpointScalarKernel.paired_endpoint_budget_contradiction
