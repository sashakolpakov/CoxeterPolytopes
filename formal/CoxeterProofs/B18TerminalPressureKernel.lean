import Std

/-!
PREPARED SOURCE ONLY: no compilation or actual axiom audit is credited.

Internal ordered-ring algebra for the B18 near3/own6 terminal theorem.
Primitive strict pressure bounds derive the entire unbounded core seal.
Primitive radical squares derive the subcritical root comparison and
the simultaneous two-sum contradiction, including zero mutual product.

External: actual normals/inventories, the original farthest vertex,
all-wall future feasibility, native/affine dispatch, terminal-sum bounds,
scalar meanings and their Real instance. No geometric instance is supplied.
-/

namespace CoxeterProofs.B18TerminalPressureKernel

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- The complete unbounded quadratic is positive, with no interval premise. -/
theorem core_quadratic_positive (y : R) :
    0 < 527 * y ^ 2 - 1360 * y + 1163 := by
  have square : 0 ≤ (527 * y - 680) ^ 2 := OrderedRing.sq_nonneg
  have identity :
      527 * (527 * y ^ 2 - 1360 * y + 1163) =
        (527 * y - 680) ^ 2 + 150501 := by grind only
  have scaled_positive : 0 < 527 * (527 * y ^ 2 - 1360 * y + 1163) := by
    rw [identity]
    grind
  grind

/-- Primitive J/delta and rho/delta bounds DERIVE the mixed quadratic term. -/
theorem mixed_pressure_product (rho delta J : R)
    (rho_positive : 0 < rho)
    (pressure_lower : 31 * delta < 100 * J)
    (correlation_upper : 5 * rho < 4 * delta) :
    155 * rho ^ 2 < 400 * rho * J := by
  have first := OrderedRing.mul_pos rho_positive
    (show 0 < 100 * J - 31 * delta by grind)
  have second := OrderedRing.mul_pos rho_positive
    (show 0 < 4 * delta - 5 * rho by grind)
  grind

/-- The actual gap is positive from primitive data, not an assumed final gain. -/
theorem supercritical_gap_positive (s rho delta H J y D : R)
    (s_positive : 0 < s) (rho_positive : 0 < rho)
    (head_lower : 13 < 20 * H)
    (pressure_lower : 31 * delta < 100 * J)
    (correlation_upper : 5 * rho < 4 * delta)
    (radical_square : 17 * s * y ^ 2 = 17 * rho ^ 2 - 9 * s)
    (actual_gap_identity : s * D = s * H + rho * J - s * y) :
    0 < D := by
  have mixed := mixed_pressure_product rho delta J rho_positive
    pressure_lower correlation_upper
  have head_product := OrderedRing.mul_pos s_positive
    (show 0 < 20 * H - 13 by grind)
  have identity :
      5 * (1360 * s * D - s * (527 * y ^ 2 - 1360 * y + 1163)) =
        340 * (20 * s * H - 13 * s) +
          17 * (400 * rho * J - 155 * rho ^ 2) := by grind only
  have scaled_gap :
      0 < 1360 * s * D - s * (527 * y ^ 2 - 1360 * y + 1163) := by
    have positive_right :
        0 < 340 * (20 * s * H - 13 * s) +
          17 * (400 * rho * J - 155 * rho ^ 2) := by grind
    rw [← identity] at positive_right
    grind
  have quadratic_product := OrderedRing.mul_pos s_positive
    (core_quadratic_positive (R := R) y)
  by_cases positive : 0 < D
  · exact positive
  · have product := OrderedRing.mul_nonneg
      (show 0 ≤ s by grind) (show 0 ≤ -D by grind)
    grind

/-- A geometrically supplied nonpositive actual gap contradicts the derived seal. -/
theorem supercritical_pressure_contradiction (s rho delta H J y D : R)
    (s_positive : 0 < s) (rho_positive : 0 < rho)
    (head_lower : 13 < 20 * H)
    (pressure_lower : 31 * delta < 100 * J)
    (correlation_upper : 5 * rho < 4 * delta)
    (radical_square : 17 * s * y ^ 2 = 17 * rho ^ 2 - 9 * s)
    (actual_gap_identity : s * D = s * H + rho * J - s * y)
    (actual_pressure : D ≤ 0) : False := by
  have positive := supercritical_gap_positive s rho delta H J y D
    s_positive rho_positive head_lower pressure_lower correlation_upper
    radical_square actual_gap_identity
  grind

/-- The base radical cap .357 is derived from its square and z<.1495. -/
theorem subcritical_base_root_bound (z a : R)
    (z_upper : 2000 * z < 299) (a_nonnegative : 0 ≤ a)
    (a_square : a ^ 2 = z - z ^ 2) :
    1000 * a < 357 := by
  have interval_product := OrderedRing.mul_pos
    (show 0 < 299 - 2000 * z by grind)
    (show 0 < 1701 - 2000 * z by grind)
  have identity :
      (299 - 2000 * z) * (1701 - 2000 * z) =
        508599 - 4000000 * a ^ 2 := by grind only
  rw [identity] at interval_product
  by_cases target : 1000 * a < 357
  · exact target
  · have product := OrderedRing.mul_nonneg
      (show 0 ≤ 1000 * a - 357 by grind)
      (show 0 ≤ 1000 * a + 357 by grind)
    grind

/-- Proper subcritical square data gives genuine convex weights 0<k<1. -/
theorem subcritical_parameter_bounds (z k : R)
    (z_lower : 1494 < 10000 * z) (z_upper : 2000 * z < 299)
    (k_lower : z ≤ k) (proper_square : k ^ 2 < z) :
    0 < z ∧ 0 < k ∧ k < 1 := by
  have z_positive : 0 < z := by grind
  have k_positive : 0 < k := by grind
  have k_upper : k < 1 := by
    by_cases target : k < 1
    · exact target
    · have product := OrderedRing.mul_nonneg
        (show 0 ≤ k - 1 by grind) (show 0 ≤ k + 1 by grind)
      grind
  exact ⟨z_positive, k_positive, k_upper⟩

/-- Both radical comparisons follow from squares; monotonicity is not a premise. -/
theorem subcritical_root_difference (z k a b : R)
    (z_lower : 1494 < 10000 * z) (z_upper : 2000 * z < 299)
    (k_lower : z ≤ k) (a_nonnegative : 0 ≤ a) (b_nonnegative : 0 ≤ b)
    (a_square : a ^ 2 = z - z ^ 2) (b_square : b ^ 2 = z - k ^ 2) :
    b ≤ a ∧ k - z ≤ 4 * (a - b) := by
  have base_cap := subcritical_base_root_bound z a z_upper a_nonnegative a_square
  have coordinate_product := OrderedRing.mul_nonneg
    (show 0 ≤ k - z by grind) (show 0 ≤ k + z by grind)
  have roots_order : b ≤ a := by
    by_cases target : b ≤ a
    · exact target
    · have product := OrderedRing.mul_pos
        (show 0 < b - a by grind) (show 0 < b + a by grind)
      grind
  have first := OrderedRing.mul_nonneg
    (show 0 ≤ k - z by grind) (show 0 ≤ 4 * (k + z) - 1 by grind)
  have second := OrderedRing.mul_nonneg
    (show 0 ≤ a - b by grind) (show 0 ≤ 1 - a - b by grind)
  have identity :
      (k - z) * (4 * (k + z) - 1) + 4 * ((a - b) * (1 - a - b)) =
        4 * (a - b) - (k - z) := by grind only
  have nonnegative : 0 ≤ 4 * (a - b) - (k - z) := by
    rw [← identity]
    grind
  exact ⟨roots_order, by grind⟩

/-- The WHOLE subcritical domain has strict pressure gap greater than .0135. -/
theorem subcritical_pressure_gap (z k a b t : R)
    (z_lower : 1494 < 10000 * z) (z_upper : 2000 * z < 299)
    (k_lower : z ≤ k) (a_nonnegative : 0 ≤ a) (b_nonnegative : 0 ≤ b)
    (a_square : a ^ 2 = z - z ^ 2) (b_square : b ^ 2 = z - k ^ 2)
    (t_lower : 4141 < 1000 * t) :
    7547 + 140 * k < 2000 * (t - b) := by
  have base_cap := subcritical_base_root_bound z a z_upper a_nonnegative a_square
  obtain ⟨_roots_order, difference⟩ := subcritical_root_difference z k a b
    z_lower z_upper k_lower a_nonnegative b_nonnegative a_square b_square
  have weaker_difference : 7 * (k - z) ≤ 100 * (a - b) := by grind
  grind

/-- The two ACTUAL sum bounds combine with derived nonnegative convex weights. -/
theorem subcritical_sum_pressure_upper (k S17 S18 : R)
    (k_positive : 0 < k) (k_upper : k < 1)
    (sum17_upper : 25 * S17 < 94) (sum18_upper : 100 * S18 < 383) :
    100 * ((1 - k) * S17 + k * S18) < 376 + 7 * k := by
  have first := OrderedRing.mul_pos
    (show 0 < 1 - k by grind) (show 0 < 94 - 25 * S17 by grind)
  have second := OrderedRing.mul_pos k_positive
    (show 0 < 383 - 100 * S18 by grind)
  grind

/-- Primitive squares, the finite pressure and both same-vector sums contradict. -/
theorem subcritical_pressure_contradiction (z k a b t S17 S18 : R)
    (z_lower : 1494 < 10000 * z) (z_upper : 2000 * z < 299)
    (k_lower : z ≤ k) (proper_square : k ^ 2 < z)
    (a_nonnegative : 0 ≤ a) (b_nonnegative : 0 ≤ b)
    (a_square : a ^ 2 = z - z ^ 2) (b_square : b ^ 2 = z - k ^ 2)
    (t_lower : 4141 < 1000 * t)
    (sum17_upper : 25 * S17 < 94) (sum18_upper : 100 * S18 < 383)
    (actual_finite_pressure : t - b ≤ (1 - k) * S17 + k * S18) : False := by
  obtain ⟨_z_positive, k_positive, k_upper⟩ :=
    subcritical_parameter_bounds z k z_lower z_upper k_lower proper_square
  have gap := subcritical_pressure_gap z k a b t z_lower z_upper k_lower
    a_nonnegative b_nonnegative a_square b_square t_lower
  have upper := subcritical_sum_pressure_upper k S17 S18 k_positive k_upper
    sum17_upper sum18_upper
  grind

end CoxeterProofs.B18TerminalPressureKernel

#print axioms CoxeterProofs.B18TerminalPressureKernel.core_quadratic_positive
#print axioms CoxeterProofs.B18TerminalPressureKernel.mixed_pressure_product
#print axioms CoxeterProofs.B18TerminalPressureKernel.supercritical_gap_positive
#print axioms CoxeterProofs.B18TerminalPressureKernel.supercritical_pressure_contradiction
#print axioms CoxeterProofs.B18TerminalPressureKernel.subcritical_base_root_bound
#print axioms CoxeterProofs.B18TerminalPressureKernel.subcritical_parameter_bounds
#print axioms CoxeterProofs.B18TerminalPressureKernel.subcritical_root_difference
#print axioms CoxeterProofs.B18TerminalPressureKernel.subcritical_pressure_gap
#print axioms CoxeterProofs.B18TerminalPressureKernel.subcritical_sum_pressure_upper
#print axioms CoxeterProofs.B18TerminalPressureKernel.subcritical_pressure_contradiction
