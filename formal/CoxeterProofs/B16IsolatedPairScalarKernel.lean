import Std

/-!
PREPARED SOURCE ONLY. No compilation or actual axiom audit is credited.

Internal logic for the isolated-pair whole-B16 arguments: odd valuation parity,
the complete small-conjugate cubic argument, and the two case-specific prefix squares.
Number-field/ideal facts, integer discreteness, the actual cubic instance,
edge pressures and the original vector's norm identity remain external.
Metric scalars are in an ordered ring; only valuations use actual integers.
-/

namespace CoxeterProofs.B16IsolatedPairScalarKernel

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- Valuations really are integers; an odd dyadic valuation cannot be a square. -/
theorem odd_dyadic_valuation_contradiction (e v k : Int)
    (valuation_identity : v = (e - 1) + 1 - 2 * e - (e - 1))
    (square_even : v = 2 * k) : False := by
  omega

/-- An exact positive product gives the strict cube-versus-square gap. -/
theorem unit_interval_cube_gap (m : R) (positive : 0 < m) (upper : m < 1) :
    0 < m ^ 2 - m ^ 3 := by
  have square_positive := OrderedRing.mul_pos positive positive
  have product := OrderedRing.mul_pos square_positive
    (show 0 < 1 - m by grind)
  grind

/-- This entire quadratic is negative, independently of a sampled point. -/
theorem six_cost_quadratic_negative (m : R) :
    -5 * m ^ 2 + 6 * m - 2 < 0 := by
  have square : 0 ≤ (5 * m - 3) ^ 2 := OrderedRing.sq_nonneg
  have identity : 5 * (-5 * m ^ 2 + 6 * m - 2) =
      -(5 * m - 3) ^ 2 - 1 := by grind only
  grind

/-- The second norm branch has its own exact negative square. -/
theorem four_cost_quadratic_negative (m : R) :
    -5 * m ^ 2 + 4 * m - 1 < 0 := by
  have square : 0 ≤ (5 * m - 2) ^ 2 := OrderedRing.sq_nonneg
  have identity : 5 * (-5 * m ^ 2 + 4 * m - 1) =
      -(5 * m - 2) ^ 2 - 1 := by grind only
  grind

/-- Ordinary data, never an instance: the actual field interpretation is external. -/
structure CubicData (R : Type u) [CommRing R] [LE R] [LT R]
    [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R] where
  w : R
  a : R
  b : R
  m : R
  S : R
  A : R
  B : R
  large_root : 1 < w
  first_small_positive : 0 < a
  small_order : a < b
  second_small_upper : b < 1
  midpoint : 2 * m = a + b
  trace_identity : S = w + a + b
  norm_identity : A = w * a * b
  shifted_norm_identity : B = (w - 1) * (1 - a) * (1 - b)
  norm_A_discrete : A = 1 ∨ 2 ≤ A
  norm_B_discrete : B = 1 ∨ 2 ≤ B
  trace_discrete : S ≤ 5 ∨ S = 6 ∨ 7 ≤ S

def polynomial (d : CubicData R) (x : R) : R :=
  x ^ 3 - d.S * x ^ 2 + (d.A + d.S - d.B - 1) * x - d.A

/-- The polynomial factorization is derived from primitive root identities. -/
theorem polynomial_factorization (d : CubicData R) (x : R) :
    polynomial d x = (x - d.w) * (x - d.a) * (x - d.b) := by
  have trace_identity := d.trace_identity
  have norm_identity := d.norm_identity
  have shifted_norm_identity := d.shifted_norm_identity
  unfold polynomial
  grind only

/-- No positive interval is assumed: the midpoint and factorization imply it. -/
theorem midpoint_polynomial_positive (d : CubicData R) :
    0 < d.m ∧ d.m < 1 ∧ 0 < polynomial d d.m := by
  have large_root := d.large_root
  have first_small_positive := d.first_small_positive
  have small_order := d.small_order
  have second_small_upper := d.second_small_upper
  have midpoint := d.midpoint
  have m_positive : 0 < d.m := by grind
  have m_upper : d.m < 1 := by grind
  have first := OrderedRing.mul_pos
    (show 0 < d.w - d.m by grind)
    (show 0 < d.m - d.a by grind)
  have product := OrderedRing.mul_pos first
    (show 0 < d.b - d.m by grind)
  have factorization := polynomial_factorization d d.m
  have positive : 0 < polynomial d d.m := by grind
  exact ⟨m_positive, m_upper, positive⟩

/-- A>=2 contradicts the full positive midpoint under the trace bound. -/
theorem norm_A_large_contradiction (d : CubicData R)
    (trace_upper : d.S ≤ 6) (large_A : 2 ≤ d.A) : False := by
  obtain ⟨m_positive, m_upper, poly_positive⟩ := midpoint_polynomial_positive d
  have discrete_B := d.norm_B_discrete
  have lower_B : 1 ≤ d.B := by grind
  have interval := OrderedRing.mul_pos m_positive
    (show 0 < 1 - d.m by grind)
  have trace_cost := OrderedRing.mul_nonneg
    (show 0 ≤ 6 - d.S by grind) (show 0 ≤ d.m * (1 - d.m) by grind)
  have norm_cost := OrderedRing.mul_nonneg
    (show 0 ≤ d.A - 2 by grind) (show 0 ≤ 1 - d.m by grind)
  have shifted_cost := OrderedRing.mul_nonneg
    (show 0 ≤ d.B - 1 by grind) (show 0 ≤ d.m by grind)
  have gap := unit_interval_cube_gap d.m m_positive m_upper
  have negative := six_cost_quadratic_negative (R := R) d.m
  have identity :
      (d.m ^ 3 - 6 * d.m ^ 2 + 6 * d.m - 2) - polynomial d d.m =
        (6 - d.S) * (d.m * (1 - d.m)) +
          (d.A - 2) * (1 - d.m) + (d.B - 1) * d.m := by
    unfold polynomial
    grind only
  grind

/-- B>=2 is excluded by the second full negative quadratic branch. -/
theorem norm_B_large_contradiction (d : CubicData R)
    (trace_upper : d.S ≤ 6) (large_B : 2 ≤ d.B) : False := by
  obtain ⟨m_positive, m_upper, poly_positive⟩ := midpoint_polynomial_positive d
  have discrete_A := d.norm_A_discrete
  have lower_A : 1 ≤ d.A := by grind
  have interval := OrderedRing.mul_pos m_positive
    (show 0 < 1 - d.m by grind)
  have trace_cost := OrderedRing.mul_nonneg
    (show 0 ≤ 6 - d.S by grind) (show 0 ≤ d.m * (1 - d.m) by grind)
  have norm_cost := OrderedRing.mul_nonneg
    (show 0 ≤ d.A - 1 by grind) (show 0 ≤ 1 - d.m by grind)
  have shifted_cost := OrderedRing.mul_nonneg
    (show 0 ≤ d.B - 2 by grind) (show 0 ≤ d.m by grind)
  have gap := unit_interval_cube_gap d.m m_positive m_upper
  have negative := four_cost_quadratic_negative (R := R) d.m
  have identity :
      (d.m ^ 3 - 6 * d.m ^ 2 + 4 * d.m - 1) - polynomial d d.m =
        (6 - d.S) * (d.m * (1 - d.m)) +
          (d.A - 1) * (1 - d.m) + (d.B - 2) * d.m := by
    unfold polynomial
    grind only
  grind

/-- When both norms equal one, every trace at most five is impossible. -/
theorem low_trace_contradiction (d : CubicData R)
    (norm_A_one : d.A = 1) (norm_B_one : d.B = 1)
    (trace_upper : d.S ≤ 5) : False := by
  obtain ⟨m_positive, m_upper, poly_positive⟩ := midpoint_polynomial_positive d
  have interval := OrderedRing.mul_pos m_positive
    (show 0 < 1 - d.m by grind)
  have trace_cost := OrderedRing.mul_nonneg
    (show 0 ≤ 5 - d.S by grind) (show 0 ≤ d.m * (1 - d.m) by grind)
  have square : 0 ≤ (2 * d.m - 1) ^ 2 := OrderedRing.sq_nonneg
  have gap := unit_interval_cube_gap d.m m_positive m_upper
  have identity : polynomial d d.m + (2 * d.m - 1) ^ 2 +
      (5 - d.S) * (d.m * (1 - d.m)) = d.m ^ 3 - d.m ^ 2 := by
    unfold polynomial
    grind only
  grind

/-- The entire cubic argument, including the trace-six endpoint, gives w>5. -/
theorem cubic_physical_root_gt_five (d : CubicData R) : 5 < d.w := by
  by_cases target : 5 < d.w
  · exact target
  · have w_upper : d.w ≤ 5 := by grind
    have small_order := d.small_order
    have small_upper := d.second_small_upper
    have trace_identity := d.trace_identity
    have trace_alternatives := d.trace_discrete
    have trace_lt_seven : d.S < 7 := by grind
    have trace_upper : d.S ≤ 6 := by grind
    have norm_A_one : d.A = 1 := by
      cases d.norm_A_discrete with
      | inl one => exact one
      | inr large => exact False.elim (norm_A_large_contradiction d trace_upper large)
    have norm_B_one : d.B = 1 := by
      cases d.norm_B_discrete with
      | inl one => exact one
      | inr large => exact False.elim (norm_B_large_contradiction d trace_upper large)
    have trace_not_low : ¬d.S ≤ 5 := by
      intro low
      exact low_trace_contradiction d norm_A_one norm_B_one low
    have trace_six : d.S = 6 := by grind
    have first := OrderedRing.mul_nonneg
      (show 0 ≤ 5 - d.w by grind) (show 0 ≤ 5 - d.a by grind)
    have product := OrderedRing.mul_nonneg first (show 0 ≤ 5 - d.b by grind)
    have factorization := polynomial_factorization d (5 : R)
    have endpoint : polynomial d (5 : R) = -1 := by
      unfold polynomial
      grind only
    grind

/-- The literal scaled vector norm is computed rather than assumed as a gain. -/
theorem cubic_literal_prefix_cost :
    (123600 : R) ^ 2 + 5 * (56006 : R) ^ 2 + 10 * (30189 : R) ^ 2 =
      (40074077390 : R) := by
  grind only

/-- Primitive actual pressures and ONE raw norm square give the exact contradiction. -/
theorem cubic_prefix_reference_contradiction (x1 S6 S16 T energy : R)
    (head_pressure : 618 < 1000 * x1)
    (prefix_pressure : 201815 < 100000 * S6)
    (short_pressure : 35276 < 10000 * S16)
    (dot_identity : T = 67594 * x1 + 25817 * S6 + 30189 * S16)
    (one_square : 0 ≤ 40000000000 * energy - 400000 * T +
      (123600 : R) ^ 2 + 5 * (56006 : R) ^ 2 + 10 * (30189 : R) ^ 2)
    (energy_upper : energy < 1) : False := by
  have head := OrderedRing.mul_pos
    (show 0 < (27037600 : R) by grind)
    (show 0 < 1000 * x1 - 618 by grind)
  have middle_gain := OrderedRing.mul_pos
    (show 0 < (103268 : R) by grind)
    (show 0 < 100000 * S6 - 201815 by grind)
  have short := OrderedRing.mul_pos
    (show 0 < (1207560 : R) by grind)
    (show 0 < 10000 * S16 - 35276 by grind)
  have pressure_identity :
      400000 * T - 80148154780 =
        27037600 * (1000 * x1 - 618) +
          103268 * (100000 * S6 - 201815) +
            1207560 * (10000 * S16 - 35276) := by
    grind only
  have norm := cubic_literal_prefix_cost (R := R)
  grind

/-- The surviving quadratic pair has a different, literally computed vector norm. -/
theorem quadratic_literal_prefix_cost :
    (122000 : R) ^ 2 + 5 * (56326 : R) ^ 2 + 10 * (30837 : R) ^ 2 =
      (40256297070 : R) := by
  grind only

/-- The exact quadratic-pair pressures contradict ONE original-vector norm. -/
theorem quadratic_prefix_reference_contradiction (x1 S6 S16 T energy : R)
    (head_pressure : 61 < 100 * x1)
    (prefix_pressure : 201815 < 100000 * S6)
    (short_pressure : 356 < 100 * S16)
    (dot_identity : T = 65674 * x1 + 25489 * S6 + 30837 * S16)
    (one_square : 0 ≤ 40000000000 * energy - 400000 * T +
      (122000 : R) ^ 2 + 5 * (56326 : R) ^ 2 + 10 * (30837 : R) ^ 2)
    (energy_upper : energy < 1) : False := by
  have head := OrderedRing.mul_pos
    (show 0 < (262696000 : R) by grind)
    (show 0 < 100 * x1 - 61 by grind)
  have middle_gain := OrderedRing.mul_pos
    (show 0 < (101956 : R) by grind)
    (show 0 < 100000 * S6 - 201815 by grind)
  have short := OrderedRing.mul_pos
    (show 0 < (123348000 : R) by grind)
    (show 0 < 100 * S16 - 356 by grind)
  have pressure_identity :
      400000 * T - 80512594140 =
        262696000 * (100 * x1 - 61) +
          101956 * (100000 * S6 - 201815) +
            123348000 * (100 * S16 - 356) := by
    grind only
  have norm := quadratic_literal_prefix_cost (R := R)
  grind

end CoxeterProofs.B16IsolatedPairScalarKernel

#print axioms CoxeterProofs.B16IsolatedPairScalarKernel.odd_dyadic_valuation_contradiction
#print axioms CoxeterProofs.B16IsolatedPairScalarKernel.unit_interval_cube_gap
#print axioms CoxeterProofs.B16IsolatedPairScalarKernel.polynomial_factorization
#print axioms CoxeterProofs.B16IsolatedPairScalarKernel.midpoint_polynomial_positive
#print axioms CoxeterProofs.B16IsolatedPairScalarKernel.norm_A_large_contradiction
#print axioms CoxeterProofs.B16IsolatedPairScalarKernel.norm_B_large_contradiction
#print axioms CoxeterProofs.B16IsolatedPairScalarKernel.low_trace_contradiction
#print axioms CoxeterProofs.B16IsolatedPairScalarKernel.cubic_physical_root_gt_five
#print axioms CoxeterProofs.B16IsolatedPairScalarKernel.cubic_literal_prefix_cost
#print axioms CoxeterProofs.B16IsolatedPairScalarKernel.cubic_prefix_reference_contradiction
#print axioms CoxeterProofs.B16IsolatedPairScalarKernel.quadratic_literal_prefix_cost
#print axioms CoxeterProofs.B16IsolatedPairScalarKernel.quadratic_prefix_reference_contradiction
