import Std

/-!
SOURCE PREPARATION ONLY: uncompiled; no actual axiom audit is credited.

Conditional ordered-ring algebra for the frozen 470-line dotted-head proof.
Analytic Taylor envelopes, value/derivative bounds, geometric norm identities,
actual reference capacities and metric interpretation are theorem PARAMETERS.
No calculus, real radical function, polytope, complete catalogue or global
dimension theorem is formalized by this scalar file.
-/

namespace CoxeterProofs.D18UniformDottedHeadScalarKernel

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

private theorem derivative_square_bound (r : R)
    (lower : -8 ≤ r) (upper : r ≤ 8) : r ^ 2 ≤ 64 := by
  have product := OrderedRing.mul_nonneg
    (show 0 ≤ 8 - r by grind) (show 0 ≤ 8 + r by grind)
  grind

/-- The F completion is polynomial algebra, with no analytic premise. -/
theorem f_completion (z r : R) :
    14000000 * (3500000 * z ^ 2 + 1000 * r * z + 102) =
      (7000000 * z + 1000 * r) ^ 2 +
        1428000000 - 1000000 * r ^ 2 := by
  grind only

private theorem f_quadratic_positive (z r : R)
    (lower : -8 ≤ r) (upper : r ≤ 8) :
    0 < 3500000 * z ^ 2 + 1000 * r * z + 102 := by
  have squared := derivative_square_bound r lower upper
  have square : 0 ≤ (7000000 * z + 1000 * r) ^ 2 := OrderedRing.sq_nonneg
  have scaled :
      0 < 14000000 * (3500000 * z ^ 2 + 1000 * r * z + 102) := by
    rw [f_completion]
    grind
  grind

/--
Taylor/value/derivative statements are explicit external analytic premises.
Here z is u-2/3 and derivative is F'(2/3); no such interpretation is assumed
by the ordered-ring type. The exact uniform F>.99564 follows by completion.
-/
theorem f_uniform_seal (z value derivative F : R)
    (value_lower : 9956502 ≤ 10000000 * value)
    (derivative_lower : -8 ≤ 10000 * derivative)
    (derivative_upper : 10000 * derivative ≤ 8)
    (analytic_taylor :
      10000000 * value + 10000000 * derivative * z + 3500000 * z ^ 2 ≤
        10000000 * F) :
    99564 < 100000 * F := by
  have positive := f_quadratic_positive z (10000 * derivative)
    derivative_lower derivative_upper
  have cleared : (9956400 : R) < 10000000 * F := by grind
  grind

/-- The tight G completion is also independent of every analytic premise. -/
theorem g_completion (z r : R) :
    100000000000 *
        (25000000000 * z ^ 2 + 10000000 * r * z + 425728) =
      (50000000000 * z + 10000000 * r) ^ 2 +
        42572800000000000 - 100000000000000 * r ^ 2 := by
  grind only

private theorem g_quadratic_positive (z r : R)
    (lower : -7 ≤ r) (upper : r ≤ 7) :
    0 < 25000000000 * z ^ 2 + 10000000 * r * z + 425728 := by
  have product := OrderedRing.mul_nonneg
    (show 0 ≤ 7 - r by grind) (show 0 ≤ 7 + r by grind)
  have squared : r ^ 2 ≤ 49 := by grind
  have square : 0 ≤ (50000000000 * z + 10000000 * r) ^ 2 :=
    OrderedRing.sq_nonneg
  have scaled :
      0 < 100000000000 *
        (25000000000 * z ^ 2 + 10000000 * r * z + 425728) := by
    rw [g_completion]
    grind
  grind

/-- z is u-.642 in the application; the global analytic envelope is external. -/
theorem g_uniform_seal (z value derivative G : R)
    (value_lower : 99200425728 ≤ 100000000000 * value)
    (derivative_lower : -7 ≤ 10000 * derivative)
    (derivative_upper : 10000 * derivative ≤ 7)
    (analytic_taylor :
      100000000000 * value + 100000000000 * derivative * z +
          25000000000 * z ^ 2 ≤ 100000000000 * G) :
    992 < 1000 * G := by
  have positive := g_quadratic_positive z (10000 * derivative)
    derivative_lower derivative_upper
  have cleared : (99200000000 : R) < 100000000000 * G := by grind
  grind

/-- The two weak-helper debits fit their OWN capacities before pressure. -/
theorem weak_helper_anchor_affordability (u pay1 pay6 capacity1 capacity6 : R)
    (u_lower : 0 ≤ u) (u_upper : u ≤ 1)
    (first_payment : 250 * pay1 = 55 - 33 * u)
    (sixth_payment : 500 * pay6 = 11 * u)
    (first_capacity : 24672 ≤ 100000 * capacity1)
    (sixth_capacity : 1023 ≤ 10000 * capacity6) :
    0 ≤ pay1 ∧ pay1 < capacity1 ∧ 0 ≤ pay6 ∧ pay6 < capacity6 := by
  grind

/-- True normalized fork pressure and its squared edge law imply C<7.6. -/
theorem fork_cosh_cap (height C : R)
    (height_nonnegative : 0 ≤ height)
    (height_upper : 1000 * height ≤ 876)
    (cosh_lower : 1 ≤ C)
    (actual_edge_pressure : C - 1 ≤ height ^ 2 * (C + 1)) :
    5 * C < 38 := by
  have product := OrderedRing.mul_nonneg
    (show 0 ≤ 876 - 1000 * height by grind)
    (show 0 ≤ 876 + 1000 * height by grind)
  have squared : 1000000 * height ^ 2 ≤ 767376 := by grind
  have multiplier := OrderedRing.mul_nonneg
    (show 0 ≤ C + 1 by grind)
    (show 0 ≤ 767376 - 1000000 * height ^ 2 by grind)
  have cleared : 1000000 * (C - 1) ≤ 767376 * (C + 1) := by grind
  grind

/-- The own floor and exact original far-next length contradict that cap. -/
theorem far_next_impossible (height C own : R)
    (height_nonnegative : 0 ≤ height)
    (height_upper : 1000 * height ≤ 876)
    (cosh_lower : 1 ≤ C)
    (actual_edge_pressure : C - 1 ≤ height ^ 2 * (C + 1))
    (native_own_lower : 809 < 1000 * own)
    (actual_length : C = 9 * own + 1) : False := by
  have cap := fork_cosh_cap height C height_nonnegative height_upper
    cosh_lower actual_edge_pressure
  grind

/-- The exact Minkowski upper bound is external; its negative sign is derived. -/
theorem near_far_native_contradiction (nearOwn farOwn mutual : R)
    (far_nonnegative : 0 ≤ farOwn)
    (near_upper : 1000 * nearOwn ≤ 912)
    (far_upper : 1000 * farOwn ≤ 1301)
    (actual_acuteness : 0 ≤ mutual)
    (minkowski_upper :
      18 * mutual ≤ 9 * farOwn * (2 * nearOwn - 1) - 10) : False := by
  have product := OrderedRing.mul_nonneg far_nonnegative
    (show 0 ≤ 824 - 1000 * (2 * nearOwn - 1) by grind)
  have bound : 1000000 * farOwn * (2 * nearOwn - 1) ≤ 1072024 := by
    grind
  grind

/-- The m6 local-A18 graft uses ONE norm and its individually funded price. -/
theorem local_a18_six_tail_contradiction (body endpoint mean total : R)
    (body_lower : 6781845 ≤ 10000000 * body)
    (endpoint_lower : 7591 ≤ 10000 * endpoint)
    (mean_upper : 10000 * mean ≤ 2352)
    (paid_reference : 10 * body + 12 * (endpoint - mean) - 3 ≤ 10 * total)
    (interior_norm : total < 1) : False := by
  grind

/-- The m5 alternative retains its own body floor and mean. -/
theorem local_a18_five_tail_contradiction (body endpoint mean total : R)
    (body_lower : 69794985 ≤ 100000000 * body)
    (endpoint_lower : 7591 ≤ 10000 * endpoint)
    (mean_upper : 1000 * mean ≤ 261)
    (paid_reference : 125 * body + 150 * (endpoint - mean) - 36 ≤ 125 * total)
    (interior_norm : total < 1) : False := by
  grind

end CoxeterProofs.D18UniformDottedHeadScalarKernel

#print axioms CoxeterProofs.D18UniformDottedHeadScalarKernel.f_uniform_seal
#print axioms CoxeterProofs.D18UniformDottedHeadScalarKernel.g_uniform_seal
#print axioms CoxeterProofs.D18UniformDottedHeadScalarKernel.weak_helper_anchor_affordability
#print axioms CoxeterProofs.D18UniformDottedHeadScalarKernel.far_next_impossible
#print axioms CoxeterProofs.D18UniformDottedHeadScalarKernel.near_far_native_contradiction
#print axioms CoxeterProofs.D18UniformDottedHeadScalarKernel.local_a18_six_tail_contradiction
#print axioms CoxeterProofs.D18UniformDottedHeadScalarKernel.local_a18_five_tail_contradiction
