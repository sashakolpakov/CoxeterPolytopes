import Std

/-!
Prepared source only: no compilation or actual axiom output is claimed.

Internal ordered-ring algebra for the tight continuous D18 centre estimate.
The integer quadratic is strictly positive by an explicit square completion.
That positivity strengthens a supplied analytic lower bound, then gives both
K10/K11 paid-reference scalar contradictions.

External: the radical chord/tangent estimates, real parameter substitution,
actual future feasible core, complete cap inventory, positive reference
coefficients and payment, and interpretation of T as the original norm.
No geometric, field, polytope or dimension conclusion is formalized.
-/

namespace CoxeterProofs.D18CentreScalarKernel

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- Exact integer completion; its positive remainder is derived algebra. -/
theorem quadratic_completion (z : R) :
    387499200 * (96874800 * z ^ 2 - 1528000 * z + 6100) =
      (193749600 * z - 1528000) ^ 2 + 28961120000 := by
  grind only

/-- No square-sign premise is supplied by the application. -/
theorem quadratic_positive (z : R) :
    0 < 96874800 * z ^ 2 - 1528000 * z + 6100 := by
  have square_nonnegative :
      0 ≤ (193749600 * z - 1528000) ^ 2 := OrderedRing.sq_nonneg
  have scaled_positive :
      0 < 387499200 * (96874800 * z ^ 2 - 1528000 * z + 6100) := by
    rw [quadratic_completion]
    grind
  by_cases positive : 0 < 96874800 * z ^ 2 - 1528000 * z + 6100
  · exact positive
  · have nonpositive : 96874800 * z ^ 2 - 1528000 * z + 6100 ≤ 0 := by
      grind
    have scaled_nonpositive := OrderedRing.mul_nonpos_of_nonneg_of_nonpos
      (show (0 : R) ≤ 387499200 by grind) nonpositive
    grind

/-- The analytic radical bound is primitive; its uniform margin is derived. -/
theorem continuous_scalar_lower_bound (z F : R)
    (analytic_lower :
      34196500 - 1528000 * z + 96874800 * z ^ 2 ≤ 1000000 * F) :
    341904 < 10000 * F := by
  have positive := quadratic_positive z
  have cleared : (34190400 : R) < 1000000 * F := by grind
  grind

/-- K10: the funded reference inequality and primitive floors force T > 1. -/
theorem k10_margin (z F base f T : R)
    (analytic_lower :
      34196500 - 1528000 * z + 96874800 * z ^ 2 ≤ 1000000 * F)
    (base_lower : 87866 ≤ 100000 * base)
    (height_upper : 1000 * f ≤ 2731)
    (paid_reference : 49 * base + F - 14 * f + 10 ≤ 49 * T) :
    1 < T := by
  have scalar := continuous_scalar_lower_bound z F analytic_lower
  have exact_gap : (4901074 : R) < 4900000 * T := by grind
  grind

/-- K11: the separate coefficient bounds give the second exact margin. -/
theorem k11_margin (z F base f T : R)
    (analytic_lower :
      34196500 - 1528000 * z + 96874800 * z ^ 2 ≤ 1000000 * F)
    (base_lower : 90246 ≤ 100000 * base)
    (height_upper : 10000 * f ≤ 28854)
    (paid_reference : 49 * base + F - 14 * f + 11 ≤ 49 * T) :
    1 < T := by
  have scalar := continuous_scalar_lower_bound z F analytic_lower
  have exact_gap : (4901534 : R) < 4900000 * T := by grind
  grind

end CoxeterProofs.D18CentreScalarKernel

#print axioms CoxeterProofs.D18CentreScalarKernel.quadratic_completion
#print axioms CoxeterProofs.D18CentreScalarKernel.quadratic_positive
#print axioms CoxeterProofs.D18CentreScalarKernel.continuous_scalar_lower_bound
#print axioms CoxeterProofs.D18CentreScalarKernel.k10_margin
#print axioms CoxeterProofs.D18CentreScalarKernel.k11_margin
