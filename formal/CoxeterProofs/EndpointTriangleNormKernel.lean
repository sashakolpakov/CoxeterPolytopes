import Std

/-!
Internal algebra and norm-budget logic for the A16 endpoint argument.
Compilation and the separate actual axiom audit are recorded in the logic map.
Source preparation alone is not verification.

The cycle expansion and integral Bezout expression are derived. Finite
products and the two-factor contraction are also derived. Actual walls,
first-type cyclic integrality, conjugate Schur bounds, reference energies,
number-field norm interpretation and integral-norm lower bounds remain
application inputs. No geometric forest or ambient bound is proved here.
-/

namespace CoxeterProofs.EndpointTriangleNormKernel

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- Cancellation of the two- and three-cycle terms is proved algebraically. -/
theorem cycle_expansion (A a b r lam U h rho v : R)
    (a_square : a ^ 2 = h) (r_square : r ^ 2 = 1 + h)
    (time_square : 17 * A ^ 2 = 8 * U + 17 * b ^ 2 - 17)
    (lambda_identity : lam = A * a - b * r)
    (cycle_identity : rho = 4 * lam ^ 2 + 8 * b * r * lam + 4 * b ^ 2)
    (v_identity : v = 4 * h) :
    17 * rho = (8 * U - 17) * v := by
  have cancellation :
      4 * (A * a - b * r) ^ 2 + 8 * b * r * (A * a - b * r) + 4 * b ^ 2 =
        4 * (A ^ 2 * a ^ 2 - b ^ 2 * r ^ 2 + b ^ 2) := by
    grind only
  rw [lambda_identity, cancellation, a_square, r_square] at cycle_identity
  rw [cycle_identity, v_identity]
  have scaled := congrArg (fun z : R => 4 * h * z) time_square
  grind only

/-- This explicit ring expression supplies the integral theta witness. -/
theorem bezout_triangle_identity (U v rho : R)
    (cycle : 17 * rho = (8 * U - 17) * v) :
    17 * (15 * rho - (7 * U - 15) * v) = U * v := by
  have scaled := congrArg (fun z : R => 15 * z) cycle
  grind

theorem physical_I_cap (U v theta : R)
    (U_positive : 0 < U) (v_positive : 0 < v)
    (U_cap : 2 * U < 17) (v_cap : v < 2)
    (identity : 17 * theta = U * v) : 0 < theta ∧ theta < 1 := by
  have positive := OrderedRing.mul_pos U_positive v_positive
  have upper := OrderedRing.mul_lt_mul_of_pos_left v_cap U_positive
  constructor <;> grind

/-- The other embeddings have -1/2 < theta < 0, without division notation. -/
theorem conjugate_cap (U v theta : R)
    (U_positive : 0 < U) (U_cap : 8 * U < 17)
    (v_lower : -4 < v) (v_negative : v < 0)
    (identity : 17 * theta = U * v) : -1 < 2 * theta ∧ theta < 0 := by
  have negative := OrderedRing.mul_pos U_positive (show 0 < -v by grind)
  have lower := OrderedRing.mul_lt_mul_of_pos_left v_lower U_positive
  constructor <;> grind

/-- The larger physical range needed for the ordinary-A1 extension. -/
theorem physical_extended_cap (U v theta : R)
    (U_positive : 0 < U) (v_positive : 0 < v)
    (U_cap : 2 * U < 17) (v_cap : 5 * v < 14)
    (identity : 17 * theta = U * v) : 0 < theta ∧ 5 * theta < 7 := by
  have positive := OrderedRing.mul_pos U_positive v_positive
  have upper := OrderedRing.mul_lt_mul_of_pos_left v_cap U_positive
  constructor <;> grind

def productN : Nat → (Nat → R) → R
  | 0, _ => 1
  | n + 1, f => productN n f * f n

/-- Every remaining embedding factor may be kept; none is silently dropped. -/
theorem product_cap (n : Nat) (f : Nat → R)
    (nonnegative : ∀ i, i < n → 0 ≤ f i)
    (bounded : ∀ i, i < n → f i ≤ 1) :
    0 ≤ productN n f ∧ productN n f ≤ 1 := by
  induction n with
  | zero => simp only [productN]; constructor <;> grind
  | succ n ih =>
    have previous := ih (fun i hi => nonnegative i (by omega))
      (fun i hi => bounded i (by omega))
    have last_nonnegative := nonnegative n (by omega)
    have last_bounded := bounded n (by omega)
    have positive := OrderedRing.mul_nonneg previous.1 last_nonnegative
    have upper := OrderedRing.mul_le_mul_of_nonneg_left last_bounded previous.1
    constructor
    · exact positive
    · change productN n f * f n ≤ 1
      grind

theorem norm_small_factor_impossible (n : Nat) (f : Nat → R) (a norm : R)
    (nonnegative : ∀ i, i < n → 0 ≤ f i)
    (bounded : ∀ i, i < n → f i ≤ 1)
    (a_positive : 0 < a) (a_small : a < 1)
    (norm_identity : norm = a * productN n f)
    (integral_norm_lower : 1 ≤ norm) : False := by
  have tail := product_cap n f nonnegative bounded
  have upper := OrderedRing.mul_le_mul_of_nonneg_left tail.2
    (show 0 ≤ a by grind)
  grind

/-- One nonphysical factor contracts a physical factor smaller than 7/5. -/
theorem expanding_pair_contracts (a b : R)
    (a_positive : 0 < a) (b_positive : 0 < b)
    (a_cap : 5 * a < 7) (b_cap : 2 * b < 1) :
    0 < a * b ∧ a * b < 1 := by
  have positive := OrderedRing.mul_pos a_positive b_positive
  have upper := OrderedRing.mul_lt_mul_of_pos_left b_cap a_positive
  constructor <;> grind

theorem norm_expanding_factor_impossible (n : Nat) (f : Nat → R)
    (a b norm : R)
    (nonnegative : ∀ i, i < n → 0 ≤ f i)
    (bounded : ∀ i, i < n → f i ≤ 1)
    (a_positive : 0 < a) (b_positive : 0 < b)
    (a_cap : 5 * a < 7) (b_cap : 2 * b < 1)
    (norm_identity : norm = (a * b) * productN n f)
    (integral_norm_lower : 1 ≤ norm) : False := by
  have pair := expanding_pair_contracts a b a_positive b_positive a_cap b_cap
  exact norm_small_factor_impossible n f (a * b) norm nonnegative bounded
    pair.1 pair.2 norm_identity integral_norm_lower

/-- Degree one is handled directly, without assuming Meyer or quadraticity. -/
theorem rational_positive_bounded_product_not_multiple (a b : Nat)
    (a_positive : 0 < a) (b_positive : 0 < b)
    (a_cap : a ≤ 8) (b_cap : b ≤ 2) (divisible : 17 ∣ a * b) : False := by
  have positive := Nat.mul_pos a_positive b_positive
  have cap : a * b ≤ 16 := Nat.mul_le_mul a_cap b_cap
  obtain ⟨k, equation⟩ := divisible
  omega

end CoxeterProofs.EndpointTriangleNormKernel

#print axioms CoxeterProofs.EndpointTriangleNormKernel.cycle_expansion
#print axioms CoxeterProofs.EndpointTriangleNormKernel.bezout_triangle_identity
#print axioms CoxeterProofs.EndpointTriangleNormKernel.physical_I_cap
#print axioms CoxeterProofs.EndpointTriangleNormKernel.conjugate_cap
#print axioms CoxeterProofs.EndpointTriangleNormKernel.physical_extended_cap
#print axioms CoxeterProofs.EndpointTriangleNormKernel.product_cap
#print axioms CoxeterProofs.EndpointTriangleNormKernel.norm_small_factor_impossible
#print axioms CoxeterProofs.EndpointTriangleNormKernel.expanding_pair_contracts
#print axioms CoxeterProofs.EndpointTriangleNormKernel.norm_expanding_factor_impossible
#print axioms CoxeterProofs.EndpointTriangleNormKernel.rational_positive_bounded_product_not_multiple
