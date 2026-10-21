import CoxeterProofs.SummedActualFootKernel

/-!
Source preparation: compilation and actual axiom outputs are recorded separately.

The continuous isolated-foot estimate used in the shorter arithmetic A proof.
The circle identity and its maximum are derived, then a two-interval scalar
criterion absorbs every loss into its positive dual weight. The final theorem
combines those local scalar obligations with the already derived finite-coordinate
summed-foot inequality. It allows adjacent feet and zero reference weights.

Applications still supply actual geometry/feasibility, complete cap inventory,
the reference and its coefficient bounds, and the three concrete rational
parameter comparisons. No dimension, field classification, or scalar conclusion
for a concrete polytope is assumed or formalized here.
-/

namespace CoxeterProofs.ContinuousLossAbsorptionKernel

open Lean.Grind
open CoxeterProofs.RightAngledCubicScalarKernel
open CoxeterProofs.RightAngledProjectionHeightKernel
open CoxeterProofs.SummedActualFootKernel

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

theorem circle_identity (v s : R) (root : s*s = v-v*v) :
    (2*v-1)*(2*v-1) + 4*s*s = 1 := by
  grind

theorem loss_square_identity (a v s : R) (root : s*s = v-v*v) :
    (2*(a*v+s)-a)*(2*(a*v+s)-a) +
      (2*v-1-2*a*s)*(2*v-1-2*a*s) = 1+a*a := by
  have circle := circle_identity v s root
  calc
    _ = (1+a*a)*((2*v-1)*(2*v-1)+4*s*s) := by grind
    _ = _ := by rw [circle]; grind

/-- A maximum on the whole continuous interval, derived without calculus. -/
theorem loss_linear_bound (a v s r : R) (root : s*s = v-v*v)
    (r_nonnegative : 0 ≤ r) (radius : 1+a*a ≤ r*r) :
    2*(a*v+s) ≤ a+r := by
  have identity := loss_square_identity a v s root
  have square_nonnegative : 0 ≤ (2*v-1-2*a*s)^2 := OrderedRing.sq_nonneg
  simp only [Semiring.pow_two] at square_nonnegative
  by_cases conclusion : 2*(a*v+s) ≤ a+r
  · exact conclusion
  · have gap : 0 < 2*(a*v+s)-a-r := by grind
    have total : 0 < 2*(a*v+s)-a+r := by grind
    have product := OrderedRing.mul_pos gap total
    grind

theorem half_circle_bound (v s : R) (root : s*s = v-v*v)
    (s_nonnegative : 0 ≤ s) : 2*s ≤ 1 := by
  have circle := circle_identity v s root
  have square_nonnegative : 0 ≤ (2*v-1)^2 := OrderedRing.sq_nonneg
  simp only [Semiring.pow_two] at square_nonnegative
  by_cases conclusion : 2*s ≤ 1
  · exact conclusion
  · have gap : 0 < 2*s-1 := by grind
    have total : 0 < 2*s+1 := by grind
    have product := OrderedRing.mul_pos gap total
    grind

/-- Both scalar tests are independent of the unknown continuous parameter v.
All factors of two/four are cleared, so no division or field instance is used. -/
theorem adjusted_loss_absorption (a b ell r beta v s : R)
    (a_nonnegative : 0 ≤ a) (b_nonnegative : 0 ≤ b)
    (weight : b ≤ beta) (v_nonnegative : 0 ≤ v) (v_upper : v ≤ 1)
    (s_nonnegative : 0 ≤ s) (root : s*s = v-v*v)
    (r_nonnegative : 0 ≤ r) (radius : 1+a*a ≤ r*r)
    (small : a+1 < 2*ell) (large : 4*(a+r-2*ell) < 3*b) :
    2*(a*v+s-ell) < beta*(2*v-v*v) := by
  have beta_nonnegative : 0 ≤ beta := by grind
  have second_factor : 0 ≤ 2-v := by grind
  have quadratic_product := OrderedRing.mul_nonneg v_nonnegative second_factor
  have quadratic_nonnegative : 0 ≤ 2*v-v*v := by grind
  by_cases first_half : 2*v ≤ 1
  · have factor : 0 ≤ 1-2*v := by grind
    have a_product := OrderedRing.mul_nonneg a_nonnegative factor
    have s_bound := half_circle_bound v s root s_nonnegative
    have target_nonnegative := OrderedRing.mul_nonneg beta_nonnegative quadratic_nonnegative
    grind
  · have left_factor : 0 ≤ 2*v-1 := by grind
    have right_factor : 0 ≤ 3-2*v := by grind
    have lower_product := OrderedRing.mul_nonneg left_factor right_factor
    have lower : 0 ≤ 4*(2*v-v*v)-3 := by grind
    have weighted_lower := OrderedRing.mul_nonneg b_nonnegative lower
    have difference : 0 ≤ beta-b := by grind
    have monotone_weight := OrderedRing.mul_nonneg difference quadratic_nonnegative
    have upper := loss_linear_bound a v s r root r_nonnegative radius
    grind

theorem absorption_energy_term (beta v loss ell : R)
    (weight : 0 ≤ beta)
    (absorbed : 2*(loss-ell) ≤ beta*(2*v-v*v)) :
    0 ≤ 2*beta*(ell-loss)+(2*v-v*v)*(beta*beta) := by
  have residual : 0 ≤ beta*(2*v-v*v)-2*(loss-ell) := by grind
  have product := OrderedRing.mul_nonneg weight residual
  grind

/-- Geometry supplies each actual-foot pressure; local scalar absorption now
implies the reference energy. Omitted directions require no loss estimate. -/
theorem absorbed_reference_energy (n d : Nat) (beta v loss : Nat → R)
    (ell : R) (omega normal : Nat → Nat → R) (x : Nat → R)
    (weight : ∀ i, i < n → 0 ≤ beta i)
    (parameter : ∀ i, i < n → 0 ≤ v i)
    (unit : ∀ i, i < n → coordinateDot d (normal i) (normal i) = 1)
    (obtuse : ∀ i, i < n → ∀ j, j < n → i ≠ j →
      coordinateDot d (normal i) (normal j) ≤ 0)
    (dual : ∀ i, i < n →
      coordinateDot d (weightedSum n beta omega) (normal i) = beta i)
    (foot : ∀ i, i < n →
      coordinateDot d (weightedSum n beta omega) (omega i) + ell - loss i ≤
        coordinateDot d x (fun c => omega i c - v i * normal i c))
    (absorbed : ∀ i, i < n → 0 < beta i →
      2*(loss i-ell) ≤ beta i*(2*v i-v i*v i)) :
    coordinateDot d (weightedSum n beta omega) (weightedSum n beta omega) ≤
      coordinateDot d x x := by
  let w := weightedSum n beta omega
  let T := fun i => coordinateDot d w (omega i) + ell - loss i
  have energy := summed_actual_foot_energy n d beta v T omega normal x
    weight parameter unit obtuse dual foot
  have nonnegative : 0 ≤ sumN n (fun i =>
      2*beta i*(T i-coordinateDot d w (omega i)) +
        (2*v i-v i*v i)*(beta i*beta i)) := by
    apply sumN_nonneg
    intro i hi
    by_cases positive : 0 < beta i
    · have row := absorption_energy_term (beta i) (v i) (loss i) ell
        (weight i hi) (absorbed i hi positive)
      change 0 ≤ 2*beta i*((coordinateDot d w (omega i)+ell-loss i)-
        coordinateDot d w (omega i))+(2*v i-v i*v i)*(beta i*beta i)
      grind
    · have zero : beta i = 0 := by have nonnegative := weight i hi; grind
      simpa only [zero, Semiring.zero_mul, Semiring.mul_zero, Semiring.add_zero]
        using (Std.le_refl (0 : R))
  change coordinateDot d w w + sumN n (fun i =>
    2*beta i*(T i-coordinateDot d w (omega i)) +
      (2*v i-v i*v i)*(beta i*beta i)) ≤ coordinateDot d x x at energy
  grind

end CoxeterProofs.ContinuousLossAbsorptionKernel

#print axioms CoxeterProofs.ContinuousLossAbsorptionKernel.circle_identity
#print axioms CoxeterProofs.ContinuousLossAbsorptionKernel.loss_square_identity
#print axioms CoxeterProofs.ContinuousLossAbsorptionKernel.loss_linear_bound
#print axioms CoxeterProofs.ContinuousLossAbsorptionKernel.half_circle_bound
#print axioms CoxeterProofs.ContinuousLossAbsorptionKernel.adjusted_loss_absorption
#print axioms CoxeterProofs.ContinuousLossAbsorptionKernel.absorption_energy_term
#print axioms CoxeterProofs.ContinuousLossAbsorptionKernel.absorbed_reference_energy
