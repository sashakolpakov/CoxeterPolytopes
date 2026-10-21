import CoxeterProofs.SummedActualFootKernel

/-!
Source preparation; compilation and actual axiom outputs are recorded separately.

A reference is formed from primitive nonnegative weights. Only positive-weight
directions require a pressure floor; omitted directions need none. Its pressure,
energy lower bound and contradiction are derived from finite coordinate sums.
Cleared two/three-slope averaging identities and the nonflat graft square are
also derived. Applications still supply actual geometry, weights, available
floors, coordinate realizations and concrete scalar bounds. No field, census,
polytope, or dimension conclusion is formalized here.
-/

namespace CoxeterProofs.OmittedPrefixReferenceKernel

open Lean.Grind
open CoxeterProofs.RightAngledCubicScalarKernel
open CoxeterProofs.RightAngledProjectionHeightKernel
open CoxeterProofs.SummedActualFootKernel

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- Zero weights need no floor; the reference pressure is proved by summing. -/
theorem supported_reference_pressure (n d : Nat) (beta : Nat → R)
    (omega : Nat → Nat → R) (x : Nat → R)
    (nonnegative : ∀ i, i < n → 0 ≤ beta i)
    (floor : ∀ i, i < n → 0 < beta i →
      coordinateDot d (weightedSum n beta omega) (omega i) ≤
        coordinateDot d x (omega i)) :
    coordinateDot d (weightedSum n beta omega) (weightedSum n beta omega) ≤
      coordinateDot d x (weightedSum n beta omega) := by
  rw [reference_energy_identity n d beta omega,
    dot_weighted_right n d beta omega x]
  apply sumN_mono
  intro i hi
  by_cases active : 0 < beta i
  · exact OrderedRing.mul_le_mul_of_nonneg_left
      (floor i hi active) (nonnegative i hi)
  · have weight_zero : beta i = 0 := by
      have weight_nonnegative := nonnegative i hi
      grind
    simpa only [weight_zero, Semiring.zero_mul] using (Std.le_refl (0 : R))

/-- Completing one coordinate square turns the derived pressure into energy. -/
theorem supported_reference_energy_lower (n d : Nat) (beta : Nat → R)
    (omega : Nat → Nat → R) (x : Nat → R)
    (nonnegative : ∀ i, i < n → 0 ≤ beta i)
    (floor : ∀ i, i < n → 0 < beta i →
      coordinateDot d (weightedSum n beta omega) (omega i) ≤
        coordinateDot d x (omega i)) :
    coordinateDot d (weightedSum n beta omega) (weightedSum n beta omega) ≤
      coordinateDot d x x := by
  have pressure := supported_reference_pressure n d beta omega x nonnegative floor
  have square := coordinate_energy_lower d x (weightedSum n beta omega)
  grind

theorem omitted_reference_contradiction (n d : Nat) (beta : Nat → R)
    (omega : Nat → Nat → R) (x : Nat → R) (E : R)
    (nonnegative : ∀ i, i < n → 0 ≤ beta i)
    (floor : ∀ i, i < n → 0 < beta i →
      coordinateDot d (weightedSum n beta omega) (omega i) ≤
        coordinateDot d x (omega i))
    (energy : E = coordinateDot d (weightedSum n beta omega)
      (weightedSum n beta omega))
    (large : 1 < E) (inside : coordinateDot d x x < 1) : False := by
  have lower := supported_reference_energy_lower n d beta omega x nonnegative floor
  grind

/-- Exact energy debit; no division or supplied variance identity is used. -/
theorem two_slope_averaging_energy (a b m : R) (average : 2 * m = a + b) :
    2 * (a*a + b*b - 2*m*m) = (a-b)*(a-b) := by
  calc
    _ = 2 * (a*a + b*b) - (2*m)*(2*m) := by grind
    _ = 2 * (a*a + b*b) - (a+b)*(a+b) := by simp only [average]
    _ = _ := by grind

theorem two_slope_averaging_nonincrease (a b m : R)
    (average : 2 * m = a + b) : 2*m*m ≤ a*a + b*b := by
  have identity := two_slope_averaging_energy a b m average
  have nonnegative : 0 ≤ (a-b)^2 := OrderedRing.sq_nonneg
  simp only [Semiring.pow_two] at nonnegative
  grind

/-- The three-slope debit is the sum of all three squared differences. -/
theorem three_slope_averaging_energy (a b c m : R)
    (average : 3 * m = a + b + c) :
    3 * (a*a + b*b + c*c - 3*m*m) =
      (a-b)*(a-b) + (a-c)*(a-c) + (b-c)*(b-c) := by
  calc
    _ = 3 * (a*a + b*b + c*c) - (3*m)*(3*m) := by grind
    _ = 3 * (a*a + b*b + c*c) - (a+b+c)*(a+b+c) := by simp only [average]
    _ = _ := by grind

theorem three_slope_averaging_nonincrease (a b c m : R)
    (average : 3 * m = a + b + c) : 3*m*m ≤ a*a + b*b + c*c := by
  have identity := three_slope_averaging_energy a b c m average
  have ab : 0 ≤ (a-b)^2 := OrderedRing.sq_nonneg
  have ac : 0 ≤ (a-c)^2 := OrderedRing.sq_nonneg
  have bc : 0 ≤ (b-c)^2 := OrderedRing.sq_nonneg
  simp only [Semiring.pow_two] at ab ac bc
  grind

/-- A nonflat graft retains its full overlap with the original reference. -/
theorem nonflat_graft_square (d : Nat) (r V P : Nat → R) (mu : R) :
    coordinateDot d (fun c => r c + mu * (V c - P c))
      (fun c => r c + mu * (V c - P c)) =
    coordinateDot d r r +
      2 * mu * (coordinateDot d r V - coordinateDot d r P) +
      mu * mu * (coordinateDot d V V - 2 * coordinateDot d V P +
        coordinateDot d P P) := by
  have cross : coordinateDot d r (fun c => mu * (V c - P c)) =
      mu * (coordinateDot d r V - coordinateDot d r P) := by
    rw [dot_scale_right, dot_sub_right]
  have correction : coordinateDot d (fun c => mu * (V c - P c))
      (fun c => mu * (V c - P c)) =
      mu * mu * (coordinateDot d V V - 2 * coordinateDot d V P +
        coordinateDot d P P) := by
    rw [dot_scale_left, dot_scale_right, dot_sub_square]
    grind
  rw [dot_add_square, cross, correction]
  grind

end CoxeterProofs.OmittedPrefixReferenceKernel

#print axioms CoxeterProofs.OmittedPrefixReferenceKernel.supported_reference_pressure
#print axioms CoxeterProofs.OmittedPrefixReferenceKernel.supported_reference_energy_lower
#print axioms CoxeterProofs.OmittedPrefixReferenceKernel.omitted_reference_contradiction
#print axioms CoxeterProofs.OmittedPrefixReferenceKernel.two_slope_averaging_energy
#print axioms CoxeterProofs.OmittedPrefixReferenceKernel.two_slope_averaging_nonincrease
#print axioms CoxeterProofs.OmittedPrefixReferenceKernel.three_slope_averaging_energy
#print axioms CoxeterProofs.OmittedPrefixReferenceKernel.three_slope_averaging_nonincrease
#print axioms CoxeterProofs.OmittedPrefixReferenceKernel.nonflat_graft_square
