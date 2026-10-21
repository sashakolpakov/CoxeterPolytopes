import CoxeterProofs.RightAngledProjectionHeightKernel

/-!
Compilation and actual axiom-output verification are recorded separately;
source preparation alone is not verification.

Finite-coordinate dual energy for a sum of actual-foot inequalities.
The weighted vector, its pairings, the shift square and the energy bound
are derived. The shift norm bound uses primitive unit normals,
nonpositive off-diagonal products and nonnegative coefficients.

An application supplies its coordinate realization, a positive dual
reference with w.n_i = beta_i, and each actual-foot comparison. Actual
feasibility, the strong-or-dotted inventory and radical/scalar estimates
remain external. No A20 or dimension conclusion is formalized here.
All index hypotheses are restricted to the finite range actually summed.
-/

namespace CoxeterProofs.SummedActualFootKernel

open Lean.Grind
open CoxeterProofs.RightAngledCubicScalarKernel
open CoxeterProofs.RightAngledProjectionHeightKernel

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

def weightedSum (n : Nat) (a : Nat → R) (v : Nat → Nat → R) : Nat → R :=
  fun c => sumN n (fun i => a i * v i c)

theorem dot_zero_left (d : Nat) (v : Nat → R) :
    coordinateDot d (fun _ => 0) v = 0 := by
  simp only [coordinateDot, Semiring.zero_mul, sumN_const, Semiring.mul_zero]

theorem dot_add_right (d : Nat) (x y z : Nat → R) :
    coordinateDot d x (fun c => y c + z c) =
      coordinateDot d x y + coordinateDot d x z := by
  calc
    _ = coordinateDot d (fun c => y c + z c) x := dot_symm d _ _
    _ = coordinateDot d y x + coordinateDot d z x := dot_add_left d y z x
    _ = _ := by rw [dot_symm d y x, dot_symm d z x]

theorem dot_scale_right (d : Nat) (a : R) (x y : Nat → R) :
    coordinateDot d x (fun c => a * y c) = a * coordinateDot d x y := by
  calc
    _ = coordinateDot d (fun c => a * y c) x := dot_symm d _ _
    _ = a * coordinateDot d y x := dot_scale_left d a y x
    _ = _ := by rw [dot_symm d y x]

theorem dot_sub_left (d : Nat) (x y z : Nat → R) :
    coordinateDot d (fun c => x c - y c) z =
      coordinateDot d x z - coordinateDot d y z := by
  unfold coordinateDot
  have entries := sumN_congr d
    (fun c => (x c - y c) * z c)
    (fun c => x c * z c - y c * z c) (by intro c hc; grind)
  rw [sumN_sub] at entries
  exact entries

theorem dot_sub_right (d : Nat) (x y z : Nat → R) :
    coordinateDot d x (fun c => y c - z c) =
      coordinateDot d x y - coordinateDot d x z := by
  calc
    _ = coordinateDot d (fun c => y c - z c) x := dot_symm d _ _
    _ = coordinateDot d y x - coordinateDot d z x := dot_sub_left d y z x
    _ = _ := by rw [dot_symm d y x, dot_symm d z x]

theorem dot_add_square (d : Nat) (x y : Nat → R) :
    coordinateDot d (fun c => x c + y c) (fun c => x c + y c) =
      coordinateDot d x x + 2 * coordinateDot d x y + coordinateDot d y y := by
  rw [dot_add_left, dot_add_right, dot_add_right, dot_symm d y x]
  grind

theorem dot_sub_square (d : Nat) (x y : Nat → R) :
    coordinateDot d (fun c => x c - y c) (fun c => x c - y c) =
      coordinateDot d x x - 2 * coordinateDot d x y + coordinateDot d y y := by
  rw [dot_sub_left, dot_sub_right, dot_sub_right, dot_symm d y x]
  grind

theorem dot_weighted_left (n d : Nat) (a : Nat → R)
    (v : Nat → Nat → R) (z : Nat → R) :
    coordinateDot d (weightedSum n a v) z =
      sumN n (fun i => a i * coordinateDot d (v i) z) := by
  induction n with
  | zero =>
    change coordinateDot d (fun _ => 0) z = 0
    exact dot_zero_left d z
  | succ n ih =>
    change coordinateDot d
      (fun c => weightedSum n a v c + a n * v n c) z =
      sumN n (fun i => a i * coordinateDot d (v i) z) +
        a n * coordinateDot d (v n) z
    rw [dot_add_left, dot_scale_left, ih]

theorem dot_weighted_right (n d : Nat) (a : Nat → R)
    (v : Nat → Nat → R) (z : Nat → R) :
    coordinateDot d z (weightedSum n a v) =
      sumN n (fun i => a i * coordinateDot d z (v i)) := by
  rw [dot_symm d z (weightedSum n a v), dot_weighted_left]
  apply sumN_congr
  intro i hi
  rw [dot_symm d (v i) z]

/-- The whole weighted Gram expansion is proved from coordinate sums. -/
theorem weighted_norm_identity (n d : Nat) (a : Nat → R)
    (v : Nat → Nat → R) :
    coordinateDot d (weightedSum n a v) (weightedSum n a v) =
      sumN n (fun i => sumN n
        (fun j => a i * a j * coordinateDot d (v i) (v j))) := by
  rw [dot_weighted_left]
  apply sumN_congr
  intro i hi
  rw [dot_weighted_right]
  calc
    a i * sumN n (fun j => a j * coordinateDot d (v i) (v j)) =
        sumN n (fun j => a i * (a j * coordinateDot d (v i) (v j))) :=
      (sumN_scale n (a i) _).symm
    _ = _ := sumN_congr n _ _ (by intro j hj; grind)

theorem sumN_mono (n : Nat) (f g : Nat → R)
    (bound : ∀ i, i < n → f i ≤ g i) : sumN n f ≤ sumN n g := by
  have nonnegative := sumN_nonneg n (fun i => g i - f i)
    (by intro i hi; have row := bound i hi; grind)
  rw [sumN_sub] at nonnegative
  grind

/-- No shift-square estimate is assumed: all cross terms are controlled
by the primitive obtuse products. The family may have any finite size. -/
theorem obtuse_weighted_norm_bound (n d : Nat) (a : Nat → R)
    (normal : Nat → Nat → R)
    (nonnegative : ∀ i, i < n → 0 ≤ a i)
    (unit : ∀ i, i < n → coordinateDot d (normal i) (normal i) = 1)
    (obtuse : ∀ i, i < n → ∀ j, j < n → i ≠ j →
      coordinateDot d (normal i) (normal j) ≤ 0) :
    coordinateDot d (weightedSum n a normal) (weightedSum n a normal) ≤
      sumN n (fun i => a i * a i) := by
  induction n with
  | zero =>
    change coordinateDot d (fun _ => 0) (fun _ => 0) ≤ 0
    simpa only [dot_zero_left] using (Std.le_refl (0 : R))
  | succ n ih =>
    have previous := ih
      (fun i hi => nonnegative i (by omega))
      (fun i hi => unit i (by omega))
      (fun i hi j hj different => obtuse i (by omega) j (by omega) different)
    have cross : coordinateDot d (weightedSum n a normal) (normal n) ≤ 0 := by
      rw [dot_weighted_left]
      apply sumN_nonpos
      intro i hi
      exact OrderedRing.mul_nonpos_of_nonneg_of_nonpos
        (nonnegative i (by omega))
        (obtuse i (by omega) n (by omega) (by omega))
    have scaledCross := OrderedRing.mul_nonpos_of_nonneg_of_nonpos
      (nonnegative n (by omega)) cross
    have lastSquare :
        coordinateDot d (fun c => a n * normal n c)
          (fun c => a n * normal n c) = a n * a n := by
      rw [dot_scale_left, dot_scale_right, unit n (by omega)]
      grind
    change coordinateDot d
      (fun c => weightedSum n a normal c + a n * normal n c)
      (fun c => weightedSum n a normal c + a n * normal n c) ≤
      sumN n (fun i => a i * a i) + a n * a n
    rw [dot_add_square, dot_scale_right, lastSquare]
    grind

/-- E=sum beta_i H_i is a theorem about the actual weighted reference. -/
theorem reference_energy_identity (n d : Nat) (beta : Nat → R)
    (omega : Nat → Nat → R) :
    coordinateDot d (weightedSum n beta omega) (weightedSum n beta omega) =
      sumN n (fun i => beta i *
        coordinateDot d (weightedSum n beta omega) (omega i)) := by
  rw [dot_weighted_left]
  apply sumN_congr
  intro i hi
  rw [dot_symm d (omega i) (weightedSum n beta omega)]

theorem weighted_feet_identity (n : Nat) (beta u : Nat → R)
    (omega normal : Nat → Nat → R) :
    weightedSum n beta (fun i c => omega i c - u i * normal i c) =
      fun c => weightedSum n beta omega c -
        weightedSum n (fun i => beta i * u i) normal c := by
  funext c
  unfold weightedSum
  have entries := sumN_congr n
    (fun i => beta i * (omega i c - u i * normal i c))
    (fun i => beta i * omega i c - (beta i * u i) * normal i c)
    (by intro i hi; grind)
  rw [sumN_sub] at entries
  exact entries

/-- Summing actual-foot pressures uses one vector, including adjacent feet. -/
theorem summed_foot_pressure (n d : Nat) (beta u T : Nat → R)
    (omega normal : Nat → Nat → R) (x : Nat → R)
    (beta_nonnegative : ∀ i, i < n → 0 ≤ beta i)
    (foot : ∀ i, i < n →
      T i ≤ coordinateDot d x (fun c => omega i c - u i * normal i c)) :
    sumN n (fun i => beta i * T i) ≤
      coordinateDot d x (fun c => weightedSum n beta omega c -
        weightedSum n (fun i => beta i * u i) normal c) := by
  rw [← weighted_feet_identity n beta u omega normal, dot_weighted_right]
  apply sumN_mono
  intro i hi
  exact OrderedRing.mul_le_mul_of_nonneg_left (foot i hi) (beta_nonnegative i hi)

theorem reference_shift_pairing (n d : Nat) (beta u : Nat → R)
    (omega normal : Nat → Nat → R)
    (dual : ∀ i, i < n →
      coordinateDot d (weightedSum n beta omega) (normal i) = beta i) :
    coordinateDot d (weightedSum n beta omega)
      (weightedSum n (fun i => beta i * u i) normal) =
      sumN n (fun i => beta i * beta i * u i) := by
  rw [dot_weighted_right]
  apply sumN_congr
  intro i hi
  rw [dual i hi]
  grind

/-- The exact shifted-square expansion is derived, not an input field. -/
theorem shifted_reference_square (n d : Nat) (beta u : Nat → R)
    (omega normal : Nat → Nat → R)
    (dual : ∀ i, i < n →
      coordinateDot d (weightedSum n beta omega) (normal i) = beta i) :
    coordinateDot d
      (fun c => weightedSum n beta omega c -
        weightedSum n (fun i => beta i * u i) normal c)
      (fun c => weightedSum n beta omega c -
        weightedSum n (fun i => beta i * u i) normal c) =
      coordinateDot d (weightedSum n beta omega) (weightedSum n beta omega) -
        2 * sumN n (fun i => beta i * beta i * u i) +
        coordinateDot d (weightedSum n (fun i => beta i * u i) normal)
          (weightedSum n (fun i => beta i * u i) normal) := by
  rw [dot_sub_square, reference_shift_pairing n d beta u omega normal dual]

theorem coordinate_energy_lower (d : Nat) (x W : Nat → R) :
    2 * coordinateDot d x W - coordinateDot d W W ≤ coordinateDot d x x := by
  have nonnegative := coordinateDot_self_nonneg d (fun c => x c - W c)
  rw [dot_sub_square] at nonnegative
  grind

/-- Generic structural certificate. Primitive normal signs give the shift
norm bound; primitive feet give the pressure; completing the coordinate
square gives the final energy. No nonadjacency or u_i<=1 is needed. -/
theorem summed_actual_foot_energy (n d : Nat) (beta u T : Nat → R)
    (omega normal : Nat → Nat → R) (x : Nat → R)
    (beta_nonnegative : ∀ i, i < n → 0 ≤ beta i)
    (u_nonnegative : ∀ i, i < n → 0 ≤ u i)
    (unit : ∀ i, i < n → coordinateDot d (normal i) (normal i) = 1)
    (obtuse : ∀ i, i < n → ∀ j, j < n → i ≠ j →
      coordinateDot d (normal i) (normal j) ≤ 0)
    (dual : ∀ i, i < n →
      coordinateDot d (weightedSum n beta omega) (normal i) = beta i)
    (foot : ∀ i, i < n →
      T i ≤ coordinateDot d x (fun c => omega i c - u i * normal i c)) :
    coordinateDot d (weightedSum n beta omega) (weightedSum n beta omega) +
      sumN n (fun i =>
        2 * beta i * (T i - coordinateDot d (weightedSum n beta omega) (omega i)) +
          (2 * u i - u i * u i) * (beta i * beta i)) ≤ coordinateDot d x x := by
  let w := weightedSum n beta omega
  let Z := weightedSum n (fun i => beta i * u i) normal
  let W := fun c => w c - Z c
  have pressure := summed_foot_pressure n d beta u T omega normal x
    beta_nonnegative foot
  have shiftBound := obtuse_weighted_norm_bound n d (fun i => beta i * u i) normal
    (fun i hi => OrderedRing.mul_nonneg (beta_nonnegative i hi) (u_nonnegative i hi))
    unit obtuse
  have square := shifted_reference_square n d beta u omega normal dual
  have energy := coordinate_energy_lower d x W
  have reference := reference_energy_identity n d beta omega
  have shiftTerms : sumN n (fun i => (beta i * u i) * (beta i * u i)) =
      sumN n (fun i => beta i * beta i * u i * u i) :=
    sumN_congr n _ _ (by intro i hi; grind)
  rw [shiftTerms] at shiftBound
  have improvement :
      sumN n (fun i =>
        2 * beta i * (T i - coordinateDot d w (omega i)) +
          (2 * u i - u i * u i) * (beta i * beta i)) =
      2 * sumN n (fun i => beta i * T i) -
        2 * sumN n (fun i => beta i * coordinateDot d w (omega i)) +
        2 * sumN n (fun i => beta i * beta i * u i) -
        sumN n (fun i => beta i * beta i * u i * u i) := by
    have entries := sumN_congr n
      (fun i => 2 * beta i * (T i - coordinateDot d w (omega i)) +
        (2 * u i - u i * u i) * (beta i * beta i))
      (fun i => 2 * (beta i * T i) -
        2 * (beta i * coordinateDot d w (omega i)) +
        2 * (beta i * beta i * u i) - beta i * beta i * u i * u i)
      (by intro i hi; grind)
    simp only [sumN_sub, sumN_add, sumN_scale] at entries
    simpa only [sumN_add] using entries
  change coordinateDot d w w +
    sumN n (fun i => 2 * beta i * (T i - coordinateDot d w (omega i)) +
      (2 * u i - u i * u i) * (beta i * beta i)) ≤ coordinateDot d x x
  rw [improvement]
  change sumN n (fun i => beta i * T i) ≤ coordinateDot d x W at pressure
  change coordinateDot d Z Z ≤
    sumN n (fun i => beta i * beta i * u i * u i) at shiftBound
  change coordinateDot d W W = coordinateDot d w w -
    2 * sumN n (fun i => beta i * beta i * u i) + coordinateDot d Z Z at square
  change coordinateDot d w w =
    sumN n (fun i => beta i * coordinateDot d w (omega i)) at reference
  grind

end CoxeterProofs.SummedActualFootKernel

#print axioms CoxeterProofs.SummedActualFootKernel.dot_weighted_left
#print axioms CoxeterProofs.SummedActualFootKernel.weighted_norm_identity
#print axioms CoxeterProofs.SummedActualFootKernel.obtuse_weighted_norm_bound
#print axioms CoxeterProofs.SummedActualFootKernel.reference_energy_identity
#print axioms CoxeterProofs.SummedActualFootKernel.weighted_feet_identity
#print axioms CoxeterProofs.SummedActualFootKernel.summed_foot_pressure
#print axioms CoxeterProofs.SummedActualFootKernel.reference_shift_pairing
#print axioms CoxeterProofs.SummedActualFootKernel.shifted_reference_square
#print axioms CoxeterProofs.SummedActualFootKernel.coordinate_energy_lower
#print axioms CoxeterProofs.SummedActualFootKernel.summed_actual_foot_energy
