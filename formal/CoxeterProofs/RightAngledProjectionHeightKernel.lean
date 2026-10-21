import CoxeterProofs.RightAngledCubicScalarKernel

/-!
Coordinate kernel for the maximum-height projection argument.
Finite Lorentz coordinates and elementary wall signs are the inputs. Unit
norm, orthogonality, strict wall feasibility and squared-height descent are
derived below; no projection or dimension conclusion is an input field.

A real application still supplies positive normalizers, actual unit walls,
finite maximum choices, the future-point/relative-interior interpretation,
orthogonal ridge-neighbor completeness and inherited face geometry.
No concrete real model, incidence theorem or dimension wrapper is supplied.
-/

namespace CoxeterProofs.RightAngledProjectionHeightKernel

open Lean.Grind
open CoxeterProofs.RightAngledCubicScalarKernel

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

abbrev Vector (R : Type u) := R × (Nat → R)

def addV (x y : Vector R) : Vector R :=
  (x.1 + y.1, fun i => x.2 i + y.2 i)

def scaleV (a : R) (x : Vector R) : Vector R :=
  (a * x.1, fun i => a * x.2 i)

def lorentz (d : Nat) (x y : Vector R) : R :=
  -(x.1 * y.1) + coordinateDot d x.2 y.2

theorem dot_symm (d : Nat) (x y : Nat → R) :
    coordinateDot d x y = coordinateDot d y x := by
  apply sumN_congr
  intro i hi
  grind

theorem dot_add_left (d : Nat) (x y z : Nat → R) :
    coordinateDot d (fun i => x i + y i) z =
      coordinateDot d x z + coordinateDot d y z := by
  unfold coordinateDot
  have entries := sumN_congr d
    (fun i => (x i + y i) * z i)
    (fun i => x i * z i + y i * z i)
    (by intro i hi; grind)
  rw [sumN_add] at entries
  exact entries

theorem dot_scale_left (d : Nat) (a : R) (x y : Nat → R) :
    coordinateDot d (fun i => a * x i) y = a * coordinateDot d x y := by
  unfold coordinateDot
  have entries := sumN_congr d
    (fun i => (a * x i) * y i) (fun i => a * (x i * y i))
    (by intro i hi; grind)
  rw [sumN_scale] at entries
  exact entries

theorem lorentz_symm (d : Nat) (x y : Vector R) :
    lorentz d x y = lorentz d y x := by
  have spatial := dot_symm d x.2 y.2
  unfold lorentz
  grind

theorem lorentz_add_left (d : Nat) (x y z : Vector R) :
    lorentz d (addV x y) z = lorentz d x z + lorentz d y z := by
  have spatial := dot_add_left d x.2 y.2 z.2
  simp only [lorentz, addV]
  grind

theorem lorentz_scale_left (d : Nat) (a : R) (x y : Vector R) :
    lorentz d (scaleV a x) y = a * lorentz d x y := by
  have spatial := dot_scale_left d a x.2 y.2
  simp only [lorentz, scaleV]
  grind

theorem lorentz_add_right (d : Nat) (x y z : Vector R) :
    lorentz d x (addV y z) = lorentz d x y + lorentz d x z := by
  calc
    _ = lorentz d (addV y z) x := lorentz_symm d x (addV y z)
    _ = lorentz d y x + lorentz d z x := lorentz_add_left d y z x
    _ = _ := by rw [lorentz_symm d y x, lorentz_symm d z x]

theorem lorentz_scale_right (d : Nat) (a : R) (x y : Vector R) :
    lorentz d x (scaleV a y) = a * lorentz d x y := by
  calc
    _ = lorentz d (scaleV a y) x := lorentz_symm d x (scaleV a y)
    _ = a * lorentz d y x := lorentz_scale_left d a y x
    _ = _ := by rw [lorentz_symm d y x]

def rawProjection (q n : Vector R) (t : R) : Vector R :=
  addV q (scaleV t n)

def normalizedProjection (q n : Vector R) (t lambda : R) : Vector R :=
  scaleV lambda (rawProjection q n t)

def height (d : Nat) (q m : Vector R) : R := -lorentz d q m

/-- All three raw projection identities follow from the coordinate form. -/
theorem raw_projection_identities (d : Nat) (q n : Vector R) (t : R)
    (q_unit : lorentz d q q = -1) (n_unit : lorentz d n n = 1)
    (incidence : lorentz d q n = -t) :
    lorentz d (rawProjection q n t) n = 0 ∧
    lorentz d (rawProjection q n t) (rawProjection q n t) = -1 - t*t ∧
    lorentz d q (rawProjection q n t) = -1 - t*t := by
  have reverse : lorentz d n q = -t := by
    rw [lorentz_symm d n q]
    exact incidence
  simp only [rawProjection, lorentz_add_left, lorentz_add_right,
    lorentz_scale_left, lorentz_scale_right]
  constructor
  · grind
  · constructor <;> grind

/-- A strict old halfspace remains strict; no maximality is used here. -/
theorem raw_wall_strict (d : Nat) (q n m : Vector R) (t : R)
    (positive : 0 < t) (inside : lorentz d q m < 0)
    (acute : lorentz d n m ≤ 0) :
    lorentz d (rawProjection q n t) m < 0 := by
  have correction := OrderedRing.mul_nonpos_of_nonneg_of_nonpos
    (Preorder.le_of_lt positive) acute
  simp only [rawProjection, lorentz_add_left, lorentz_scale_left]
  grind

/-- Negative pairing with the old unit point is proved, rather than a
future-point predicate being assumed. The real interpretation stays external. -/
theorem normalized_projection (d : Nat) (q n : Vector R) (t lambda : R)
    (q_unit : lorentz d q q = -1) (n_unit : lorentz d n n = 1)
    (incidence : lorentz d q n = -t) (positive : 0 < lambda)
    (normalizing : lambda*lambda*(1+t*t) = 1) :
    lorentz d (normalizedProjection q n t lambda)
      (normalizedProjection q n t lambda) = -1 ∧
    lorentz d (normalizedProjection q n t lambda) n = 0 ∧
    lorentz d q (normalizedProjection q n t lambda) < 0 := by
  obtain ⟨orthogonal, square, pairing⟩ :=
    raw_projection_identities d q n t q_unit n_unit incidence
  have nonnegative : 0 ≤ t*t := by
    have squareNonnegative : 0 ≤ t^2 := OrderedRing.sq_nonneg
    simpa only [Semiring.pow_two] using squareNonnegative
  have onePositive : (0 : R) < 1 := OrderedRing.zero_lt_one
  have pastPairing : lorentz d q (rawProjection q n t) < 0 := by grind
  have scaledPairing := OrderedRing.mul_neg_of_pos_of_neg positive pastPairing
  simp only [normalizedProjection, lorentz_scale_left, lorentz_scale_right]
  constructor
  · grind
  · constructor <;> grind

theorem normalized_wall_strict (d : Nat) (q n m : Vector R) (t lambda : R)
    (t_positive : 0 < t) (lambda_positive : 0 < lambda)
    (inside : lorentz d q m < 0) (acute : lorentz d n m ≤ 0) :
    lorentz d (normalizedProjection q n t lambda) m < 0 := by
  have raw := raw_wall_strict d q n m t t_positive inside acute
  have scaled := OrderedRing.mul_neg_of_pos_of_neg lambda_positive raw
  simpa only [normalizedProjection, lorentz_scale_left] using scaled

/-- Restriction to an orthogonal wall gives an exact cleared height identity. -/
theorem orthogonal_height_identity (d : Nat) (q n m : Vector R) (t lambda : R)
    (orthogonal : lorentz d n m = 0)
    (normalizing : lambda*lambda*(1+t*t) = 1) :
    height d (normalizedProjection q n t lambda) m = lambda * height d q m ∧
    (1+t*t) * height d (normalizedProjection q n t lambda) m *
      height d (normalizedProjection q n t lambda) m =
      height d q m * height d q m := by
  have identity : height d (normalizedProjection q n t lambda) m =
      lambda * height d q m := by
    simp only [height, normalizedProjection, rawProjection,
      lorentz_scale_left, lorentz_add_left, orthogonal]
    grind
  refine ⟨identity, ?_⟩
  rw [identity]
  grind

/-- A maximum old height bounds every new orthogonal-neighbor height.
Actual ridge-neighbor completeness is not supplied by this algebra theorem. -/
theorem maximal_height_square_step (d : Nat) (q n m : Vector R) (t lambda : R)
    (orthogonal : lorentz d n m = 0) (lambda_positive : 0 < lambda)
    (normalizing : lambda*lambda*(1+t*t) = 1)
    (old_positive : 0 < height d q m) (maximal : height d q m ≤ t) :
    0 < height d (normalizedProjection q n t lambda) m ∧
    (1+t*t) * height d (normalizedProjection q n t lambda) m *
      height d (normalizedProjection q n t lambda) m ≤ t*t := by
  obtain ⟨identity, squared⟩ :=
    orthogonal_height_identity d q n m t lambda orthogonal normalizing
  have positive := OrderedRing.mul_pos lambda_positive old_positive
  have lower : 0 ≤ t - height d q m := by grind
  have upper : 0 ≤ t + height d q m := by grind
  have difference := OrderedRing.mul_nonneg lower upper
  constructor
  · rw [identity]
    exact positive
  · grind

/-- The strict descent invariant uses only a local cleared recurrence.
The polynomial identity in the proof avoids inversion and square roots. -/
theorem descent_step_strict (k : Nat) (a b : R) (nonnegative : 0 ≤ a)
    (step : (1+a)*b ≤ a) (previous : (k : R)*a < 1) :
    ((k : R)+1)*b < 1 := by
  have onePositive : (0 : R) < 1 := OrderedRing.zero_lt_one
  have kNonnegative : 0 ≤ (k : R) := OrderedRing.natCast_nonneg
  have coefficient : 0 ≤ (k : R)+1 := by grind
  have residual : 0 ≤ a - (1+a)*b := by grind
  have product := OrderedRing.mul_nonneg coefficient residual
  have identity : (1+a)*(1-((k : R)+1)*b) =
      (1-(k : R)*a) + ((k : R)+1)*(a-(1+a)*b) := by grind
  have strictProduct : 0 < (1+a)*(1-((k : R)+1)*b) := by grind
  by_cases result : ((k : R)+1)*b < 1
  · exact result
  · have factor : 0 ≤ 1+a := by grind
    have gap : 1-((k : R)+1)*b ≤ 0 := by grind
    have forbidden := OrderedRing.mul_nonpos_of_nonneg_of_nonpos factor gap
    grind

/-- Every finite chain is covered; no recurrence beyond `steps` is required.
The j=0 assertion is the harmless strict base 0<1. -/
theorem descent_strict (steps : Nat) (a : Nat → R)
    (nonnegative : ∀ j, j ≤ steps → 0 ≤ a j)
    (step : ∀ j, j < steps → (1+a j)*a (j+1) ≤ a j) :
    ∀ j, j ≤ steps → (j : R)*a j < 1 := by
  intro j
  induction j with
  | zero =>
    intro bound
    simpa only [Semiring.natCast_zero, Semiring.zero_mul] using
      (OrderedRing.zero_lt_one (R := R))
  | succ j ih =>
    intro bound
    have earlier : j ≤ steps := by omega
    have advance := descent_step_strict j (a j) (a (j+1))
      (nonnegative j earlier) (step j (by omega)) (ih earlier)
    simpa only [Semiring.natCast_succ] using advance

/-- The height-to-radial strict gap is scalar algebra. Actual radial vectors
and their pair-product inequality still require the coordinate interpretation. -/
theorem radial_gap (p hi hj ri rj z : R) (p_positive : 0 < p)
    (hi_nonnegative : 0 ≤ hi) (hj_nonnegative : 0 ≤ hj)
    (ri_positive : 0 < ri) (rj_positive : 0 < rj)
    (ri_unit : ri*ri = 1+hi*hi) (rj_unit : rj*rj = 1+hj*hj)
    (hi_small : (p-1)*hi*hi < 1) (hj_small : (p-1)*hj*hj < 1)
    (allwall : ri*rj*z ≤ hi*hj) : p*z < 1 := by
  have onePositive : (0 : R) < 1 := OrderedRing.zero_lt_one
  have boundi : p*hi*hi < ri*ri := by grind
  have boundj : p*hj*hj < rj*rj := by grind
  have pi := OrderedRing.mul_nonneg (Preorder.le_of_lt p_positive) hi_nonnegative
  have ii := OrderedRing.mul_nonneg pi hi_nonnegative
  have rjSquare := OrderedRing.mul_pos rj_positive rj_positive
  have weak := OrderedRing.mul_le_mul_of_nonneg_left
    (Preorder.le_of_lt boundj) ii
  have strict := OrderedRing.mul_lt_mul_of_pos_right boundi rjSquare
  have squares : (p*hi*hj)*(p*hi*hj) < (ri*rj)*(ri*rj) := by grind
  have xNonnegative := OrderedRing.mul_nonneg pi hj_nonnegative
  have radiiPositive := OrderedRing.mul_pos ri_positive rj_positive
  have comparison : p*hi*hj < ri*rj := by
    by_cases result : p*hi*hj < ri*rj
    · exact result
    · have difference : 0 ≤ p*hi*hj - ri*rj := by grind
      have sum : 0 ≤ p*hi*hj + ri*rj := by grind
      have forbidden := OrderedRing.mul_nonneg difference sum
      grind
  have scaled := OrderedRing.mul_le_mul_of_nonneg_left allwall
    (Preorder.le_of_lt p_positive)
  have upper : ri*rj*(p*z) < ri*rj := by grind
  by_cases result : p*z < 1
  · exact result
  · have lower : 1 ≤ p*z := by grind
    have forbidden := OrderedRing.mul_le_mul_of_nonneg_left lower
      (Preorder.le_of_lt radiiPositive)
    grind

end CoxeterProofs.RightAngledProjectionHeightKernel

#print axioms CoxeterProofs.RightAngledProjectionHeightKernel.raw_projection_identities
#print axioms CoxeterProofs.RightAngledProjectionHeightKernel.normalized_projection
#print axioms CoxeterProofs.RightAngledProjectionHeightKernel.normalized_wall_strict
#print axioms CoxeterProofs.RightAngledProjectionHeightKernel.orthogonal_height_identity
#print axioms CoxeterProofs.RightAngledProjectionHeightKernel.maximal_height_square_step
#print axioms CoxeterProofs.RightAngledProjectionHeightKernel.descent_step_strict
#print axioms CoxeterProofs.RightAngledProjectionHeightKernel.descent_strict
#print axioms CoxeterProofs.RightAngledProjectionHeightKernel.radial_gap
