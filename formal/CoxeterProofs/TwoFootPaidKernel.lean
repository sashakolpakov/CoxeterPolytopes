import CoxeterProofs.MultiAnchorPaidFootKernel

/-!
Source preparation only: bounded compilation and actual axiom outputs are
required before describing these declarations as kernel checked.

Primitive capacities pay TWO actual-foot pressures on the SAME coordinate
vector. The remaining coefficients, payment identity, paid pressure and
full square are derived. Neither residual is assumed orthogonal to the
reference, and the mutual residual product is retained. Actual geometry,
available floors, reference realization, residual estimates and numerical
strict comparisons remain application inputs. No B20 or dimension theorem
is formalized here.
-/

namespace CoxeterProofs.TwoFootPaidKernel

open Lean.Grind
open CoxeterProofs.RightAngledCubicScalarKernel
open CoxeterProofs.RightAngledProjectionHeightKernel
open CoxeterProofs.SummedActualFootKernel
open CoxeterProofs.AnchorAndPaidFootKernel

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

def footResidual (n : Nat) (a : Nat → R) (omega : Nat → Nat → R)
    (Z : Nat → R) : Nat → R :=
  fun c => Z c - weightedSum n a omega c

def twoPaidVector (n : Nat) (beta a b : Nat → R)
    (omega : Nat → Nat → R) (mu nu : R) (Z U : Nat → R) : Nat → R :=
  fun c => weightedSum n beta omega c +
    mu * footResidual n a omega Z c + nu * footResidual n b omega U c

/-- Joint affordability derives the old/residual signs and support inclusion. -/
theorem two_payment_coefficients (n : Nat) (beta a b : Nat → R) (mu nu : R)
    (mu_nonnegative : 0 ≤ mu) (nu_nonnegative : 0 ≤ nu)
    (a_nonnegative : ∀ i, i < n → 0 ≤ a i)
    (b_nonnegative : ∀ i, i < n → 0 ≤ b i)
    (capacity : ∀ i, i < n → mu * a i + nu * b i ≤ beta i) :
    ∀ i, i < n → 0 ≤ beta i ∧ 0 ≤ beta i - mu * a i - nu * b i ∧
      (0 < beta i - mu * a i - nu * b i → 0 < beta i) := by
  intro i hi
  have first := OrderedRing.mul_nonneg mu_nonnegative (a_nonnegative i hi)
  have second := OrderedRing.mul_nonneg nu_nonnegative (b_nonnegative i hi)
  have available := capacity i hi
  constructor
  · grind
  · constructor
    · grind
    · intro active
      grind

/-- Both anchor costs are removed from the actual weighted reference. -/
theorem two_weighted_payment_decomposition (n : Nat) (beta a b : Nat → R)
    (omega : Nat → Nat → R) (mu nu : R) :
    ∀ c, weightedSum n beta omega c =
      weightedSum n (fun i => beta i - mu * a i - nu * b i) omega c +
        mu * weightedSum n a omega c + nu * weightedSum n b omega c := by
  intro c
  unfold weightedSum
  have entries := sumN_congr n
    (fun i => beta i * omega i c)
    (fun i => (beta i - mu * a i - nu * b i) * omega i c +
      mu * (a i * omega i c) + nu * (b i * omega i c))
    (by intro i hi; grind)
  rw [sumN_add, sumN_add, sumN_scale, sumN_scale] at entries
  exact entries

/-- The full square retains BOTH reference products and the mutual product. -/
theorem two_shift_square (d : Nat) (w r s : Nat → R) (mu nu : R) :
    coordinateDot d (fun c => w c + mu * r c + nu * s c)
      (fun c => w c + mu * r c + nu * s c) =
      coordinateDot d w w + 2 * mu * coordinateDot d w r +
        2 * nu * coordinateDot d w s + mu * mu * coordinateDot d r r +
        2 * mu * nu * coordinateDot d r s + nu * nu * coordinateDot d s s := by
  simp only [dot_add_left, dot_add_right, dot_scale_left, dot_scale_right]
  rw [dot_symm d r w, dot_symm d s w, dot_symm d s r]
  grind

/-- Each original floor is used only with its nonnegative remaining weight. -/
theorem two_foot_paid_pressure (n d : Nat) (beta a b : Nat → R)
    (omega : Nat → Nat → R) (x Z U : Nat → R) (mu nu h t : R)
    (mu_nonnegative : 0 ≤ mu) (nu_nonnegative : 0 ≤ nu)
    (a_nonnegative : ∀ i, i < n → 0 ≤ a i)
    (b_nonnegative : ∀ i, i < n → 0 ≤ b i)
    (capacity : ∀ i, i < n → mu * a i + nu * b i ≤ beta i)
    (floors : ∀ i, i < n → 0 < beta i →
      coordinateDot d (weightedSum n beta omega) (omega i) ≤
        coordinateDot d x (omega i))
    (first_foot : h ≤ coordinateDot d x Z)
    (second_foot : t ≤ coordinateDot d x U) :
    coordinateDot d (weightedSum n beta omega) (weightedSum n beta omega) +
      mu * (h - coordinateDot d (weightedSum n beta omega) (weightedSum n a omega)) +
      nu * (t - coordinateDot d (weightedSum n beta omega) (weightedSum n b omega)) ≤
        coordinateDot d x (twoPaidVector n beta a b omega mu nu Z U) := by
  let w := weightedSum n beta omega
  let remaining := weightedSum n (fun i => beta i - mu * a i - nu * b i) omega
  have coefficients := two_payment_coefficients n beta a b mu nu
    mu_nonnegative nu_nonnegative a_nonnegative b_nonnegative capacity
  have retained := remaining_weight_pressure n d
    (fun i => beta i - mu * a i - nu * b i) omega x w
    (fun i hi => (coefficients i hi).2.1)
    (fun i hi active => floors i hi ((coefficients i hi).2.2 active))
  have decomposition : w = fun c => remaining c +
      mu * weightedSum n a omega c + nu * weightedSum n b omega c := by
    funext c
    exact two_weighted_payment_decomposition n beta a b omega mu nu c
  have reference_identity : coordinateDot d w w =
      coordinateDot d w remaining + mu * coordinateDot d w (weightedSum n a omega) +
        nu * coordinateDot d w (weightedSum n b omega) := by
    calc
      coordinateDot d w w = coordinateDot d w (fun c => remaining c +
          mu * weightedSum n a omega c + nu * weightedSum n b omega c) :=
        congrArg (coordinateDot d w) decomposition
      _ = _ := by rw [dot_add_right, dot_add_right, dot_scale_right, dot_scale_right]
  have paid_identity : twoPaidVector n beta a b omega mu nu Z U =
      fun c => remaining c + mu * Z c + nu * U c := by
    funext c
    have payment := two_weighted_payment_decomposition n beta a b omega mu nu c
    change weightedSum n beta omega c +
      mu * (Z c - weightedSum n a omega c) +
      nu * (U c - weightedSum n b omega c) =
        remaining c + mu * Z c + nu * U c
    change weightedSum n beta omega c = remaining c +
      mu * weightedSum n a omega c + nu * weightedSum n b omega c at payment
    grind
  have first_pressure := OrderedRing.mul_le_mul_of_nonneg_left first_foot mu_nonnegative
  have second_pressure := OrderedRing.mul_le_mul_of_nonneg_left second_foot nu_nonnegative
  change coordinateDot d w w +
    mu * (h - coordinateDot d w (weightedSum n a omega)) +
    nu * (t - coordinateDot d w (weightedSum n b omega)) ≤
      coordinateDot d x (twoPaidVector n beta a b omega mu nu Z U)
  rw [paid_identity, dot_add_right, dot_add_right, dot_scale_right, dot_scale_right]
  change coordinateDot d w remaining ≤ coordinateDot d x remaining at retained
  grind

/-- Completing ONE coordinate square gives the full nonorthogonal energy gain. -/
theorem two_foot_paid_energy (n d : Nat) (beta a b : Nat → R)
    (omega : Nat → Nat → R) (x Z U : Nat → R) (mu nu h t : R)
    (mu_nonnegative : 0 ≤ mu) (nu_nonnegative : 0 ≤ nu)
    (a_nonnegative : ∀ i, i < n → 0 ≤ a i)
    (b_nonnegative : ∀ i, i < n → 0 ≤ b i)
    (capacity : ∀ i, i < n → mu * a i + nu * b i ≤ beta i)
    (floors : ∀ i, i < n → 0 < beta i →
      coordinateDot d (weightedSum n beta omega) (omega i) ≤
        coordinateDot d x (omega i))
    (first_foot : h ≤ coordinateDot d x Z)
    (second_foot : t ≤ coordinateDot d x U) :
    coordinateDot d (weightedSum n beta omega) (weightedSum n beta omega) +
      2 * mu * (h - coordinateDot d (weightedSum n beta omega) Z) +
      2 * nu * (t - coordinateDot d (weightedSum n beta omega) U) -
      mu * mu * coordinateDot d (footResidual n a omega Z) (footResidual n a omega Z) -
      2 * mu * nu * coordinateDot d (footResidual n a omega Z) (footResidual n b omega U) -
      nu * nu * coordinateDot d (footResidual n b omega U) (footResidual n b omega U) ≤
        coordinateDot d x x := by
  let w := weightedSum n beta omega
  let r := footResidual n a omega Z
  let s := footResidual n b omega U
  have pressure := two_foot_paid_pressure n d beta a b omega x Z U mu nu h t
    mu_nonnegative nu_nonnegative a_nonnegative b_nonnegative capacity floors
    first_foot second_foot
  have square := two_shift_square d w r s mu nu
  have completed := coordinate_energy_lower d x (twoPaidVector n beta a b omega mu nu Z U)
  have first_pairing : coordinateDot d w r =
      coordinateDot d w Z - coordinateDot d w (weightedSum n a omega) := by
    exact dot_sub_right d w Z (weightedSum n a omega)
  have second_pairing : coordinateDot d w s =
      coordinateDot d w U - coordinateDot d w (weightedSum n b omega) := by
    exact dot_sub_right d w U (weightedSum n b omega)
  change coordinateDot d (twoPaidVector n beta a b omega mu nu Z U)
    (twoPaidVector n beta a b omega mu nu Z U) =
      coordinateDot d w w + 2 * mu * coordinateDot d w r +
        2 * nu * coordinateDot d w s + mu * mu * coordinateDot d r r +
        2 * mu * nu * coordinateDot d r s + nu * nu * coordinateDot d s s at square
  change coordinateDot d w w +
    mu * (h - coordinateDot d w (weightedSum n a omega)) +
    nu * (t - coordinateDot d w (weightedSum n b omega)) ≤
      coordinateDot d x (twoPaidVector n beta a b omega mu nu Z U) at pressure
  change coordinateDot d w w + 2 * mu * (h - coordinateDot d w Z) +
    2 * nu * (t - coordinateDot d w U) - mu * mu * coordinateDot d r r -
    2 * mu * nu * coordinateDot d r s - nu * nu * coordinateDot d s s ≤
      coordinateDot d x x
  have twice := OrderedRing.mul_le_mul_of_nonneg_left pressure
    (show (0 : R) ≤ 2 by grind)
  calc
    _ = 2 * (coordinateDot d w w +
          mu * (h - coordinateDot d w (weightedSum n a omega)) +
          nu * (t - coordinateDot d w (weightedSum n b omega))) -
        coordinateDot d (twoPaidVector n beta a b omega mu nu Z U)
          (twoPaidVector n beta a b omega mu nu Z U) := by
      rw [square, first_pairing, second_pairing]
      grind
    _ ≤ 2 * coordinateDot d x (twoPaidVector n beta a b omega mu nu Z U) -
        coordinateDot d (twoPaidVector n beta a b omega mu nu Z U)
          (twoPaidVector n beta a b omega mu nu Z U) := by grind
    _ ≤ coordinateDot d x x := completed

end CoxeterProofs.TwoFootPaidKernel

#print axioms CoxeterProofs.TwoFootPaidKernel.two_payment_coefficients
#print axioms CoxeterProofs.TwoFootPaidKernel.two_weighted_payment_decomposition
#print axioms CoxeterProofs.TwoFootPaidKernel.two_shift_square
#print axioms CoxeterProofs.TwoFootPaidKernel.two_foot_paid_pressure
#print axioms CoxeterProofs.TwoFootPaidKernel.two_foot_paid_energy
