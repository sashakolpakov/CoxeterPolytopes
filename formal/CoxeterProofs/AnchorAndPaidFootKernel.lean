import CoxeterProofs.SummedActualFootKernel

/-!
Source preparation only. Actual bounded compilation and axiom-audit records
are required before calling these declarations kernel checked.

Two logical steps used in the fresh A18 work: a saved-domain match plus
reflected Cartan obstruction forces a unique weak index; and a paid actual
foot gives an energy bound from primitive old floors on the SAME vector.
The pressure proof derives the subtraction from nonnegative remaining
weights. It does not subtract an unknown height from a known lower bound.

Actual cap inventory, correctness/output of the rational enumerations,
arithmetic integrality, geometric feasibility, the concrete N19 reference
and its scalar bounds remain explicit inputs, not proved polytope geometry.
-/

namespace CoxeterProofs.AnchorAndPaidFootKernel

open Lean.Grind
open CoxeterProofs.RightAngledCubicScalarKernel
open CoxeterProofs.RightAngledProjectionHeightKernel
open CoxeterProofs.SummedActualFootKernel

/-- A finite kernel computation of the exact N19 reflection identity. -/
theorem same_k_reflection : ∀ i j : Fin 19,
    i.val * (19-i.val) = j.val * (19-j.val) →
    i = j ∨ i.val + j.val = 19 := by
  decide

/-- Script matches are necessary; they are not asserted to realize caps. -/
theorem weak_index_unique {Parameter : Type} (weak : Fin 19 → Prop)
    (parameter : Fin 19 → Parameter) (anchor : Fin 19)
    (anchor_weak : weak anchor)
    (necessary_matches : ∀ i, weak i → parameter i = parameter anchor ∧
      i.val * (19-i.val) = anchor.val * (19-anchor.val))
    (reflected_impossible : ∀ i j, weak i → weak j →
      parameter i = parameter j → i ≠ j → i.val+j.val = 19 → False) :
    ∀ i, weak i → i = anchor := by
  intro i hi
  have matched := necessary_matches i hi
  have alternatives := same_k_reflection i anchor matched.2
  cases alternatives with
  | inl same => exact same
  | inr reflected =>
    by_cases same : i = anchor
    · exact same
    · exact False.elim (reflected_impossible i anchor hi anchor_weak
        matched.1 same reflected)

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

def paidVector (n : Nat) (b : Nat → R) (omega : Nat → Nat → R)
    (mu : R) (V : Nat → R) : Nat → R :=
  fun c => weightedSum n b omega c + mu * V c

/-- The floors concern the original reference, with only retained weights. -/
theorem remaining_weight_pressure (n d : Nat) (b : Nat → R)
    (omega : Nat → Nat → R) (x w : Nat → R)
    (weight : ∀ i, i < n → 0 ≤ b i)
    (floor : ∀ i, i < n → 0 < b i →
      coordinateDot d w (omega i) ≤ coordinateDot d x (omega i)) :
    coordinateDot d w (weightedSum n b omega) ≤
      coordinateDot d x (weightedSum n b omega) := by
  rw [dot_weighted_right, dot_weighted_right]
  apply sumN_mono
  intro i hi
  by_cases active : 0 < b i
  · exact OrderedRing.mul_le_mul_of_nonneg_left (floor i hi active) (weight i hi)
  · have zero : b i = 0 := by have nonnegative := weight i hi; grind
    simpa only [zero, Semiring.zero_mul] using (Std.le_refl (0 : R))

/-- The payment identity is supplied from actual nonnegative coefficients. -/
theorem paid_reference_pressure (n d : Nat) (b : Nat → R)
    (omega : Nat → Nat → R) (x w P V : Nat → R) (mu a H E g : R)
    (weight : ∀ i, i < n → 0 ≤ b i)
    (floor : ∀ i, i < n → 0 < b i →
      coordinateDot d w (omega i) ≤ coordinateDot d x (omega i))
    (mu_nonnegative : 0 ≤ mu)
    (payment : ∀ c, w c = weightedSum n b omega c + mu*a*P c)
    (energy : coordinateDot d w w = E)
    (height : coordinateDot d w P = H)
    (foot : g ≤ coordinateDot d x V) :
    E + mu*(g-a*H) ≤ coordinateDot d x (paidVector n b omega mu V) := by
  have retained := remaining_weight_pressure n d b omega x w weight floor
  have decomposition : w = fun c => weightedSum n b omega c + mu*a*P c :=
    funext payment
  have reference_identity :
      E = coordinateDot d w (weightedSum n b omega) + mu*a*H := by
    calc
      E = coordinateDot d w w := energy.symm
      _ = coordinateDot d w (fun c => weightedSum n b omega c + mu*a*P c) :=
        congrArg (coordinateDot d w) decomposition
      _ = _ := by rw [dot_add_right, dot_scale_right, height]
  have new_pressure := OrderedRing.mul_le_mul_of_nonneg_left foot mu_nonnegative
  change E + mu*(g-a*H) ≤ coordinateDot d x
    (fun c => weightedSum n b omega c + mu*V c)
  rw [dot_add_right, dot_scale_right]
  grind

/-- One completed coordinate square, including the full graft overlap. -/
theorem paid_reference_energy (n d : Nat) (b : Nat → R)
    (omega : Nat → Nat → R) (x w P V : Nat → R) (mu a H E g S : R)
    (weight : ∀ i, i < n → 0 ≤ b i)
    (floor : ∀ i, i < n → 0 < b i →
      coordinateDot d w (omega i) ≤ coordinateDot d x (omega i))
    (mu_nonnegative : 0 ≤ mu)
    (payment : ∀ c, w c = weightedSum n b omega c + mu*a*P c)
    (energy : coordinateDot d w w = E)
    (height : coordinateDot d w P = H)
    (foot : g ≤ coordinateDot d x V)
    (flat : coordinateDot d w (fun c => V c-a*P c) = 0)
    (tail_square : coordinateDot d (fun c => V c-a*P c)
      (fun c => V c-a*P c) = S) :
    E+2*mu*(g-a*H)-mu*mu*S ≤ coordinateDot d x x := by
  have pressure := paid_reference_pressure n d b omega x w P V mu a H E g
    weight floor mu_nonnegative payment energy height foot
  have same_vector : paidVector n b omega mu V =
      fun c => w c + mu*(V c-a*P c) := by
    funext c
    simp only [paidVector, payment c]
    grind
  have norm_identity : coordinateDot d (paidVector n b omega mu V)
      (paidVector n b omega mu V) = E+mu*mu*S := by
    rw [same_vector, dot_add_square, dot_scale_right, flat,
      dot_scale_left, dot_scale_right, tail_square, energy]
    grind
  have square := coordinate_energy_lower d x (paidVector n b omega mu V)
  grind

theorem paid_reference_contradiction (n d : Nat) (b : Nat → R)
    (omega : Nat → Nat → R) (x w P V : Nat → R) (mu a H E g S : R)
    (weight : ∀ i, i < n → 0 ≤ b i)
    (floor : ∀ i, i < n → 0 < b i →
      coordinateDot d w (omega i) ≤ coordinateDot d x (omega i))
    (mu_nonnegative : 0 ≤ mu)
    (payment : ∀ c, w c = weightedSum n b omega c + mu*a*P c)
    (energy : coordinateDot d w w = E)
    (height : coordinateDot d w P = H)
    (foot : g ≤ coordinateDot d x V)
    (flat : coordinateDot d w (fun c => V c-a*P c) = 0)
    (tail_square : coordinateDot d (fun c => V c-a*P c)
      (fun c => V c-a*P c) = S)
    (large : 1 < E+2*mu*(g-a*H)-mu*mu*S)
    (inside : coordinateDot d x x < 1) : False := by
  have lower := paid_reference_energy n d b omega x w P V mu a H E g S
    weight floor mu_nonnegative payment energy height foot flat tail_square
  grind

end CoxeterProofs.AnchorAndPaidFootKernel

#print axioms CoxeterProofs.AnchorAndPaidFootKernel.same_k_reflection
#print axioms CoxeterProofs.AnchorAndPaidFootKernel.weak_index_unique
#print axioms CoxeterProofs.AnchorAndPaidFootKernel.remaining_weight_pressure
#print axioms CoxeterProofs.AnchorAndPaidFootKernel.paid_reference_pressure
#print axioms CoxeterProofs.AnchorAndPaidFootKernel.paid_reference_energy
#print axioms CoxeterProofs.AnchorAndPaidFootKernel.paid_reference_contradiction
