import CoxeterProofs.AnchorAndPaidFootKernel

/-!
Source preparation only: bounded compilation and actual axiom outputs are
required before describing these declarations as kernel checked.

Primitive coefficient capacities pay any finite number of old anchors for
ONE new foot pressure on the SAME coordinate vector. The residual weights,
payment decomposition and joint old height are derived, not input fields.
Actual geometry, available old floors, foot feasibility, reference flatness,
residual square and concrete strict scalar comparisons remain application
inputs. No cap inventory, polytope or dimension theorem is formalized here.
-/

namespace CoxeterProofs.MultiAnchorPaidFootKernel

open Lean.Grind
open CoxeterProofs.RightAngledCubicScalarKernel
open CoxeterProofs.RightAngledProjectionHeightKernel
open CoxeterProofs.SummedActualFootKernel
open CoxeterProofs.AnchorAndPaidFootKernel

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

def jointAnchorHeight (n d : Nat) (beta a : Nat → R)
    (omega : Nat → Nat → R) : R :=
  sumN n (fun i => a i * coordinateDot d (weightedSum n beta omega) (omega i))

/-- Affordability also derives nonnegativity of the old coefficients. -/
theorem payment_coefficients (n : Nat) (beta a : Nat → R) (mu : R)
    (mu_nonnegative : 0 ≤ mu)
    (cost_nonnegative : ∀ i, i < n → 0 ≤ a i)
    (capacity : ∀ i, i < n → mu * a i ≤ beta i) :
    ∀ i, i < n → 0 ≤ beta i ∧ 0 ≤ beta i - mu * a i ∧
      (0 < beta i - mu * a i → 0 < beta i) := by
  intro i hi
  have product := OrderedRing.mul_nonneg mu_nonnegative (cost_nonnegative i hi)
  have available := capacity i hi
  constructor
  · grind
  · constructor
    · grind
    · intro active
      grind

/-- Payment is an identity of the actual weighted reference, at every coordinate. -/
theorem weighted_payment_decomposition (n : Nat) (beta a : Nat → R)
    (omega : Nat → Nat → R) (mu : R) :
    ∀ c, weightedSum n beta omega c =
      weightedSum n (fun i => beta i - mu * a i) omega c +
        mu * weightedSum n a omega c := by
  intro c
  unfold weightedSum
  have entries := sumN_congr n
    (fun i => beta i * omega i c)
    (fun i => (beta i - mu * a i) * omega i c + mu * (a i * omega i c))
    (by intro i hi; grind)
  rw [sumN_add, sumN_scale] at entries
  exact entries

/-- The joint height is computed from primitive anchor pairings. -/
theorem joint_anchor_height (n d : Nat) (beta a : Nat → R)
    (omega : Nat → Nat → R) :
    coordinateDot d (weightedSum n beta omega) (weightedSum n a omega) =
      jointAnchorHeight n d beta a omega := by
  unfold jointAnchorHeight
  exact dot_weighted_right n d a omega (weightedSum n beta omega)

/-- Original floors are needed only where the original reference has support. -/
theorem multi_anchor_paid_pressure (n d : Nat) (beta a : Nat → R)
    (omega : Nat → Nat → R) (x V : Nat → R) (mu g : R)
    (mu_nonnegative : 0 ≤ mu)
    (cost_nonnegative : ∀ i, i < n → 0 ≤ a i)
    (capacity : ∀ i, i < n → mu * a i ≤ beta i)
    (floors : ∀ i, i < n → 0 < beta i →
      coordinateDot d (weightedSum n beta omega) (omega i) ≤
        coordinateDot d x (omega i))
    (foot : g ≤ coordinateDot d x V) :
    coordinateDot d (weightedSum n beta omega) (weightedSum n beta omega) +
      mu * (g - jointAnchorHeight n d beta a omega) ≤
        coordinateDot d x
          (paidVector n (fun i => beta i - mu * a i) omega mu V) := by
  have coefficients := payment_coefficients n beta a mu
    mu_nonnegative cost_nonnegative capacity
  have payment : ∀ c, weightedSum n beta omega c =
      weightedSum n (fun i => beta i - mu * a i) omega c +
        mu * 1 * weightedSum n a omega c := by
    intro c
    simpa only [Semiring.mul_one] using
      weighted_payment_decomposition n beta a omega mu c
  have pressure := paid_reference_pressure n d (fun i => beta i - mu * a i)
    omega x (weightedSum n beta omega) (weightedSum n a omega) V mu 1
    (jointAnchorHeight n d beta a omega)
    (coordinateDot d (weightedSum n beta omega) (weightedSum n beta omega)) g
    (fun i hi => (coefficients i hi).2.1)
    (fun i hi active => floors i hi ((coefficients i hi).2.2 active))
    mu_nonnegative payment rfl (joint_anchor_height n d beta a omega) foot
  simpa only [Semiring.one_mul] using pressure

/-- One paid reference gives one norm gain, irrespective of the number of anchors. -/
theorem multi_anchor_paid_energy (n d : Nat) (beta a : Nat → R)
    (omega : Nat → Nat → R) (x V : Nat → R) (mu g S : R)
    (mu_nonnegative : 0 ≤ mu)
    (cost_nonnegative : ∀ i, i < n → 0 ≤ a i)
    (capacity : ∀ i, i < n → mu * a i ≤ beta i)
    (floors : ∀ i, i < n → 0 < beta i →
      coordinateDot d (weightedSum n beta omega) (omega i) ≤
        coordinateDot d x (omega i))
    (foot : g ≤ coordinateDot d x V)
    (flat : coordinateDot d (weightedSum n beta omega)
      (fun c => V c - weightedSum n a omega c) = 0)
    (tail_square : coordinateDot d (fun c => V c - weightedSum n a omega c)
      (fun c => V c - weightedSum n a omega c) = S) :
    coordinateDot d (weightedSum n beta omega) (weightedSum n beta omega) +
      2 * mu * (g - jointAnchorHeight n d beta a omega) - mu * mu * S ≤
        coordinateDot d x x := by
  have coefficients := payment_coefficients n beta a mu
    mu_nonnegative cost_nonnegative capacity
  have payment : ∀ c, weightedSum n beta omega c =
      weightedSum n (fun i => beta i - mu * a i) omega c +
        mu * 1 * weightedSum n a omega c := by
    intro c
    simpa only [Semiring.mul_one] using
      weighted_payment_decomposition n beta a omega mu c
  have lower := paid_reference_energy n d (fun i => beta i - mu * a i)
    omega x (weightedSum n beta omega) (weightedSum n a omega) V mu 1
    (jointAnchorHeight n d beta a omega)
    (coordinateDot d (weightedSum n beta omega) (weightedSum n beta omega)) g S
    (fun i hi => (coefficients i hi).2.1)
    (fun i hi active => floors i hi ((coefficients i hi).2.2 active))
    mu_nonnegative payment rfl (joint_anchor_height n d beta a omega) foot
    (by simpa only [Semiring.one_mul] using flat)
    (by simpa only [Semiring.one_mul] using tail_square)
  simpa only [Semiring.one_mul] using lower

theorem multi_anchor_paid_contradiction (n d : Nat) (beta a : Nat → R)
    (omega : Nat → Nat → R) (x V : Nat → R) (mu g S : R)
    (mu_nonnegative : 0 ≤ mu)
    (cost_nonnegative : ∀ i, i < n → 0 ≤ a i)
    (capacity : ∀ i, i < n → mu * a i ≤ beta i)
    (floors : ∀ i, i < n → 0 < beta i →
      coordinateDot d (weightedSum n beta omega) (omega i) ≤
        coordinateDot d x (omega i))
    (foot : g ≤ coordinateDot d x V)
    (flat : coordinateDot d (weightedSum n beta omega)
      (fun c => V c - weightedSum n a omega c) = 0)
    (tail_square : coordinateDot d (fun c => V c - weightedSum n a omega c)
      (fun c => V c - weightedSum n a omega c) = S)
    (large : 1 <
      coordinateDot d (weightedSum n beta omega) (weightedSum n beta omega) +
        2 * mu * (g - jointAnchorHeight n d beta a omega) - mu * mu * S)
    (inside : coordinateDot d x x < 1) : False := by
  have lower := multi_anchor_paid_energy n d beta a omega x V mu g S
    mu_nonnegative cost_nonnegative capacity floors foot flat tail_square
  grind

end CoxeterProofs.MultiAnchorPaidFootKernel

#print axioms CoxeterProofs.MultiAnchorPaidFootKernel.payment_coefficients
#print axioms CoxeterProofs.MultiAnchorPaidFootKernel.weighted_payment_decomposition
#print axioms CoxeterProofs.MultiAnchorPaidFootKernel.joint_anchor_height
#print axioms CoxeterProofs.MultiAnchorPaidFootKernel.multi_anchor_paid_pressure
#print axioms CoxeterProofs.MultiAnchorPaidFootKernel.multi_anchor_paid_energy
#print axioms CoxeterProofs.MultiAnchorPaidFootKernel.multi_anchor_paid_contradiction
