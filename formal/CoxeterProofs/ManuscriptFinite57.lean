import CoxeterProofs.GuardedRank58
import CoxeterProofs.ManuscriptDataBounds

/-! Conditional finite57 composition from separate exact finite/ideal rows,
actual-domain domination, the internal shape correction and the older58
tail. The final dimension bound is not a hypothesis. -/
namespace CoxeterProofs.ManuscriptFinite57

open ManuscriptDataBounds

structure Certificate where
  scale : Nat
  finiteOffers : List Nat
  terminalIdeals : List Charge

/-- Exact script outcomes after one common positive denominator is cleared. -/
structure ScriptChecks (c : Certificate) : Prop where
  scalePositive : 0 < c.scale
  finiteUpper : ∀ a ∈ c.finiteOffers, a ≤ c.scale * 377
  terminalSurplus :
    ∀ a ∈ c.terminalIdeals, a.cost + c.scale * 36 ≤ a.curvature

structure ActualCharges where
  finiteCharges : List Nat
  idealCharges : List Charge
  deficit : Nat

/-- Schema for the published exact index identity at the actual chamber.
The value must be instantiated by its Khovanskii index expression;
this file does not define the geometric index or prove that identity. -/
structure PublishedIndex (a : ActualCharges) where
  value : Nat
  exactDeficit : a.deficit = value

/-- Internal written geometry. Surplus domination expresses concentration
and actual factor-fiber coverage, not equality of actual cusps with
terminal table rows. The4205 correction includes the two actual extrema. -/
structure WrittenGeometry (c : Certificate) (a : ActualCharges)
    (index : PublishedIndex a) : Prop where
  finiteDomination :
    ∀ x ∈ a.finiteCharges, ∃ y ∈ c.finiteOffers, x ≤ y
  idealSurplusDomination :
    ∀ x ∈ a.idealCharges, ∃ y ∈ c.terminalIdeals,
      x.cost + y.curvature ≤ x.curvature + y.cost
  shapeLower :
    (c.scale * 377) * a.finiteCharges.length +
      curvatures a.idealCharges + c.scale * 4205 ≤ index.value
  paymentUpper : a.deficit ≤ a.finiteCharges.sum + costs a.idealCharges

theorem no_rank58_charges (c : Certificate) (script : ScriptChecks c)
    (a : ActualCharges) (index : PublishedIndex a)
    (geometry : WrittenGeometry c a index) : False := by
  have finite : ∀ x ∈ a.finiteCharges, x ≤ c.scale * 377 := by
    intro x hx
    obtain ⟨y, hy, hxy⟩ := geometry.finiteDomination x hx
    exact Nat.le_trans hxy (script.finiteUpper y hy)
  have finiteSum := sum_le_cap_mul_length a.finiteCharges (c.scale * 377) finite
  have ideal : ∀ x ∈ a.idealCharges, x.cost + c.scale * 36 ≤ x.curvature := by
    intro x hx
    obtain ⟨y, hy, domination⟩ := geometry.idealSurplusDomination x hx
    have gap := script.terminalSurplus y hy
    omega
  have idealSum := sum_charge_gaps a.idealCharges (c.scale * 36) ideal
  have shape := geometry.shapeLower
  rw [← index.exactDeficit] at shape
  have upper := geometry.paymentUpper
  have strict : 0 < c.scale * 4205 :=
    Nat.mul_pos script.scalePositive (by decide)
  omega

/-- The internally completed higher-dimensional tail is retained as an
explicit premise. It is not relabeled a published theorem. -/
structure Earlier58Tail (n : Nat) : Prop where
  atMost58 : n ≤ 58

/-- Actual rank58 geometry, not an assertion that n≤57. -/
structure Rank58Bridge (n : Nat) (c : Certificate) : Prop where
  chargesAt58 : n = 58 →
    ∃ a : ActualCharges, ∃ index : PublishedIndex a,
      WrittenGeometry c a index

theorem finite_volume_dimension_at_most57 (n : Nat) (c : Certificate)
    (script : ScriptChecks c) (earlier : Earlier58Tail n)
    (geometry : Rank58Bridge n c) : n ≤ 57 := by
  have excluded : n ≠ 58 := by
    intro at58
    obtain ⟨a, index, actual⟩ := geometry.chargesAt58 at58
    exact no_rank58_charges c script a index actual
  exact GuardedRank58.ceiling57_from_tail_and_exclusion n earlier.atMost58 excluded

end CoxeterProofs.ManuscriptFinite57

#print axioms CoxeterProofs.ManuscriptFinite57.no_rank58_charges
#print axioms CoxeterProofs.ManuscriptFinite57.finite_volume_dimension_at_most57
