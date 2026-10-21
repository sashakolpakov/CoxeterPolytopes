import CoxeterProofs.CertificateBounds

/-!
The complete integer case analysis in the second-volume argument.

`FacetCase` describes the conclusions of the geometric/three-dimensional
facet reductions. Those reductions are not yet formalized; this module does
not turn them into axioms or claim an end-to-end geometric formalization.
-/
namespace CoxeterProofs.IdealRA4

inductive FacetCase (m q s : Nat) : Prop where
  | antiprism (minimum : 2 ≤ m) (faces : q = 2) (cusps : s = m + 6)
  | defectThree (defect : m = 3) (faces : q = 3) (cusps : s = 3)
  | defectFour (defect : m = 4) (faces : 4 ≤ q)
  | largeDefect (defect : 5 ≤ m) (faces : 3 ≤ q)

theorem case_contribution_at_least_twenty
    {m q s : Nat} (shape : FacetCase m q s) :
    20 ≤ m * (1 + q + s) := by
  cases shape with
  | antiprism hm hq hs =>
    have hf : 11 ≤ 1 + q + s := by omega
    have lower := Nat.mul_le_mul hm hf
    omega
  | defectThree hm hq hs =>
    subst m q s
    decide
  | defectFour hm hq =>
    have hf : 5 ≤ 1 + q + s := by omega
    have lower := Nat.mul_le_mul (Nat.le_refl 4) hf
    subst m
    omega
  | largeDefect hm hq =>
    have hf : 4 ≤ 1 + q + s := by omega
    have lower := Nat.mul_le_mul hm hf
    omega

theorem defect_gap
    {W m q s : Nat}
    (divisibility : W % 3 = 0)
    (budget : m * (1 + q + s) ≤ 2 * W)
    (shape : FacetCase m q s) : 12 ≤ W := by
  have lower := case_contribution_at_least_twenty shape
  omega

theorem defect_gap_of_partition
    {W m q s : Nat} (selected rest : List Nat)
    (divisibility : W % 3 = 0)
    (partition : selected.sum + rest.sum = 2 * W)
    (count : selected.length = 1 + q + s)
    (minimum : ∀ a ∈ selected, m ≤ a)
    (shape : FacetCase m q s) : 12 ≤ W := by
  exact defect_gap divisibility
    (CoxeterProofs.selected_defect_budget W m q s selected rest partition count minimum)
    shape

theorem intermediate_defects_excluded
    {W m q s : Nat}
    (divisibility : W % 3 = 0)
    (budget : m * (1 + q + s) ≤ 2 * W)
    (shape : FacetCase m q s) : W ≠ 3 ∧ W ≠ 6 ∧ W ≠ 9 := by
  have gap := defect_gap divisibility budget shape
  omega

/- `numerator / 3` is volume in units pi^2; the geometric volume formula
   supplies W = 3*(numerator-4). No real-analysis axiom is added here. -/
theorem volume_numerator_at_least_eight
    {W numerator : Nat} (gap : 12 ≤ W)
    (volumeIdentity : W = 3 * (numerator - 4)) : 8 ≤ numerator := by
  omega

theorem double_has_gap_value : 3 * (42 - 38) = 12 := by decide
theorem double_has_volume_numerator : 42 - 38 + 4 = 8 := by decide

end CoxeterProofs.IdealRA4
