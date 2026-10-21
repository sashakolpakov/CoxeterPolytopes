import Std

/-!
Conditional integer logic for the dimension-ten right-angled argument.
Actual face geometry, the external compact/ideal bounds, the real subset
inequalities and their exact script certificates remain explicit interfaces.
No assertion that real-valued ratios or rows are integers is made here.
The analytic transport field is their cleared aggregate count consequence.
-/
namespace CoxeterProofs.RightAngledMixed10

/-- Actual ordinary counts in a hypothetical mixed eleven-polytope.
`defect` and `balance` are signed quantities; neither is assumed nonnegative.
The local geometric arguments supply the stated identities and bounds. -/
structure Counts11 where
  finiteVertices : Int
  cusps : Int
  twoFaces : Int
  threeFaces : Int
  bipyramids : Int
  excess : Int
  defect : Int
  exceptional : Int
  balance : Int
  noncanonicalIdealCorners : Int
  finitePositive : 0 < finiteVertices
  excessNonnegative : 0 ≤ excess
  strictWeighted : 25 * finiteVertices + 3860 * cusps < 2 * twoFaces
  cofaceEuler : 4 * threeFaces = 18 * twoFaces - 165 * finiteVertices - 23040 * cusps
  excessIdentity : excess = bipyramids - 10 * threeFaces + 9 * twoFaces + 11520 * cusps
  defectIdentity : defect = 8 * threeFaces - 9 * twoFaces - 2 * bipyramids
  exceptionalUpper : exceptional ≤ 18 * excess + 4 * defect
  balanceLower : -30 * excess ≤ balance
  noncanonicalUpper : noncanonicalIdealCorners ≤ 2 * excess
  summedCornerBalance :
    16 * (165 * finiteVertices - 2 * bipyramids) + 3 * balance
      - 16 * noncanonicalIdealCorners ≤ 31 * exceptional

theorem balanced_exceptional_lower (d : Counts11) :
    16 * (165 * d.finiteVertices - 2 * d.bipyramids) - 122 * d.excess
      ≤ 31 * d.exceptional := by
  have corner := d.summedCornerBalance
  have balance := d.balanceLower
  have noncanonical := d.noncanonicalUpper
  omega

theorem exact_defect_identity (d : Counts11) :
    8 * d.defect = -10 * d.bipyramids - 165 * d.finiteVertices
      + 46080 * d.cusps - 6 * d.excess := by
  have coface := d.cofaceEuler
  have excess := d.excessIdentity
  have defect := d.defectIdentity
  omega

theorem strict_cusp_row (d : Counts11) :
    720 * d.cusps < 2 * d.bipyramids - 75 * d.finiteVertices - 2 * d.excess := by
  have weighted := d.strictWeighted
  have coface := d.cofaceEuler
  have excess := d.excessIdentity
  omega

theorem strict_opposition_lower (d : Counts11) :
    159195 * d.finiteVertices + 2794 * d.excess < 3722 * d.bipyramids := by
  have lower := balanced_exceptional_lower d
  have upper := d.exceptionalUpper
  have defect := exact_defect_identity d
  have cusp := strict_cusp_row d
  omega

/-- A cleared ordinary-count consequence of the real analytic proof.
It follows from every subset-size contract and pairing reverse actual
opposition ratios. Its proof is written, not formalized in this module. -/
structure Data11 where
  counts : Counts11
  subsetTransport :
    8 * counts.bipyramids ≤ 336 * counts.finiteVertices
      + 3 * counts.noncanonicalIdealCorners

theorem opposition_upper (d : Data11) :
    4 * d.counts.bipyramids ≤ 168 * d.counts.finiteVertices + 3 * d.counts.excess := by
  have transport := d.subsetTransport
  have noncanonical := d.counts.noncanonicalUpper
  omega

theorem exact_comparison :
    (2 * 159195 - 1861 * 168 : Int) = 5742 ∧
      (2 * 2794 - 1861 * 3 : Int) = 5 ∧
      (0 : Int) < 5742 ∧ (0 : Int) < 5 := by decide

theorem no_data11 (d : Data11) : False := by
  have lower := strict_opposition_lower d.counts
  have upper := opposition_upper d
  have positive := d.counts.finitePositive
  have nonnegative := d.counts.excessNonnegative
  omega

/-- The interface includes face heredity and the published smaller
compact/entirely-ideal bounds: an eleven-face would be mixed and would
provide actual data. No previous internal ceiling is needed as a premise. -/
theorem dimension_at_most10 (n : Nat)
    (elevenFaceGeometry : 11 ≤ n → Nonempty Data11) : n ≤ 10 := by
  have excluded : ¬ 11 ≤ n := by
    intro large
    obtain ⟨d⟩ := elevenFaceGeometry large
    exact no_data11 d
  omega

end CoxeterProofs.RightAngledMixed10

#print axioms CoxeterProofs.RightAngledMixed10.balanced_exceptional_lower
#print axioms CoxeterProofs.RightAngledMixed10.exact_defect_identity
#print axioms CoxeterProofs.RightAngledMixed10.strict_cusp_row
#print axioms CoxeterProofs.RightAngledMixed10.strict_opposition_lower
#print axioms CoxeterProofs.RightAngledMixed10.opposition_upper
#print axioms CoxeterProofs.RightAngledMixed10.exact_comparison
#print axioms CoxeterProofs.RightAngledMixed10.no_data11
#print axioms CoxeterProofs.RightAngledMixed10.dimension_at_most10
