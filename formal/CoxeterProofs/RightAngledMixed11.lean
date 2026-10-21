import Std

/-!
Final count logic for the reviewed positive-part opposition argument.
Actual Lorentz geometry, real-valued row deletion/reciprocal summation,
the script-supported full-interval polynomial majorant and the published
dimension12 tail remain explicit premises. This module does not claim
that actual analytic rows are integer-valued, or formalize their geometry.
-/
namespace CoxeterProofs.RightAngledMixed11

/-- All fields are actual integer face counts or the specified signed
excess/defect expressions. The exceptional-incidence comparisons are
cleared forms of the complete written local arguments. -/
structure Counts12 where
  finiteVertices : Int
  cusps : Int
  twoFaces : Int
  threeFaces : Int
  bipyramids : Int
  excess : Int
  defect : Int
  exceptional : Int
  finiteNonnegative : 0 ≤ finiteVertices
  strictWeighted : 15 * finiteVertices + 4246 * cusps < twoFaces
  cofaceEuler : threeFaces = 5 * twoFaces - 55 * finiteVertices - 14080 * cusps
  excessIdentity : excess = bipyramids - 10 * threeFaces + 10 * twoFaces + 28160 * cusps
  defectIdentity : defect = 8 * threeFaces - 10 * twoFaces - 2 * bipyramids
  excessNonnegative : 0 ≤ excess
  exceptionalUpper : exceptional ≤ 22 * excess + 4 * defect
  exceptionalLower : 13 * (220 * finiteVertices - 2 * bipyramids) ≤ 27 * exceptional

theorem polygon_coupling (d : Counts12) :
    5130 * d.twoFaces ≤
      101 * d.bipyramids + 69080 * d.finiteVertices + 22049280 * d.cusps := by
  have a := d.cofaceEuler
  have e := d.excessIdentity
  have r := d.defectIdentity
  have upper := d.exceptionalUpper
  have lower := d.exceptionalLower
  omega

theorem canonical_sum_lower (d : Counts12) :
    50 * d.finiteVertices + 880 * d.cusps < d.bipyramids := by
  have weighted := d.strictWeighted
  have incidence := d.cofaceEuler
  have identity := d.excessIdentity
  have nonnegative := d.excessNonnegative
  omega

theorem strict_opposition_lower (d : Counts12) :
    184460 * d.finiteVertices < 1619 * (2 * d.bipyramids) := by
  have coupling := polygon_coupling d
  have weighted := d.strictWeighted
  have canonical := canonical_sum_lower d
  omega

/-- No geometry or script output is assumed by these constant identities. -/
theorem exact_comparison :
    (440 * 12500 : Int) = 5500000 ∧
      (184460 * 12500 - 1420561 * 1619 : Int) = 5861741 ∧
      (0 : Int) < 5861741 := by decide

/-- The premise is the cleared aggregate analytic comparison. Its written
proof uses the independently replayed quartic majorant, all220 formal
triple moments, nonnegative deletion and actual undirected reciprocal
pairing. It is not an assumption of a dimension conclusion. -/
structure Data12 where
  counts : Counts12
  positivePartTransport :
    12500 * (2 * counts.bipyramids) ≤ 1420561 * counts.finiteVertices

theorem no_data12 (d : Data12) : False := by
  have lower := strict_opposition_lower d.counts
  have upper := d.positivePartTransport
  have nonnegative := d.counts.finiteNonnegative
  omega

theorem dimension_at_most11 (n : Nat)
    (publishedDufour : n ≤ 12)
    (actualGeometry : n = 12 → Nonempty Data12) : n ≤ 11 := by
  have excluded : n ≠ 12 := by
    intro equal
    obtain ⟨d⟩ := actualGeometry equal
    exact no_data12 d
  omega

end CoxeterProofs.RightAngledMixed11

#print axioms CoxeterProofs.RightAngledMixed11.polygon_coupling
#print axioms CoxeterProofs.RightAngledMixed11.canonical_sum_lower
#print axioms CoxeterProofs.RightAngledMixed11.strict_opposition_lower
#print axioms CoxeterProofs.RightAngledMixed11.exact_comparison
#print axioms CoxeterProofs.RightAngledMixed11.no_data12
#print axioms CoxeterProofs.RightAngledMixed11.dimension_at_most11
