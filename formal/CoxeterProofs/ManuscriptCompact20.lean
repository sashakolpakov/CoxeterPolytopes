import CoxeterProofs.CompactArithmetic20
import CoxeterProofs.ManuscriptDataBounds

/-!
Conditional compact20 composition. The exact rational certificate is
represented after clearing one positive denominator. Its finite rows,
actual-star coverage and the older field proofs are separate parameters.
The theorem does not assume the desired dimension bound.
-/
namespace CoxeterProofs.ManuscriptCompact20

open CompactHVector ManuscriptDataBounds

inductive GroundField21 where
  | rational | sqrtTwo | sqrtThree | sqrtFive | sqrtSix
  | sqrtTwoThree | sqrtTwoFive
  | cubicSeven | cyclotomicNine | cyclotomicEleven
  | cyclotomicFifteen | cyclotomicSixteen | cyclotomicTwenty
  deriving DecidableEq, Repr

/-- Finite exact outputs of the saved rank21 charge certificate.
The forest offers are complete-recurrence allowances, not only the
finite support rows selected by a numerical optimizer. -/
structure Certificate where
  denominator : Nat
  nativePayments : List Nat
  forestOffers : List Nat

/-- Admitted exact script conclusions, distinct from actual-domain coverage. -/
structure ScriptChecks (c : Certificate) : Prop where
  denominatorPositive : 0 < c.denominator
  nativePaid : ∀ p ∈ c.nativePayments, c.denominator ≤ p
  forestUpper : ∀ w ∈ c.forestOffers, w ≤ c.denominator * 40

/-- Written geometry/h-vector interfaces for one actual golden rank21
chamber. Native coverage, whole components, the D5/H4 exceptions and
angle double counting are required to construct this witness. -/
structure Rank21Geometry (c : Certificate) where
  tail : List Term
  vertexCharges : List Nat
  pairedCoefficientLower :
    ∀ t ∈ tail, 40 * t.multiplicity ≤ t.coefficient
  vertexCount : vertexCharges.length = vertices (CompactArithmetic20.endpoint21 :: tail)
  completeForestDomination :
    ∀ a ∈ vertexCharges, ∃ b ∈ c.forestOffers, a ≤ b
  coveredFaceIncidence :
    (∀ p ∈ c.nativePayments, c.denominator ≤ p) →
    c.denominator * weighted (CompactArithmetic20.endpoint21 :: tail) ≤ vertexCharges.sum

theorem no_rank21_geometry (c : Certificate) (script : ScriptChecks c)
    (geometry : Rank21Geometry c) : False := by
  have charges : ∀ a ∈ geometry.vertexCharges, a ≤ c.denominator * 40 := by
    intro a ha
    obtain ⟨b, hb, hab⟩ := geometry.completeForestDomination a ha
    exact Nat.le_trans hab (script.forestUpper b hb)
  have upper := sum_le_cap_mul_length geometry.vertexCharges
    (c.denominator * 40) charges
  rw [geometry.vertexCount, Nat.mul_assoc] at upper
  exact CompactArithmetic20.no_rank21_weight_model geometry.tail
    c.denominator geometry.vertexCharges.sum script.denominatorPositive
    geometry.pairedCoefficientLower
    (geometry.coveredFaceIncidence script.nativePaid) upper

/-- Published compact29 and the thirteen-field classification valid at21.
The predicate must mean equality of the actual arithmetic ground field,
not containment in a larger field. -/
structure PublishedInputs (n : Nat) (actualField : GroundField21 → Prop) : Prop where
  compactCeiling29 : n ≤ 29
  thirteenFieldsAt21 : n = 21 → ∃ k : GroundField21, actualField k

/-- Previously completed internal proofs. These are explicitly hypotheses,
not published facts: the24 higher branches and12 nongolden branches. -/
structure EarlierInternalProofs (n : Nat) (actualField : GroundField21 → Prop) : Prop where
  higherBranches : 22 ≤ n → n ≤ 29 → False
  nongoldenAt21 :
    n = 21 → ∀ k, actualField k → k ≠ .sqrtFive → False

/-- The new actual-star/whole-forest bridge is needed only for the golden
dimension21 branch. Its conclusion is charge data, not a dimension bound. -/
structure GoldenGeometryBridge (n : Nat) (actualField : GroundField21 → Prop)
    (c : Certificate) : Prop where
  actualRank21Data : n = 21 → actualField .sqrtFive → Nonempty (Rank21Geometry c)

theorem compact_arithmetic_dimension_at_most20
    (n : Nat) (actualField : GroundField21 → Prop) (c : Certificate)
    (published : PublishedInputs n actualField)
    (earlier : EarlierInternalProofs n actualField)
    (script : ScriptChecks c)
    (geometry : GoldenGeometryBridge n actualField c) : n ≤ 20 := by
  by_cases low : n ≤ 20
  · exact low
  · have contradiction : False := by
      by_cases at21 : n = 21
      · obtain ⟨k, hk⟩ := published.thirteenFieldsAt21 at21
        by_cases golden : k = .sqrtFive
        · subst k
          obtain ⟨witness⟩ := geometry.actualRank21Data at21 hk
          exact no_rank21_geometry c script witness
        · exact earlier.nongoldenAt21 at21 k hk golden
      · exact earlier.higherBranches (by omega) published.compactCeiling29
    exact False.elim contradiction

end CoxeterProofs.ManuscriptCompact20

#print axioms CoxeterProofs.ManuscriptCompact20.no_rank21_geometry
#print axioms CoxeterProofs.ManuscriptCompact20.compact_arithmetic_dimension_at_most20
