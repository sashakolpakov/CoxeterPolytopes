import Std

/-!
Conditional logic for the sharp right-angled ceilings eight and four.

Only actual facet counts are natural numbers. `heightBound k F q` is
an ABSTRACT REAL-GEOMETRIC predicate: for k >= 1 all intrinsic unit-wall
heights are strictly below 1/sqrt(k). No real height, tensor norm or
polynomial pair sum is replaced by an integer variable.

The module checks small exact polynomial coefficients, the concrete
four-vector Gram contradiction, finite cardinality logic, composition
of three actual projection steps, and the two dimension implications.
Lorentz geometry, actual face inheritance, the tensor-square identity
and its real equality case remain explicit written interfaces.
The published Dufour input is separate from Musin/Kolpakov, which
are used only by the entirely-ideal ceiling-four implication.

Sources: right_angled_dimension8_descending_height_proof_20261020.md;
review_right_angled_dimension8_descending_height_batch_guard_20261020.md.
Source prepared without compilation. Only the coordinating agent may
schedule the serial resource-gated module compilation and axiom audit.
-/

namespace CoxeterProofs.RightAngledSharpDimensions

universe u v

/- Small exact coefficient certificates, in increasing degree order. -/
def polynomialAdd : List Int → List Int → List Int
  | [], q => q
  | p, [] => p
  | a :: p, b :: q => (a + b) :: polynomialAdd p q

def polynomialScale (a : Int) (p : List Int) : List Int :=
  p.map (fun b => a * b)

def polynomialMul : List Int → List Int → List Int
  | [], _ => []
  | a :: p, q =>
      polynomialAdd (polynomialScale a q) (0 :: polynomialMul p q)

/-- 16 F(x) = (2x+1)^2 (4x-1) = 16x^3+12x^2-1. -/
theorem cubic_factor_coefficients :
    polynomialMul (polynomialMul [1, 2] [1, 2]) [-1, 4]
      = [-1, 0, 12, 16] := by decide

/-- With P2=(6x^2-1)/5 and P3=(8x^3-3x)/5,
16F = P0+6P1+10P2+10P3. The real PSD assertion is separate. -/
theorem cubic_kernel_coefficients :
    polynomialAdd [1]
      (polynomialAdd (polynomialScale 6 [0, 1])
        (polynomialAdd (polynomialScale 2 [-1, 0, 6])
          (polynomialScale 2 [0, -3, 0, 8])))
      = [-1, 0, 12, 16] := by decide

theorem cubic_diagonal_and_positive_coefficients :
    ([-1, 0, 12, 16] : List Int).sum = 27 ∧
      (0 : Int) < 1 ∧ (0 : Int) < 6 ∧ (0 : Int) < 10 := by decide

/-- Twice the Gram matrix of four unit vectors with pair products -1/2. -/
def doubledFourHalfGram (i j : Nat) : Int := if i = j then 2 else -1

def fourHalfAllOnes : Int :=
  ((List.range 4).map (fun i =>
    ((List.range 4).map (fun j => doubledFourHalfGram i j)).sum)).sum

theorem four_half_all_ones_value : fourHalfAllOnes = -4 := by decide

/-- Real Euclidean Gram positivity supplies this particular integer
quadratic-form consequence. It is not an integer proxy for a real norm. -/
structure FourHalfGramPositive : Prop where
  allOnesNonnegative : 0 ≤ fourHalfAllOnes

theorem no_four_half_gram_positive (h : FourHalfGramPositive) : False := by
  have nonnegative := h.allOnesNonnegative
  have value := four_half_all_ones_value
  omega

/-- Exact integer consequences of the WRITTEN real tensor proof for
an N-point strict S5 code. The second field is the Gram-PSD consequence
of selecting four points if equality in the real pair sum holds.
Neither the final code bound nor the dimension conclusion is assumed. -/
structure StrictCodeSixConsequences (N : Nat) : Prop where
  tensorCount : N * N ≤ 27 * N
  equalityFourGram : N = 27 → FourHalfGramPositive

theorem strict_code_count_le27 (N : Nat) (h : StrictCodeSixConsequences N) :
    N ≤ 27 := by
  by_cases zero : N = 0
  · omega
  · exact Nat.le_of_mul_le_mul_right h.tensorCount (by omega)

theorem strict_code_count_le26 (N : Nat) (h : StrictCodeSixConsequences N) :
    N ≤ 26 := by
  have bound := strict_code_count_le27 N h
  have endpoint : N ≠ 27 := by
    intro equal
    exact no_four_half_gram_positive (h.equalityFourGram equal)
  omega

/-- Typed geometric predicates, to be instantiated by actual hyperbolic
faces and finite points. The subscript in a height bound counts descents,
not a dimension or a discretized real height. -/
structure GeometricData (Face : Type u) (Point : Type v) where
  finiteVolumeRightAngled : Face → Nat → Prop
  entirelyIdeal : Face → Prop
  relativeInterior : Face → Point → Prop
  actualFacet : Face → Face → Prop
  heightBound : Nat → Face → Point → Prop
  radialCodeFour : Face → Point → Prop
  radialCodeSix : Face → Point → Prop
  facetCount : Face → Nat
  regular24 : Face → Prop
  atRegularCenter : Face → Point → Prop
  allHeightsExactlyOne : Face → Point → Prop

/-- Written actual-wall and real-height inputs for the farthest-facet
lemma. `farthestProjection` supplies a REAL projection, its actual facet,
strict interior, inherited finite volume/right angles and the next
strict real bound. It is used only in dimensions at least four. -/
structure ProjectionGeometry {Face : Type u} {Point : Type v}
    (d : GeometricData Face Point) : Prop where
  interiorExists : ∀ F n, d.finiteVolumeRightAngled F n →
    ∃ q, d.relativeInterior F q
  finiteStart : ∀ F q, d.relativeInterior F q → d.heightBound 0 F q
  farthestProjection : ∀ k n F q, 3 ≤ n →
    d.finiteVolumeRightAngled F (n + 1) →
    d.relativeInterior F q → d.heightBound k F q →
    ∃ Q p, d.actualFacet F Q ∧ d.finiteVolumeRightAngled Q n ∧
      d.relativeInterior Q p ∧ d.heightBound (k + 1) Q p
  idealFacetInheritance : ∀ F Q, d.actualFacet F Q →
    d.entirelyIdeal F → d.entirelyIdeal Q

/-- Intrinsic unit-normal and all-wall sign reductions. The code
predicates include actual-neighbor cardinality and injective directions. -/
structure RadialGeometry {Face : Type u} {Point : Type v}
    (d : GeometricData Face Point) : Prop where
  fourCode : ∀ F q, d.finiteVolumeRightAngled F 4 →
    d.relativeInterior F q → d.heightBound 1 F q → d.radialCodeFour F q
  sixCode : ∀ F q, d.finiteVolumeRightAngled F 6 →
    d.relativeInterior F q → d.heightBound 3 F q → d.radialCodeSix F q

/-- The universal real tensor-square proof and its strict endpoint.
No optimizer, spherical-code classification or published kissing bound
is used for this six-dimensional code input. -/
structure WrittenCubicGeometry {Face : Type u} {Point : Type v}
    (d : GeometricData Face Point) : Prop where
  tensorAndStrictEquality : ∀ F q, d.radialCodeSix F q →
    StrictCodeSixConsequences (d.facetCount F)

/-- Dufour2010, Lemma1 and Proposition4: all finite-volume RA6
polytopes have at least27 facets, irrespective of vertex types. -/
structure PublishedDufour {Face : Type u} {Point : Type v}
    (d : GeometricData Face Point) : Prop where
  minimumSixFacets : ∀ F, d.finiteVolumeRightAngled F 6 → 27 ≤ d.facetCount F

/-- Musin cardinality and Kolpakov's journal Theorem4, used ONLY by
the entirely-ideal bound. Metric equality is distinct from code uniqueness. -/
structure PublishedIdealFour {Face : Type u} {Point : Type v}
    (d : GeometricData Face Point) : Prop where
  musin : ∀ F q, d.radialCodeFour F q → d.facetCount F ≤ 24
  kolpakovMinimum : ∀ F, d.finiteVolumeRightAngled F 4 →
    d.entirelyIdeal F → 24 ≤ d.facetCount F
  kolpakovEquality : ∀ F, d.finiteVolumeRightAngled F 4 →
    d.entirelyIdeal F → d.facetCount F = 24 → d.regular24 F

/-- The explicit opposite-normal-pair metric argument for the regular
24-cell, and the incompatibility of its unit heights with a strict bound. -/
structure WrittenRegularCellGeometry {Face : Type u} {Point : Type v}
    (d : GeometricData Face Point) : Prop where
  oppositePairsForceCenter : ∀ F q, d.regular24 F →
    d.relativeInterior F q → d.heightBound 1 F q → d.atRegularCenter F q
  centerHasUnitHeights : ∀ F q, d.atRegularCenter F q → d.allHeightsExactlyOne F q
  strictVersusUnit : ∀ F q, d.heightBound 1 F q →
    d.allHeightsExactlyOne F q → False

theorem actual_three_projection_flag {Face : Type u} {Point : Type v}
    (d : GeometricData Face Point) (projection : ProjectionGeometry d)
    (F : Face) (nine : d.finiteVolumeRightAngled F 9) :
    ∃ R z, d.finiteVolumeRightAngled R 6 ∧
      d.relativeInterior R z ∧ d.heightBound 3 R z := by
  obtain ⟨q, interior⟩ := projection.interiorExists F 9 nine
  have start := projection.finiteStart F q interior
  obtain ⟨F8, q8, _, eight, interior8, bound8⟩ :=
    projection.farthestProjection 0 8 F q (by decide) nine interior start
  obtain ⟨F7, q7, _, seven, interior7, bound7⟩ :=
    projection.farthestProjection 1 7 F8 q8 (by decide) eight interior8 bound8
  obtain ⟨F6, q6, _, six, interior6, bound6⟩ :=
    projection.farthestProjection 2 6 F7 q7 (by decide) seven interior7 bound7
  exact ⟨F6, q6, six, interior6, bound6⟩

theorem no_nine_polytope {Face : Type u} {Point : Type v}
    (d : GeometricData Face Point) (projection : ProjectionGeometry d)
    (radial : RadialGeometry d) (cubic : WrittenCubicGeometry d)
    (published : PublishedDufour d) (F : Face)
    (nine : d.finiteVolumeRightAngled F 9) : False := by
  obtain ⟨R, z, six, interior, bound⟩ := actual_three_projection_flag d projection F nine
  have code := radial.sixCode R z six interior bound
  have upper := strict_code_count_le26 _ (cubic.tensorAndStrictEquality R z code)
  have lower := published.minimumSixFacets R six
  omega

theorem all_right_angled_dimension_at_most_eight
    {Face : Type u} {Point : Type v}
    (d : GeometricData Face Point) (projection : ProjectionGeometry d)
    (radial : RadialGeometry d) (cubic : WrittenCubicGeometry d)
    (published : PublishedDufour d) (dimension : Nat)
    (nineFaceInheritance : 9 ≤ dimension →
      ∃ F, d.finiteVolumeRightAngled F 9) : dimension ≤ 8 := by
  have excluded : ¬ 9 ≤ dimension := by
    intro large
    obtain ⟨F, nine⟩ := nineFaceInheritance large
    exact no_nine_polytope d projection radial cubic published F nine
  omega

theorem no_ideal_five_polytope {Face : Type u} {Point : Type v}
    (d : GeometricData Face Point) (projection : ProjectionGeometry d)
    (radial : RadialGeometry d) (published : PublishedIdealFour d)
    (metric : WrittenRegularCellGeometry d) (F : Face)
    (five : d.finiteVolumeRightAngled F 5) (ideal : d.entirelyIdeal F) : False := by
  obtain ⟨q, interior⟩ := projection.interiorExists F 5 five
  have start := projection.finiteStart F q interior
  obtain ⟨Q, p, facet, four, interior4, bound⟩ :=
    projection.farthestProjection 0 4 F q (by decide) five interior start
  have ideal4 := projection.idealFacetInheritance F Q facet ideal
  have code := radial.fourCode Q p four interior4 bound
  have upper := published.musin Q p code
  have lower := published.kolpakovMinimum Q four ideal4
  have count : d.facetCount Q = 24 := by omega
  have regular := published.kolpakovEquality Q four ideal4 count
  have center := metric.oppositePairsForceCenter Q p regular interior4 bound
  exact metric.strictVersusUnit Q p bound (metric.centerHasUnitHeights Q p center)

theorem entirely_ideal_dimension_at_most_four
    {Face : Type u} {Point : Type v}
    (d : GeometricData Face Point) (projection : ProjectionGeometry d)
    (radial : RadialGeometry d) (published : PublishedIdealFour d)
    (metric : WrittenRegularCellGeometry d) (dimension : Nat)
    (fiveFaceInheritance : 5 ≤ dimension →
      ∃ F, d.finiteVolumeRightAngled F 5 ∧ d.entirelyIdeal F) : dimension ≤ 4 := by
  have excluded : ¬ 5 ≤ dimension := by
    intro large
    obtain ⟨F, five, ideal⟩ := fiveFaceInheritance large
    exact no_ideal_five_polytope d projection radial published metric F five ideal
  omega

end CoxeterProofs.RightAngledSharpDimensions

#print axioms CoxeterProofs.RightAngledSharpDimensions.cubic_factor_coefficients
#print axioms CoxeterProofs.RightAngledSharpDimensions.cubic_kernel_coefficients
#print axioms CoxeterProofs.RightAngledSharpDimensions.cubic_diagonal_and_positive_coefficients
#print axioms CoxeterProofs.RightAngledSharpDimensions.four_half_all_ones_value
#print axioms CoxeterProofs.RightAngledSharpDimensions.no_four_half_gram_positive
#print axioms CoxeterProofs.RightAngledSharpDimensions.strict_code_count_le27
#print axioms CoxeterProofs.RightAngledSharpDimensions.strict_code_count_le26
#print axioms CoxeterProofs.RightAngledSharpDimensions.actual_three_projection_flag
#print axioms CoxeterProofs.RightAngledSharpDimensions.no_nine_polytope
#print axioms CoxeterProofs.RightAngledSharpDimensions.all_right_angled_dimension_at_most_eight
#print axioms CoxeterProofs.RightAngledSharpDimensions.no_ideal_five_polytope
#print axioms CoxeterProofs.RightAngledSharpDimensions.entirely_ideal_dimension_at_most_four
