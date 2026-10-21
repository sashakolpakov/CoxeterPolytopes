import CoxeterProofs.RightAngledCuspPropagation

/-!
Current conditional arithmetic behind the ideal5/6 manuscript.
The actual small-facet bridge, cusp geometry, spinor defect, eleven-vertex
equality exclusion, facet injections and reciprocal transport are explicit
inputs, not published axioms. No full hyperbolic geometry is formalized.
Only the earlier count structures and primitive propagation implication
are reused; none of their superseded final numerical bounds is invoked.
-/
namespace CoxeterProofs.RightAngledElevenVertexCounts

/-- The inherited structure contains the actual positive carrier, its
three-facet with at least11 vertices, the +3 defect counts and propagation
counts. This NEW internal geometric input excludes its equality case36
when the selected three-facet has exactly11 vertices. It does not assume
the new37/63/211 conclusions. -/
structure Ideal5Data where
  propagation : RightAngledCuspPropagation.PropagatingIdeal5Data
  elevenVertexEqualityExcluded :
    propagation.base.base.threeFacetVertices = 11 →
      propagation.base.base.carrierVertices ≠ 36

abbrev Ideal5Data.counts (d : Ideal5Data) : RightAngledCuspCapacity.Ideal5Data :=
  d.propagation.base.base

theorem carrier_bounds (d : Ideal5Data) :
    37 ≤ d.counts.carrierVertices ∧ 25 ≤ d.counts.carrierFacets := by
  have defect : 3 * d.counts.threeFacetVertices + 3 ≤ d.counts.carrierVertices :=
    (RightAngledCuspDefect.stronger_carrier_bound d.propagation.base).2
  have large := d.counts.largeThreeFacet
  have injection := d.counts.carrierInjection
  constructor
  · by_cases h : d.counts.threeFacetVertices = 11
    · have excluded : d.counts.carrierVertices ≠ 36 :=
        d.elevenVertexEqualityExcluded h
      omega
    · omega
  · omega

theorem ideal5_counts (d : Ideal5Data) :
    63 ≤ d.counts.facetCounts.length ∧
      631 ≤ 3 * d.counts.vertices ∧ 211 ≤ d.counts.vertices := by
  have carrier := carrier_bounds d
  have injection := d.counts.actualFacetInjection
  have propagation := RightAngledCuspPropagation.facet_vertex_propagation
    d.counts.vertices d.counts.carrierVertices d.counts.carrierFacets
    d.propagation.noNeighbors d.propagation.exceptionalEdges
    d.propagation.smallCorners d.propagation.goodCusps
    d.propagation.outsideIdentity d.propagation.polygonCornerLower
    d.propagation.pairedCornerUpper d.propagation.actualOppositeIncidence
  omega

theorem sum_ideal5_facet_minima (xs : List Ideal5Data) :
    63 * xs.length ≤ (xs.map (fun d => d.counts.facetCounts.length)).sum := by
  induction xs with
  | nil => simp
  | cons a xs ih =>
    have bound := (ideal5_counts a).1
    simp only [List.length_cons, List.map_cons, List.sum_cons]
    omega

/-- Every actual five-facet has its own explicit interface above. The
injection and coface equality concern ALL these facets. Reciprocal transport
is the internal written universal good-edge/103/132 argument after clearing
denominators, not a finite sample or an assumed final size bound. -/
structure Ideal6Data where
  fiveFacets : List Ideal5Data
  vertices : Nat
  fourFaces : Nat
  actualFacetInjection : ∃ d ∈ fiveFacets,
    1 + d.counts.facetCounts.length + d.counts.vertices ≤ fiveFacets.length
  cofaceIncidence :
    (fiveFacets.map (fun d => d.counts.facetCounts.length)).sum = 2 * fourFaces
  reciprocalTransport : 99 * fourFaces ≤ 99 * fiveFacets.length + 136 * vertices

theorem ideal6_counts (d : Ideal6Data) :
    275 ≤ d.fiveFacets.length ∧
      6039 * d.fiveFacets.length ≤ 272 * d.vertices ∧ 6106 ≤ d.vertices := by
  obtain ⟨q, _, injection⟩ := d.actualFacetInjection
  have sizes := ideal5_counts q
  have cofaceLower := sum_ideal5_facet_minima d.fiveFacets
  have incidence := d.cofaceIncidence
  have transport := d.reciprocalTransport
  omega

theorem exact_rounding :
    (13 * 37 + 6 * 25 : Nat) = 631 ∧
      3 * 210 < 631 ∧ 631 ≤ 3 * 211 ∧
      (6039 * 275 : Nat) = 1660725 ∧
      272 * 6105 < 1660725 ∧ 1660725 ≤ 272 * 6106 := by decide

theorem sum_ideal6_facet_minima (xs : List Ideal6Data) :
    275 * xs.length ≤ (xs.map (fun d => d.fiveFacets.length)).sum := by
  induction xs with
  | nil => simp
  | cons a xs ih =>
    have bound := (ideal6_counts a).1
    simp only [List.length_cons, List.map_cons, List.sum_cons]
    omega

theorem sum_minimum (xs : List Nat) (m : Nat)
    (minimum : ∀ a ∈ xs, m ≤ a) : m * xs.length ≤ xs.sum := by
  induction xs with
  | nil => simp
  | cons a xs ih =>
    have head := minimum a (by simp)
    have tail := ih (fun b hb => minimum b (by simp [hb]))
    simp only [List.length_cons, List.sum_cons, Nat.mul_add, Nat.mul_one]
    omega

/-- Intended for the complete face lists of a mixed right-angled polytope
in ambient dimension11 or12. Each ideal face has the internal structure
above. The other face minima16/27 and strict NK averages20/42 are the
separately stated PUBLISHED inputs, with their actual-face inheritance and
dimension range supplied by the manuscript. List completeness/realization
is a geometric premise; it is not proved by the arithmetic type. -/
structure MixedFaceData where
  idealFiveFaces : List Ideal5Data
  otherFiveFacetCounts : List Nat
  generalFiveMinimum : ∀ a ∈ otherFiveFacetCounts, 16 ≤ a
  fiveStrictAverage :
    (idealFiveFaces.map (fun d => d.counts.facetCounts.length)).sum +
      otherFiveFacetCounts.sum < 20 * (idealFiveFaces.length + otherFiveFacetCounts.length)
  idealSixFaces : List Ideal6Data
  otherSixFacetCounts : List Nat
  generalSixMinimum : ∀ a ∈ otherSixFacetCounts, 27 ≤ a
  sixStrictAverage :
    (idealSixFaces.map (fun d => d.fiveFacets.length)).sum +
      otherSixFacetCounts.sum < 42 * (idealSixFaces.length + otherSixFacetCounts.length)

/-- Exact cleared forms of I5/N5 < 4/47 and I6/N6 < 15/248.
Positive total counts are derived, so the fractional interpretation has
no zero-denominator exception. No dimension exclusion follows. -/
theorem mixed_density (d : MixedFaceData) :
    0 < d.idealFiveFaces.length + d.otherFiveFacetCounts.length ∧
      0 < d.idealSixFaces.length + d.otherSixFacetCounts.length ∧
      47 * d.idealFiveFaces.length <
        4 * (d.idealFiveFaces.length + d.otherFiveFacetCounts.length) ∧
      248 * d.idealSixFaces.length <
        15 * (d.idealSixFaces.length + d.otherSixFacetCounts.length) := by
  have idealFive := sum_ideal5_facet_minima d.idealFiveFaces
  have otherFive := sum_minimum d.otherFiveFacetCounts 16 d.generalFiveMinimum
  have averageFive := d.fiveStrictAverage
  have idealSix := sum_ideal6_facet_minima d.idealSixFaces
  have otherSix := sum_minimum d.otherSixFacetCounts 27 d.generalSixMinimum
  have averageSix := d.sixStrictAverage
  omega

end CoxeterProofs.RightAngledElevenVertexCounts

#print axioms CoxeterProofs.RightAngledElevenVertexCounts.carrier_bounds
#print axioms CoxeterProofs.RightAngledElevenVertexCounts.ideal5_counts
#print axioms CoxeterProofs.RightAngledElevenVertexCounts.sum_ideal5_facet_minima
#print axioms CoxeterProofs.RightAngledElevenVertexCounts.ideal6_counts
#print axioms CoxeterProofs.RightAngledElevenVertexCounts.exact_rounding
#print axioms CoxeterProofs.RightAngledElevenVertexCounts.sum_ideal6_facet_minima
#print axioms CoxeterProofs.RightAngledElevenVertexCounts.mixed_density
