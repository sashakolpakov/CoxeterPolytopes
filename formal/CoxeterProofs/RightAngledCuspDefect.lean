import CoxeterProofs.RightAngledCuspCapacity

/-!
Conditional counting logic for the new positive ideal-four cusp defect.
The spinor polygon lemma, actual opposite-wall geometry, and the exact
triangle/outgoing-edge correspondence are explicit count interfaces below.
No final facet/vertex bound and no dimension exclusion is an axiom.
-/
namespace CoxeterProofs.RightAngledCuspDefect

open RightAngledCuspCapacity

/-- `small` counts small polygon corners, `good` counts cusps at which both
paired width coefficients are small. The lower bound uses the polygon
lemma and at least8 triangular faces, at most2 per exceptional edge. -/
theorem positive_defect_count (m exceptional small good noNeighbors : Nat)
    (smallLower : 2 * m + 4 + (8 - 2 * exceptional) ≤ small)
    (smallUpper : small ≤ 2 * m + 2 * good)
    (oppositeCapacity : good ≤ 2 * noNeighbors) :
    1 ≤ noNeighbors ∧ 6 ≤ 2 * noNeighbors + exceptional ∧
      3 ≤ noNeighbors + exceptional := by omega

/-- The selected four-carrier in an ideal five-polytope, with its actual
outside-vertex counts and polygon-corner interfaces. `base` retains the
published minimum and complete actual facet/Euler data, not a hypothesis
that the stronger62/189 counts already hold. -/
structure DefectIdeal5Data where
  base : Ideal5Data
  outside : Nat
  noNeighbors : Nat
  oneNeighbor : Nat
  twoNeighbors : Nat
  exceptionalEdges : Nat
  smallCorners : Nat
  goodCusps : Nat
  verticesPartition : base.carrierVertices = base.threeFacetVertices + outside
  outsidePartition : outside = noNeighbors + oneNeighbor + twoNeighbors
  outgoingEdges : oneNeighbor + 2 * twoNeighbors = 4 * base.threeFacetVertices
  triangularBijection : exceptionalEdges + twoNeighbors = 2 * base.threeFacetVertices
  polygonCornerLower :
    2 * base.threeFacetVertices + 4 + (8 - 2 * exceptionalEdges) ≤ smallCorners
  pairedCornerUpper : smallCorners ≤ 2 * base.threeFacetVertices + 2 * goodCusps
  actualOppositeCapacity : goodCusps ≤ 2 * noNeighbors

theorem stronger_carrier_bound (d : DefectIdeal5Data) :
    1 ≤ d.noNeighbors ∧
      3 * d.base.threeFacetVertices + 3 ≤ d.base.carrierVertices := by
  have defect := positive_defect_count d.base.threeFacetVertices d.exceptionalEdges
    d.smallCorners d.goodCusps d.noNeighbors d.polygonCornerLower
    d.pairedCornerUpper d.actualOppositeCapacity
  have vertices := d.verticesPartition
  have outside := d.outsidePartition
  have outgoing := d.outgoingEdges
  have triangles := d.triangularBijection
  omega

theorem ideal5_counts (d : DefectIdeal5Data) :
    62 ≤ d.base.facetCounts.length ∧
      3 * d.base.facetCounts.length + 3 ≤ d.base.vertices ∧
      189 ≤ d.base.vertices := by
  have carrier := stronger_carrier_bound d
  have oldCounts := RightAngledCuspCapacity.ideal5_counts d.base
  have carrierFacets := (carrier_bounds d.base).2
  have large := d.base.largeThreeFacet
  have injection := d.base.actualFacetInjection
  omega

theorem sum_ideal5_facet_minima (xs : List DefectIdeal5Data) :
    62 * xs.length ≤ (xs.map (fun d => d.base.facetCounts.length)).sum := by
  induction xs with
  | nil => simp
  | cons a xs ih =>
    have bound := (ideal5_counts a).1
    simp only [List.length_cons, List.map_cons, List.sum_cons]
    omega

theorem ideal6_counts
    (fiveFacets : List DefectIdeal5Data) (vertices fourFaces : Nat)
    (actualFacetInjection : ∃ d ∈ fiveFacets,
      1 + d.base.facetCounts.length + d.base.vertices ≤ fiveFacets.length)
    (cofaceIncidence :
      (fiveFacets.map (fun d => d.base.facetCounts.length)).sum = 2 * fourFaces)
    (reciprocalTransport :
      99 * fourFaces ≤ 99 * fiveFacets.length + 136 * vertices) :
    252 ≤ fiveFacets.length ∧
      1485 * fiveFacets.length ≤ 68 * vertices ∧ 5504 ≤ vertices := by
  obtain ⟨d, _, injection⟩ := actualFacetInjection
  have sizes := ideal5_counts d
  have cofaceLower := sum_ideal5_facet_minima fiveFacets
  omega

theorem exact_ideal6_rounding :
    (1485 * 252 : Nat) = 374220 ∧
      68 * 5503 < 374220 ∧ 374220 ≤ 68 * 5504 := by decide

theorem ideal5_density11_12 (faces ideal : Nat)
    (facetAverage : 16 * faces + 46 * ideal < 20 * faces) :
    23 * ideal < 2 * faces := by omega

theorem ideal6_density11_12 (faces ideal : Nat)
    (facetAverage : 27 * faces + 225 * ideal < 42 * faces) :
    15 * ideal < faces := by omega

end CoxeterProofs.RightAngledCuspDefect

#print axioms CoxeterProofs.RightAngledCuspDefect.positive_defect_count
#print axioms CoxeterProofs.RightAngledCuspDefect.stronger_carrier_bound
#print axioms CoxeterProofs.RightAngledCuspDefect.ideal5_counts
#print axioms CoxeterProofs.RightAngledCuspDefect.ideal6_counts
#print axioms CoxeterProofs.RightAngledCuspDefect.exact_ideal6_rounding
#print axioms CoxeterProofs.RightAngledCuspDefect.ideal5_density11_12
#print axioms CoxeterProofs.RightAngledCuspDefect.ideal6_density11_12
