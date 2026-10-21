import CoxeterProofs.RightAngledCuspDefect

/-!
Conditional logic for the reviewed ideal5 facet propagation theorem.
The exact outside/triangle identity, spinor corner counts and actual
good-opposite incidence remain explicit geometric inputs. The final
206/269/5875 bounds are derived, not assumed.
-/
namespace CoxeterProofs.RightAngledCuspPropagation

open RightAngledCuspCapacity RightAngledCuspDefect

theorem facet_vertex_propagation
    (vertices facetVertices facetFacets noNeighbors exceptional small good : Nat)
    (outsideIdentity : vertices = 5 * facetVertices + noNeighbors + exceptional)
    (polygonCorners : 6 * facetVertices + 6 * facetFacets ≤ small + 3 * exceptional)
    (pairedCorners : small ≤ 8 * facetVertices + 4 * good)
    (actualOppositeIncidence : 9 * good ≤ 2 * noNeighbors) :
    13 * facetVertices + 6 * facetFacets ≤ 3 * vertices := by omega

/-- The selected positive carrier is an actual four-facet. Its older
defect interface supplies36/25, while these new fields describe its
embedding in the actual ideal five-polytope. -/
structure PropagatingIdeal5Data where
  base : DefectIdeal5Data
  noNeighbors : Nat
  exceptionalEdges : Nat
  smallCorners : Nat
  goodCusps : Nat
  outsideIdentity : base.base.vertices =
    5 * base.base.carrierVertices + noNeighbors + exceptionalEdges
  polygonCornerLower :
    6 * base.base.carrierVertices + 6 * base.base.carrierFacets ≤
      smallCorners + 3 * exceptionalEdges
  pairedCornerUpper : smallCorners ≤ 8 * base.base.carrierVertices + 4 * goodCusps
  actualOppositeIncidence : 9 * goodCusps ≤ 2 * noNeighbors

theorem ideal5_counts (d : PropagatingIdeal5Data) :
    62 ≤ d.base.base.facetCounts.length ∧ 206 ≤ d.base.base.vertices := by
  have oldCounts := RightAngledCuspDefect.ideal5_counts d.base
  have carrier := stronger_carrier_bound d.base
  have large := d.base.base.largeThreeFacet
  have carrierFacets := (carrier_bounds d.base.base).2
  have propagation := facet_vertex_propagation d.base.base.vertices
    d.base.base.carrierVertices d.base.base.carrierFacets d.noNeighbors
    d.exceptionalEdges d.smallCorners d.goodCusps d.outsideIdentity
    d.polygonCornerLower d.pairedCornerUpper d.actualOppositeIncidence
  omega

theorem sum_ideal5_facet_minima (xs : List PropagatingIdeal5Data) :
    62 * xs.length ≤ (xs.map (fun d => d.base.base.facetCounts.length)).sum := by
  induction xs with
  | nil => simp
  | cons a xs ih =>
    have bound := (ideal5_counts a).1
    simp only [List.length_cons, List.map_cons, List.sum_cons]
    omega

theorem ideal6_counts
    (fiveFacets : List PropagatingIdeal5Data) (vertices fourFaces : Nat)
    (actualFacetInjection : ∃ d ∈ fiveFacets,
      1 + d.base.base.facetCounts.length + d.base.base.vertices ≤ fiveFacets.length)
    (cofaceIncidence :
      (fiveFacets.map (fun d => d.base.base.facetCounts.length)).sum = 2 * fourFaces)
    (reciprocalTransport :
      99 * fourFaces ≤ 99 * fiveFacets.length + 136 * vertices) :
    269 ≤ fiveFacets.length ∧
      1485 * fiveFacets.length ≤ 68 * vertices ∧ 5875 ≤ vertices := by
  obtain ⟨d, _, injection⟩ := actualFacetInjection
  have sizes := ideal5_counts d
  have cofaceLower := sum_ideal5_facet_minima fiveFacets
  omega

theorem exact_ideal6_rounding :
    (1485 * 269 : Nat) = 399465 ∧
      68 * 5874 < 399465 ∧ 399465 ≤ 68 * 5875 := by decide

theorem ideal6_density11_12 (faces ideal : Nat)
    (facetAverage : 27 * faces + 242 * ideal < 42 * faces) :
    242 * ideal < 15 * faces := by omega

end CoxeterProofs.RightAngledCuspPropagation

#print axioms CoxeterProofs.RightAngledCuspPropagation.facet_vertex_propagation
#print axioms CoxeterProofs.RightAngledCuspPropagation.ideal5_counts
#print axioms CoxeterProofs.RightAngledCuspPropagation.ideal6_counts
#print axioms CoxeterProofs.RightAngledCuspPropagation.exact_ideal6_rounding
#print axioms CoxeterProofs.RightAngledCuspPropagation.ideal6_density11_12
