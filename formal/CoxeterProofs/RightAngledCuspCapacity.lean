import Std

/-!
Conditional integer logic for the reviewed outside-cusp capacity theorem.
The geometric capacity-two lemma, Figure7/small-facet bridge, actual
facet injection, published minimum24, box/Euler incidences, and reciprocal
transport remain explicit inputs. No hyperbolic realization or dimension
exclusion is asserted. No final size bound is an assumed axiom.
-/
namespace CoxeterProofs.RightAngledCuspCapacity

/-- Exact deficit in an ideal four-polytope. The three classes count outside
vertices having zero, one, or two neighbors in the fixed facet. -/
theorem dimension4_capacity_deficit
    (vertices facetVertices outside zeroCount one two : Nat)
    (partition : vertices = facetVertices + outside)
    (outsidePartition : outside = zeroCount + one + two)
    (edgeCount : 4 * facetVertices = one + 2 * two) :
    2 * vertices = 6 * facetVertices + 2 * zeroCount + one ∧
      3 * facetVertices ≤ vertices := by omega

theorem sum_minimum24 (xs : List Nat)
    (minimum : ∀ a ∈ xs, 24 ≤ a) : 24 * xs.length ≤ xs.sum := by
  induction xs with
  | nil => simp
  | cons a xs ih =>
    have head := minimum a (by simp)
    have tail := ih (fun b hb => minimum b (by simp [hb]))
    simp only [List.length_cons, List.sum_cons]
    omega

theorem sum_strict_minimum24 (xs : List Nat)
    (minimum : ∀ a ∈ xs, 24 ≤ a)
    (large : ∃ a ∈ xs, 25 ≤ a) : 24 * xs.length < xs.sum := by
  induction xs with
  | nil => simp at large
  | cons a xs ih =>
    have head := minimum a (by simp)
    have tailMinimum : ∀ b ∈ xs, 24 ≤ b :=
      fun b hb => minimum b (by simp [hb])
    obtain ⟨b, hb, big⟩ := large
    rcases List.mem_cons.mp hb with same | rest
    · have tail := sum_minimum24 xs tailMinimum
      simp only [List.length_cons, List.sum_cons]
      omega
    · have tail := ih tailMinimum ⟨b, rest, big⟩
      simp only [List.length_cons, List.sum_cons]
      omega

/-- Geometric interfaces for an actual ideal five-polytope.
`threeFacetVertices` belongs to the positive-kappa carrier. Its lower11
comes from the reviewed restricted Figure7 theorem plus the complete
published small-three-polyhedron classification; it is not a census of
four- or five-polytopes. `facetCounts` lists the facet numbers of ALL
actual four-facets, including the selected carrier. -/
structure Ideal5Data where
  vertices : Nat
  threeFaces : Nat
  facetCounts : List Nat
  threeFacetVertices : Nat
  carrierVertices : Nat
  carrierFacets : Nat
  largeThreeFacet : 11 ≤ threeFacetVertices
  outsideCapacity : 3 * threeFacetVertices ≤ carrierVertices
  carrierInjection : 2 * threeFacetVertices + 3 ≤ carrierFacets
  actualFacetInjection :
    1 + carrierFacets + carrierVertices ≤ facetCounts.length
  actualCarrier : carrierFacets ∈ facetCounts
  publishedFacetMinimum : ∀ a ∈ facetCounts, 24 ≤ a
  cofaceIncidence : facetCounts.sum = 2 * threeFaces
  boxEuler : 3 * vertices + 3 * facetCounts.length = threeFaces + 6

theorem carrier_bounds (d : Ideal5Data) :
    33 ≤ d.carrierVertices ∧ 25 ≤ d.carrierFacets := by
  have h := d.largeThreeFacet
  have capacity := d.outsideCapacity
  have injection := d.carrierInjection
  omega

theorem ideal5_counts (d : Ideal5Data) :
    59 ≤ d.facetCounts.length ∧
      3 * d.facetCounts.length + 3 ≤ d.vertices ∧ 180 ≤ d.vertices := by
  obtain ⟨carrierVertices, carrierFacets⟩ := carrier_bounds d
  have injection := d.actualFacetInjection
  have strictSum := sum_strict_minimum24 d.facetCounts
    d.publishedFacetMinimum ⟨d.carrierFacets, d.actualCarrier, carrierFacets⟩
  have incidence := d.cofaceIncidence
  have euler := d.boxEuler
  omega

theorem sum_ideal5_facet_minima (xs : List Ideal5Data) :
    59 * xs.length ≤ (xs.map (fun d => d.facetCounts.length)).sum := by
  induction xs with
  | nil => simp
  | cons a xs ih =>
    have bound := (ideal5_counts a).1
    simp only [List.length_cons, List.map_cons, List.sum_cons]
    omega

/-- The universal reciprocal transport premise is cleared over integers,
before any facet-minimum substitution. Every five-facet has its own
explicit ideal5 interface, so its new59 minimum is derived above. -/
theorem ideal6_counts
    (fiveFacets : List Ideal5Data) (vertices fourFaces : Nat)
    (actualFacetInjection : ∃ d ∈ fiveFacets,
      1 + d.facetCounts.length + d.vertices ≤ fiveFacets.length)
    (cofaceIncidence :
      (fiveFacets.map (fun d => d.facetCounts.length)).sum = 2 * fourFaces)
    (reciprocalTransport :
      99 * fourFaces ≤ 99 * fiveFacets.length + 136 * vertices) :
    240 ≤ fiveFacets.length ∧
      5643 * fiveFacets.length ≤ 272 * vertices ∧ 4980 ≤ vertices := by
  obtain ⟨d, _, injection⟩ := actualFacetInjection
  have sizes := ideal5_counts d
  have cofaceLower := sum_ideal5_facet_minima fiveFacets
  omega

theorem exact_ideal6_rounding :
    (5643 * 240 : Nat) = 1354320 ∧
      272 * 4979 < 1354320 ∧ 1354320 ≤ 272 * 4980 := by decide

/-- NK averages and the ideal minima enter through explicit cleared
incidence premises. These strict density inequalities do not exclude
ideal dimension6 or mixed dimension12. -/
theorem ideal5_density11_12 (faces ideal : Nat)
    (facetAverage : 16 * faces + 43 * ideal < 20 * faces) :
    43 * ideal < 4 * faces := by omega

theorem ideal5_density9_10 (faces ideal : Nat)
    (facetAverage : 16 * faces + 43 * ideal < 30 * faces) :
    43 * ideal < 14 * faces := by omega

theorem ideal6_density11_12 (faces ideal : Nat)
    (facetAverage : 27 * faces + 213 * ideal < 42 * faces) :
    71 * ideal < 5 * faces := by omega

theorem mixed_vertex_ratio11 (faces ideal finiteVertices idealVertices : Nat)
    (density : 43 * ideal < 4 * faces)
    (finiteCover : faces ≤ ideal + 462 * finiteVertices)
    (edgeAverage : 77 * finiteVertices + 7168 * idealVertices < 22 * faces) :
    13312 * idealVertices < 20669 * finiteVertices := by omega

theorem mixed_vertex_ratio12 (faces ideal finiteVertices idealVertices : Nat)
    (density : 43 * ideal < 4 * faces)
    (finiteCover : faces ≤ ideal + 792 * finiteVertices)
    (edgeAverage : 6 * finiteVertices + 1024 * idealVertices < faces) :
    6656 * idealVertices < 5637 * finiteVertices := by omega

end CoxeterProofs.RightAngledCuspCapacity

#print axioms CoxeterProofs.RightAngledCuspCapacity.dimension4_capacity_deficit
#print axioms CoxeterProofs.RightAngledCuspCapacity.sum_strict_minimum24
#print axioms CoxeterProofs.RightAngledCuspCapacity.ideal5_counts
#print axioms CoxeterProofs.RightAngledCuspCapacity.ideal6_counts
#print axioms CoxeterProofs.RightAngledCuspCapacity.exact_ideal6_rounding
#print axioms CoxeterProofs.RightAngledCuspCapacity.ideal5_density11_12
#print axioms CoxeterProofs.RightAngledCuspCapacity.ideal5_density9_10
#print axioms CoxeterProofs.RightAngledCuspCapacity.ideal6_density11_12
#print axioms CoxeterProofs.RightAngledCuspCapacity.mixed_vertex_ratio11
#print axioms CoxeterProofs.RightAngledCuspCapacity.mixed_vertex_ratio12
