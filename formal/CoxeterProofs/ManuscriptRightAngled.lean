import CoxeterProofs.FiveFaceCarrier
import CoxeterProofs.RightAngledDensity
import CoxeterProofs.ManuscriptDataBounds

/-! Conditional right-angled necessary conditions. The real Lorentz,
reciprocal-row and cusp incidence arguments remain explicit interfaces;
these arithmetic statements do not exclude dimensions12 or6. -/
namespace CoxeterProofs.ManuscriptRightAngled

open ManuscriptDataBounds

/-- No published result or script is assumed by this finite-count lemma. -/
theorem chain_kappa_nonpositive (m kappa : Int) (positive : 1 ≤ m)
    (chainCounts : kappa = 3 * (18 * m + 6) - 4 * (14 * m + 10) + 24) :
    kappa ≤ 0 := by omega

/-- The identity is the ideal5 box/Euler incidence interface, not a
numerical observation about one construction. The witness is an actual
listed four-facet once the finite list is instantiated geometrically. -/
theorem ideal5_has_nonchain_kappa (kappas : List Int)
    (boxEulerIdentity : kappas.sum = 48) :
    ∃ k ∈ kappas, ¬ ∃ m : Int, 1 ≤ m ∧
      k = 3 * (18 * m + 6) - 4 * (14 * m + 10) + 24 := by
  obtain ⟨k, hk, positive⟩ := FiveFaceCarrier.positive_carrier kappas boxEulerIdentity
  refine ⟨k, hk, ?_⟩
  intro witness
  obtain ⟨m, hm, counts⟩ := witness
  have nonpositive := chain_kappa_nonpositive m k hm counts
  omega

/-- The additional geometric bridge from the Oct/D/A4 three-face spectrum
to chain counts is explicit. It is the new written exclusion/decomposition
theorem, not a published assertion. -/
theorem ideal5_excludes_small_three_face_spectrum
    {Facet : Type} (facets : List Facet) (kappa : Facet → Int)
    (SmallSpectrum : Prop)
    (boxEulerIdentity : (facets.map kappa).sum = 48)
    (writtenChainBridge : SmallSpectrum → ∀ q ∈ facets,
      ∃ m : Int, 1 ≤ m ∧
        kappa q = 3 * (18 * m + 6) - 4 * (14 * m + 10) + 24) :
    ¬ SmallSpectrum := by
  intro spectrum
  obtain ⟨k, hk, nonchain⟩ :=
    ideal5_has_nonchain_kappa (facets.map kappa) boxEulerIdentity
  obtain ⟨q, hq, value⟩ := List.mem_map.mp hk
  apply nonchain
  rw [← value]
  exact writtenChainBridge spectrum q hq

/-- Analytic transport has already supplied this integer count comparison:
small deficit1/8 and exceptional allowance48. It is not claimed that the
actual reciprocal rows themselves are integer-valued. -/
theorem carrier_excess_fraction_from_transport (finite bad excess : Int)
    (transportComparison : finite - bad ≤ 384 * bad)
    (largeCarrierCover : bad ≤ 14 * excess) :
    finite ≤ 5390 * excess := by omega

/-- Dimension12 necessary condition. The Schur premise is the exact
integer clearing of 2e+boundary/4≤(516992/3125)*vertices. Its geometric
row proof and the true-twin capacity are not proved by this arithmetic. -/
theorem nine_part_boundary (vertices internal boundary : Int)
    (degree : 190 * vertices = 2 * internal + boundary)
    (schur : 25000 * internal + 3125 * boundary ≤ 2067968 * vertices) :
    102344 * vertices ≤ 3125 * boundary := by omega

theorem nine_part_further_support
    (vertices internal boundary complete further : Int)
    (degree : 190 * vertices = 2 * internal + boundary)
    (schur : 25000 * internal + 3125 * boundary ≤ 2067968 * vertices)
    (partition : boundary = complete + further)
    (trueTwins : complete ≤ 20 * vertices) :
    39844 * vertices ≤ 3125 * further := by
  have boundaryLower := nine_part_boundary vertices internal boundary degree schur
  omega

/-- The inherited ideal6 size consequence retains its published/geometry
parameters instead of asserting a universal dimension exclusion. -/
theorem ideal6_necessary_size
    (facetFacets facetVertices facets vertices : Int)
    (facetSize : 56 ≤ facetFacets) (vertexSize : 170 ≤ facetVertices)
    (actualInjection : 1 + facetFacets + facetVertices ≤ facets)
    (boxIncidence : 170 * facets ≤ 10 * vertices) :
    227 ≤ facets ∧ 3859 ≤ vertices :=
  RightAngledDensity.ideal6_size facetFacets facetVertices facets vertices
    facetSize vertexSize actualInjection boxIncidence

/-- New sharp reciprocal/cusp-incidence consequence. The premise is the
cleared geometric inequality, not a sampled estimate or a formalization
of the actual reciprocal transport. Integer rounding is proved here. -/
theorem ideal6_sharp_vertex_lower4462 (facets vertices : Int)
    (facetMinimum : 227 ≤ facets)
    (reciprocalIncidence : 2673 * facets ≤ 136 * vertices) :
    4462 ≤ vertices := by omega

theorem ideal6_sharp_rounding_exact :
    (2673 * 227 : Int) = 606771 ∧
      136 * 4461 < 606771 ∧ 606771 ≤ 136 * 4462 := by decide

end CoxeterProofs.ManuscriptRightAngled

#print axioms CoxeterProofs.ManuscriptRightAngled.ideal5_excludes_small_three_face_spectrum
#print axioms CoxeterProofs.ManuscriptRightAngled.carrier_excess_fraction_from_transport
#print axioms CoxeterProofs.ManuscriptRightAngled.nine_part_further_support
#print axioms CoxeterProofs.ManuscriptRightAngled.ideal6_necessary_size
#print axioms CoxeterProofs.ManuscriptRightAngled.ideal6_sharp_vertex_lower4462
#print axioms CoxeterProofs.ManuscriptRightAngled.ideal6_sharp_rounding_exact
