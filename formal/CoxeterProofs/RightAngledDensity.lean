import Std

/-! Linear consequences of the written carrier, face-injection and published
average bounds. Geometric incidence assertions are explicit hypotheses. -/
namespace CoxeterProofs.RightAngledDensity

theorem carrier_vertices (facets vertices : Int)
    (facetMinimum : 25 ≤ facets) (defect : 5 ≤ vertices-facets) :
    30 ≤ vertices := by omega

theorem ideal5_facets (facets vertices ambientFacets : Int)
    (facetMinimum : 25 ≤ facets) (defect : 5 ≤ vertices-facets)
    (injection : 1+facets+vertices ≤ ambientFacets) :
    56 ≤ ambientFacets := by omega

theorem ideal5_vertices (vertices facets threeFaces : Int)
    (euler : threeFaces = 3*vertices+3*facets-6)
    (facetIncidence : 24*facets ≤ 2*threeFaces)
    (facetMinimum : 56 ≤ facets) :
    3*facets+2 ≤ vertices ∧ 170 ≤ vertices := by omega

theorem ideal6_size (facetFacets facetVertices facets vertices : Int)
    (facetSize : 56 ≤ facetFacets) (vertexSize : 170 ≤ facetVertices)
    (injection : 1+facetFacets+facetVertices ≤ facets)
    (incidence : 170*facets ≤ 10*vertices) :
    227 ≤ facets ∧ 3859 ≤ vertices := by omega

theorem ideal5_density_9_10 (faces ideal : Int)
    (average : 16*faces+40*ideal < 30*faces) :
    20*ideal < 7*faces := by omega

theorem ideal5_density_11_12 (faces ideal : Int)
    (average : 16*faces+40*ideal < 20*faces) :
    10*ideal < faces := by omega

theorem ideal6_density_11_12 (faces ideal : Int)
    (average : 27*faces+200*ideal < 42*faces) :
    40*ideal < 3*faces := by omega

theorem vertex_ratio9 (faces ideal finiteVertices idealVertices : Int)
    (average : 16*faces+40*ideal < 30*faces)
    (cover : faces-ideal ≤ 126*finiteVertices)
    (edgeAverage : 9*finiteVertices+256*idealVertices < 18*faces) :
    3328*idealVertices < 45243*finiteVertices := by omega

theorem vertex_ratio10 (faces ideal finiteVertices idealVertices : Int)
    (average : 16*faces+40*ideal < 30*faces)
    (cover : faces-ideal ≤ 252*finiteVertices)
    (edgeAverage : 5*finiteVertices+256*idealVertices < 5*faces) :
    3328*idealVertices < 25135*finiteVertices := by omega

theorem vertex_ratio11 (faces ideal finiteVertices idealVertices : Int)
    (average : 16*faces+40*ideal < 20*faces)
    (cover : faces-ideal ≤ 462*finiteVertices)
    (edgeAverage : 77*finiteVertices+7168*idealVertices < 22*faces) :
    3072*idealVertices < 4807*finiteVertices := by omega

theorem vertex_ratio12 (faces ideal finiteVertices idealVertices : Int)
    (average : 16*faces+40*ideal < 20*faces)
    (cover : faces-ideal ≤ 792*finiteVertices)
    (edgeAverage : 6*finiteVertices+1024*idealVertices < faces) :
    512*idealVertices < 437*finiteVertices := by omega

end CoxeterProofs.RightAngledDensity
