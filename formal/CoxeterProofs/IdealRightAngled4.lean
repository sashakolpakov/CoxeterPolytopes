import Std

/-!
Conditional logic for the maximum-cusp-height ideal right-angled bound.

This module does not formalize Lorentz space, real heights, spherical codes,
Musin's theorem, Kolpakov's theorem, face inheritance, or finite graph
connectedness for convex polytopes. These are separate explicit inputs.
In particular `Height` is an arbitrary type, not an integer substitute for
real height. Geometry supplies its actual maximum and the local normalized
height bounds. The kernel checks the cardinality squeeze, the consequence
of metric equality, propagation along finite paths, and the contradiction.

Written proof: ideal_right_angled_maximum_cusp_height_kissing_candidate_
20261020.md, independently reviewed in review_ideal_right_angled_maximum_
height_kissing_batch_guard_20261020.md. In the journal version Kolpakov's
minimum-facet/metric-equality result is Theorem 4, not Theorem 2.

SOURCE PREPARED ONLY. Compilation and actual axiom audit are scheduled
separately by the coordinating agent under the serial resource gates.
-/

namespace CoxeterProofs.IdealRightAngled4

universe u v

/-- A finite sequence of actual ridge-adjacency steps. -/
inductive FacetPath {F : Type u} (adjacent : F → F → Prop) : F → F → Prop where
  | nil (a : F) : FacetPath adjacent a a
  | step {a b c : F} :
      adjacent a b → FacetPath adjacent b c → FacetPath adjacent a c

theorem property_along_path {F : Type u} {adjacent : F → F → Prop}
    (good : F → Prop)
    (closed : ∀ a b, adjacent a b → good a → good b)
    {a b : F} (path : FacetPath adjacent a b) : good a → good b := by
  induction path with
  | nil _ => exact fun h => h
  | step edge _ ih =>
      intro h
      exact ih (closed _ _ edge h)

/-- Finite actual-facet data for a hypothetical ideal right-angled H5
polytope. The abstract predicates stand for the named real geometric
objects; they are not axioms or conclusions about existence of a polytope.

`actualFoot f`: the explicitly constructed future unit point lies in facet f.
`boundedHeights f`: all intrinsic neighboring temporal coefficients lie in [0,1].
`radialCode f`: the distinct actual neighbors form an injective 60-degree S3 code.
`idealFour f`: the facet is an ideal finite-volume right-angled H4 polytope.
`regular24 f`: it is metrically the regular ideal right-angled 24-cell.
`centered f`: the actual foot equals its regular-cell center.

`maximumHeight` and its nonzero value come from the finite positive maximum
in the written geometry, not an assertion that a final dimension bound holds.
-/
structure FiveData (facetCount : Nat) (Height : Type v) where
  height : Fin facetCount → Height
  zeroHeight : Height
  maximumHeight : Height
  maximumNonzero : maximumHeight ≠ zeroHeight
  maximalFacet : Fin facetCount
  maximalFacetHeight : height maximalFacet = maximumHeight
  cuspFacet : Fin facetCount
  cuspFacetHeight : height cuspFacet = zeroHeight
  adjacent : Fin facetCount → Fin facetCount → Prop
  neighbors : Fin facetCount → List (Fin facetCount)
  neighborsNodup : ∀ f, (neighbors f).Nodup
  neighborsExact : ∀ f g, g ∈ neighbors f ↔ adjacent f g
  connected : ∀ f g, FacetPath adjacent f g
  actualFoot : Fin facetCount → Prop
  boundedHeights : Fin facetCount → Prop
  radialCode : Fin facetCount → Prop
  idealFour : Fin facetCount → Prop
  regular24 : Fin facetCount → Prop
  centered : Fin facetCount → Prop

/-- External theorems, explicitly separated from written chamber geometry.
Musin supplies only cardinality. Kolpakov supplies the minimum and METRIC
equality, so uniqueness of a kissing configuration is not an assumption. -/
structure PublishedInputs {n : Nat} {Height : Type v}
    (d : FiveData n Height) : Prop where
  musinCardinality : ∀ f, d.radialCode f → (d.neighbors f).length ≤ 24
  kolpakovMinimum : ∀ f, d.idealFour f → 24 ≤ (d.neighbors f).length
  kolpakovMetricEquality :
    ∀ f, d.idealFour f → (d.neighbors f).length = 24 → d.regular24 f

/-- Written real geometry, not yet formalized here.
The foot uses all-wall acuteness; the code uses intrinsic unit normals.
The last two fields are respectively the explicit opposite-normal-pair
center proof and its exact height-ratio consequence at each adjacent facet.
No global propagation or final nonexistence is assumed in this structure. -/
structure WrittenGeometry {n : Nat} {Height : Type v}
    (d : FiveData n Height) : Prop where
  idealFacetInheritance : ∀ f, d.idealFour f
  maximumFoot :
    ∀ f, d.height f = d.maximumHeight → d.actualFoot f
  maximumNormalizedBounds :
    ∀ f, d.height f = d.maximumHeight → d.actualFoot f → d.boundedHeights f
  radialCodeReduction :
    ∀ f, d.actualFoot f → d.boundedHeights f → d.radialCode f
  oppositePairsForceCenter :
    ∀ f, d.regular24 f → d.actualFoot f → d.boundedHeights f → d.centered f
  centeredNeighborHeight :
    ∀ f, d.height f = d.maximumHeight → d.centered f →
      ∀ g, d.adjacent f g → d.height g = d.maximumHeight

theorem maximal_local_data {n : Nat} {Height : Type v}
    (d : FiveData n Height) (geometry : WrittenGeometry d)
    (f : Fin n) (maximum : d.height f = d.maximumHeight) :
    d.actualFoot f ∧ d.boundedHeights f ∧ d.radialCode f := by
  have foot := geometry.maximumFoot f maximum
  have bounds := geometry.maximumNormalizedBounds f maximum foot
  exact ⟨foot, bounds, geometry.radialCodeReduction f foot bounds⟩

theorem maximal_neighbor_count_eq24 {n : Nat} {Height : Type v}
    (d : FiveData n Height) (published : PublishedInputs d)
    (geometry : WrittenGeometry d) (f : Fin n)
    (maximum : d.height f = d.maximumHeight) :
    (d.neighbors f).length = 24 := by
  obtain ⟨_, _, code⟩ := maximal_local_data d geometry f maximum
  have upper := published.musinCardinality f code
  have lower := published.kolpakovMinimum f (geometry.idealFacetInheritance f)
  omega

theorem maximal_facet_regular24 {n : Nat} {Height : Type v}
    (d : FiveData n Height) (published : PublishedInputs d)
    (geometry : WrittenGeometry d) (f : Fin n)
    (maximum : d.height f = d.maximumHeight) : d.regular24 f := by
  exact published.kolpakovMetricEquality f (geometry.idealFacetInheritance f)
    (maximal_neighbor_count_eq24 d published geometry f maximum)

theorem maximal_neighbor_height {n : Nat} {Height : Type v}
    (d : FiveData n Height) (published : PublishedInputs d)
    (geometry : WrittenGeometry d) (f g : Fin n)
    (maximum : d.height f = d.maximumHeight) (edge : d.adjacent f g) :
    d.height g = d.maximumHeight := by
  obtain ⟨foot, bounds, _⟩ := maximal_local_data d geometry f maximum
  have regular := maximal_facet_regular24 d published geometry f maximum
  have center := geometry.oppositePairsForceCenter f regular foot bounds
  exact geometry.centeredNeighborHeight f maximum center g edge

theorem all_facets_maximal {n : Nat} {Height : Type v}
    (d : FiveData n Height) (published : PublishedInputs d)
    (geometry : WrittenGeometry d) (f : Fin n) :
    d.height f = d.maximumHeight := by
  have closed : ∀ a b, d.adjacent a b →
      d.height a = d.maximumHeight → d.height b = d.maximumHeight := by
    intro a b edge maximum
    exact maximal_neighbor_height d published geometry a b maximum edge
  exact property_along_path (fun a => d.height a = d.maximumHeight) closed
    (d.connected d.maximalFacet f) d.maximalFacetHeight

theorem no_ideal_five_data {n : Nat} {Height : Type v}
    (d : FiveData n Height) (published : PublishedInputs d)
    (geometry : WrittenGeometry d) : False := by
  have maximal := all_facets_maximal d published geometry d.cuspFacet
  exact d.maximumNonzero (maximal.symm.trans d.cuspFacetHeight)

/-- All fields of a proposed five-face instance, including the explicit
external theorem hypotheses. This is data, not an axiom that an instance
exists. The preceding theorem proves that no instance can exist. -/
structure FiveFaceCertificate (n : Nat) (Height : Type v) where
  data : FiveData n Height
  published : PublishedInputs data
  geometry : WrittenGeometry data

theorem no_five_face_certificate {n : Nat} {Height : Type v}
    (certificate : FiveFaceCertificate n Height) : False :=
  no_ideal_five_data certificate.data certificate.published certificate.geometry

/-- Actual face inheritance supplies a five-dimensional certificate from
any ideal finite-volume right-angled polytope of dimension at least five.
Its real geometric construction remains an explicit hypothesis here. -/
theorem dimension_at_most_four {Height : Type v} (dimension : Nat)
    (fiveFaceInheritance : 5 ≤ dimension →
      ∃ facetCount : Nat, Nonempty (FiveFaceCertificate facetCount Height)) :
    dimension ≤ 4 := by
  have excluded : ¬ 5 ≤ dimension := by
    intro large
    obtain ⟨_, ⟨certificate⟩⟩ := fiveFaceInheritance large
    exact no_five_face_certificate certificate
  omega

end CoxeterProofs.IdealRightAngled4

#print axioms CoxeterProofs.IdealRightAngled4.property_along_path
#print axioms CoxeterProofs.IdealRightAngled4.maximal_neighbor_count_eq24
#print axioms CoxeterProofs.IdealRightAngled4.maximal_facet_regular24
#print axioms CoxeterProofs.IdealRightAngled4.maximal_neighbor_height
#print axioms CoxeterProofs.IdealRightAngled4.all_facets_maximal
#print axioms CoxeterProofs.IdealRightAngled4.no_ideal_five_data
#print axioms CoxeterProofs.IdealRightAngled4.no_five_face_certificate
#print axioms CoxeterProofs.IdealRightAngled4.dimension_at_most_four
