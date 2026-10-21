# Intrinsic golden cusp construction and higher-base compact search

No new finite-volume or compact hyperbolic Coxeter polytope, and no new dimension, is established. The main positive result is an intrinsic **587-facet H20 chamber**, with exact Coxeter angles and 252 certified finite vertices, accompanied by an exact proof that its volume is infinite. New constructions deliberately leave the previously excluded rational ambient-cap route.

The [dimension restrictions](arithmetic_dimension_bounds_and_effective_high_dimension_search.md) remain decisive: noncompact arithmetic examples have dimension at most 21 and cannot have dimension 20. Compact arithmetic examples are not subject to the 21 bound; the general compact bound is 29. The working targets are intrinsic nonarithmetic H20 constructions and compact H9 incidence constructions. Allcock's published edge-walking algorithm operates over integer lattices; the present geometric adaptations require their own exact metric and global-volume certificates.

## Intrinsic H20 construction

The retained cusp diagram is affine **A11 × A2 × E6**, with 22 facets and a 19-dimensional Euclidean link. The E6 factor corrects an earlier descriptive D6 name; the old literal Gram and rank calculations are unchanged. Set its three factor inradii to `(1,phi,1)`, with `phi²=phi+1`.

The [complete endpoint profile census](fresh_h20_intrinsic_caps_20261014_endpoint.md) allows arbitrary real new normals. Its 260 finite-pool profiles comprise 144 finite and 116 semidefinite cases, expanding to 16,272 literal occurrences, together with all isolated dihedral tails. H3 and H4 attachments are included. The [separate review](fresh_h20_intrinsic_caps_20261014_review.md) independently checks all accepted records and transports, and source-reviews the exhaustion argument.

The construction proceeds through three frozen exact checkpoints:

| Chamber | Constructive result | Global result |
| --- | --- | --- |
| [275 facets](fresh_h20_intrinsic_caps_20261014_golden_seed.md) | Central radius-2 cap and 252 radius-1 vertex caps; all original cusp edges end finitely | Strict half-vertex escape |
| [527 facets](fresh_h20_intrinsic_caps_20261014_golden_extension.md) | 252 radius-1 half-vertex caps; all old facets and certified finite vertices survive | Strict quarter-vertex escape |
| [587 facets](fresh_h20_intrinsic_caps_20261014_quarter_face_exact.md) | 60 jointly compatible radius-3/2 face-center caps cover all 252 quarter probes | Strict eighth-vertex escape of square 1/2 |
| [587 facets after replacement](fresh_h20_intrinsic_caps_20261014_radius4_replacement.md) | A radius-4 endpoint cap replaces one radius-1 cap and cuts the annular escape; 251 old finite vertices survive and one new finite endpoint is certified | Eighth-vertex escape survives |

Every inequality in these constructions defines an actual facet. The [independent FLINT matrix replay](fresh_h20_intrinsic_caps_20261014_matrix_review.md) checks all 171,991 pairs in the final chamber, its common strict timelike center, all 252 finite endpoints and their positive-definite rank-20 incident Grams, and witnesses for all 587 actual facets. It also verifies the strict escape. It does not claim a complete corner census. A [separate symbolic source review](fresh_h20_intrinsic_caps_20261014_golden_review.md) checks the geometric formulas of the first two layers.

The [independent replacement review](fresh_h20_intrinsic_caps_20261014_replacement_review.md) subsequently checks all 586 changed pairs and the new endpoint of norm `-8/3`. It inherits unchanged geometry from the full matrix replay. This replacement is a concrete step outside the small-radius obstruction; it is still not a finite-volume example.

Its new parallel pair also has an [exact open-cusp diagnostic](fresh_h20_intrinsic_caps_20261014_parallel_cusp_gap.md): a feasible null ray has 17 incident walls of Gram rank 16, leaving three unbounded horospherical directions. This is an additional explicit unfinished end, not a bounded cusp.

The golden central wall has label-5 contacts. Its irrational Cartan cycle and rational indefinite conjugate block exclude quasi-arithmeticity of any hypothetical cofinite completion retaining this core, as detailed in the 275-facet report. This is conditional, not an existence result.

## Exact restrictions that redirect the search

The [all-radius endpoint bound](fresh_h20_intrinsic_caps_20261014_endpoint_bounds.md) proves that every cap supplying a finite or ideal endpoint on an original cusp edge misses the same half-vertex point, by ball power greater than 1/4. It includes all 16,272 finite-pool occurrences, every real radius, every isolated dihedral label and the parallel limits. Endpoint caps alone cannot complete this seed. Interior caps are necessary if it is retained.

The [separate bound review](fresh_h20_intrinsic_caps_20261014_endpoint_bounds_review.md) independently checks all 260 radius maxima, exceptional right-contact cases, critical margins and tail polynomial identities. It source-reviews the transported full census and the annular argument, with explicit scope limits.

The [two-probe obstruction](fresh_h20_intrinsic_caps_20261014_quarter_gate.md) exactly excludes the uniform 252-cap quarter-vertex bank of radius 3/2. It explains the move to face-center placements; it does not exclude those placements.

The [common annular escape](fresh_h20_intrinsic_caps_20261014_annular_gate.md) is stronger for the subsequent small-radius search. The point `3z_0/4` escapes every combination of radius-at-most-one old vertex caps, radius-at-most-one caps centered in half the alcove, radius-at-most-3/2 caps centered in a quarter of the alcove, and every radius-at-most-two cap meeting the central sphere at any finite acute angle or tangency. Arbitrary selection and replacement within those families do not help. Larger endpoint caps lie outside this theorem.

The numerical [221,592-cap endpoint scout](fresh_h20_intrinsic_caps_20261014_completion_scout.json), [quarter-face selection](fresh_h20_intrinsic_caps_20261014_quarter_face_scout.json), and [central-contact scout](fresh_h20_intrinsic_caps_20261014_central_contact_contiguous_scout.json) are proposal mechanisms only. The selected 60 quarter caps have subsequently been certified exactly. The central scout finds 8,262 individual candidates compatible with all 587 walls, but none cuts the eighth-vertex probes; the annular theorem rules out completion using their entire radius family. A bounded joint replacement MILP reached its 40-round limit and remains computationally inconclusive. Its family is nevertheless excluded from finite-volume completion by the separate geometric theorem.

The [constructive scout report](fresh_h20_intrinsic_caps_20261014_constructive_scouts.md) records these searches and subsequent larger endpoint-cap tests. Two finite radius dictionaries fail numerical compatibility-and-coverage models requiring strict coverage of old alcove vertices. Those solver results do not exclude ideal completions with additional incident walls and bounded links, nor the complete continuous endpoint family. No finite-bank numerical infeasibility is promoted to a general theorem.

The [scale benchmark](fresh_h20_intrinsic_caps_20261014_scale_scout.md) checks 64 central-angle choices for the same three factors. A quarter-vertex escape persists even using the necessary acute upper bound on a uniform vertex-cap radius. The explicit parallel-neighbor radius is smaller; it is not claimed to be a maximum. Merely rescaling these factors does not repair this construction template.

The [new cusp benchmark](fresh_h20_intrinsic_cusp_benchmark_20261014.md) reconstructs 38 affine ADE factor geometries and all 607 rank-19 products with at most four factors. Among the 229 products with at most three factors, **A6 × A6 × A7** uniquely minimizes the largest vertex squared distance at unit inradii: all 392 vertices have square 196. This offers a balanced replacement for the old A11 direction. The four-factor comparison favors A4 × A5³, with 1,080 vertices of square 125. These are geometric benchmarks, not finite-volume constructions or optima over all scales and cusp types.

## Compact H9 from higher-dimensional bases

The [higher-base report](fresh_h9_highbase_construction_20261014.md) constructs three cyclic-base multiwedge families and fifteen five-wedge recipes on the simple 120-cell. All are explicitly polytopal before angle tests.

All 3,367 weight compositions over the three cyclic bases have a missing face of order above five. The fourteen tetrahedron-cluster recipes and the separate pentagon recipe on the 120-cell fail exact necessary Lannér and finite-vertex labeling conditions. One primary search hit a state limit and is preserved as inconclusive; a separate reordered exact join excludes that same recipe without truncation.

An independent checker reconstructs all 267,350 actual vertices across the fifteen 120-cell recipes and their ridge incidences. The [separate label review](fresh_h9_highbase_construction_20261014_review.md) verifies the all-label triangle compression and source logic; it does not claim a second whole label-join implementation. These exclusions concern the stated families only.

## Validation and continuing work

Exact producers and their saved-result replays pass. Numerical scouts retain their numerical scope. The known compact H8 calibration from the [preceding checkpoint](fresh_cusp_coverage_calibration_20261014.json) remains a separate verified reference; it has not been counted as a new example. The archive [audit](intrinsic_golden_cusp_construction_and_highbase_search_2026_10_14_audit.json) checks source hashes and local report links, not mathematical truth.

The next constructive step is to use larger endpoint caps jointly with interior caps, replacing conflicting old walls, or start from the balanced A6² × A7 cusp. Merely adding more members of the excluded small-radius families cannot close the chamber. A success must include a global finite-volume certificate and a literature novelty check.
