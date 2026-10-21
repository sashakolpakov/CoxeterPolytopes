# Finite-angle restriction and algebraic cusp directions

No new compact or finite-volume Coxeter polytope is established. This checkpoint narrows the literal Borcherds-chamber section route, supplies a complete finite algebraic census for its remaining rank-three one-omitted cusp branch, and excludes two more compact-H9 incidence families. The broader construction goal remains open.

The literature restriction remains decisive: noncompact arithmetic reflection groups occur only in dimensions at most 19 or in dimension 21. Searching integral lattice hosts for noncompact arithmetic examples in dimension 20 or above 21 cannot succeed. Compact arithmetic examples instead have the upper bound 29. The [dimension audit](arithmetic_dimension_bounds_and_effective_high_dimension_search.md) and [algorithm audit](vinberg_allcock_search_strategy_20261013.md) retain the primary references and explain why Allcock's published integral-lattice algorithm must be distinguished from a nonarithmetic diagram-and-metric search inspired by edge walking.

## Exact results and their scope

| Route | Result | Scope |
|---|---|---|
| H20 sections, arbitrary real normal | Every actual finite ridge must have label 2, 3, 4, or 6 | The literal 210-wall source; no arithmeticity conclusion |
| One omitted cusp facet, affine rank at least four | All 175 omission/direction pairs, representing 119 spatial planes, are excluded for all real admissible parameters | All 45 such affine factors at the eighteen recorded cusp representatives |
| One omitted cusp facet, affine rank three | Complete census of 91 exact local direction records | Global ridges remain untested for this census |
| One-contact A2-to-B2 edge sections | Eight exact closing-contact recipes; three complete continuous families excluded exactly | Broad 5,048-family scout is numerical |
| Compact H9, prior sixteen proposals | All have type `Delta2^3 × I^3`, excluded by the published product classification | Every angle labeling and dotted distance for this type |
| Compact H9, sixteen-facet nine-gon multiwedges | All 51 admissible dihedral weight types fail exact necessary diagram constraints | This particular globally polytopal incidence family |

The [fresh Bugaenko calibration](fresh_frontier_calibration_20261014.json) again passes: compact H8, eleven facets, 41 finite vertices, zero cusps, with the exact Vinberg completion checked. It calibrates the verification machinery; it is a known example.

## H20 finite angles

The [arbitrary-normal angle theorem](fresh_h20_crystallographic_angle_gate_20261014.md) extends the preceding one-contact restrictions to every actual finite ridge of a hypothetical finite-volume Coxeter section. For a label outside `{2,3,4,6}`, the dihedral face is compact. Its possible boundary contacts reduce to the literal positive-definite three-root catalog, together with the preceding complete H3 obstruction. The surviving possibilities force a smaller compact face whose bounding source reflections preserve a rational integral lattice, contradicting rational isotropy in that dimension.

The [separate review](fresh_h20_crystallographic_angle_gate_20261014_review.md) found no gap and independently reconstructs the six new exact interval tables and algebraic endpoint inequalities. It explicitly retains the earlier H3 proofs as dependencies. The theorem does not make the full section arithmetic and does not constrain arbitrary dotted distances. Deleting source walls or changing the ambient chamber falls outside its premise.

## Cusp directions

The [rank-at-least-four classification and obstruction](fresh_h20_omitted_factor_classification_20261014.md) examines 387 possible omissions and 93,921 contact supports. Its positive-kernel and boundary-angle equations give 175 accepted omission/direction pairs. Of the 119 spatial planes, 56 match earlier excluded planes. The remaining 63 include 48 directions over `Q(sqrt(6))`. Exact interpolated ambient face points establish genuine section ridges throughout the parameter intervals, with two isolated actual-vertex certificates covering the exceptional parameter and endpoint. One angle-gap case invokes the new finite-angle theorem.

A [separate local verifier](fresh_h20_omitted_factor_classification_20261014_local_review.md) reconstructs all 175 accepted records and their shift intervals. The [separate global witness review](fresh_h20_omitted_factor_classification_20261014_global_review.md) reconstructs 126 ambient face points, 208 polynomial sign certificates, endpoint redundancy, and both isolated actual vertices. Exhaustiveness of the local support enumeration and the geometric interpretation also receive source/proof review; the independent checkers do not repeat every earlier theorem on which the exclusion depends.

The [rank-three calculation](fresh_h20_rank3_omitted_sections_20261014.md) uses a different reduction because affine A2 triangles allow three contacts of the same sign. Prescribing the planar rays of an affine A2, C2, or G2 triangle gives two rational conics in positive relative scales. Exact resultants and real-root isolation solve all 120 systems. The output contains 16 A2, 19 C2, and 56 G2 records; 44 use quartic representations, fourteen quadratic ones, and 33 rational ones. These count omission/target/direction records, not projectively deduplicated mirrors or new polytopes.

Every rank-three record contains the algebraic normal, positive kernel, null shift, and full local interval `[-L,0)`. This is the next concrete section branch to test globally. Its quartic records explain why rational and quadratic direction banks alone cannot establish completeness here. Other local facet counts and sections avoiding the recorded cusps remain separate.

The [independent rank-three verifier](fresh_h20_rank3_omitted_sections_20261014_review.md) passes generation and saved replay. It rederives the 120 conics and uses direct polynomial saturation, independently of the producer's resultant construction, to check every positive branch and all 91 normal/shift identities. It imports no producer code, while sharing SymPy for rational computer algebra. Global ridge validity remains outside both checks.

## Finite edge closure and compact incidences

The [B2 closure calculation](fresh_h20_b2_edge_sections_20261014.md) proves that a further facet of the B2 face must realize one of eight contacts, each allowing at most two positions in the old edge's Lorentzian plane. Three representative families are checked exhaustively with exact quadratic-tower arithmetic: all 1,062 transverse recipes fail at actual ridges. The [source review](fresh_h20_b2_edge_sections_20261014_review.md) explains the finite-volume premise and the role of the numerical witness packet.

The broad sweep finishes all 5,048 anchored families: 870,063 numerical defects and 237 remaining proposals. All 237 remainders are exactly old-endpoint double roots. This leaves no numerical survivor, but does **not** certify the completeness or every rejection of the floating-point sweep. It must not be reported as an exact exclusion of all 5,048 families.

For compact H9, the [incidence report](fresh_h9_nonproduct_incidence_20261014.md) identifies the previous sixteen proposals as a forbidden product type, using equality of their entire maximal-face lists. Its different constructive branch builds actual Euclidean nine-polytopes by seven wedges over a rational nine-gon. All 6,435 positive weight compositions reduce to 51 dihedral types after the Lannér missing-face bound. Exact mandatory Lannér joins are empty for every type, without a state truncation, metric grid, or SAT premise. The [separate review](fresh_h9_nonproduct_incidence_20261014_multiwedge_review.md) checks the all-label relaxation and enumeration scope; it is not a second join implementation.

## What remains open

The immediate new work is the global ridge test for the 91 rank-three algebraic records. Remaining multi-contact edge profiles and constructions that change the ambient assembly are not excluded by the angle theorem or the cusp cases proved here. Compact-H9 work must change the starting incidence family again. No general H20 or H9 nonexistence theorem follows from this checkpoint.

The [preceding checkpoint](cusp_classification_and_noncrystallographic_face_obstruction_2026_10_14.md) remains the source for the all-retained cusp classification and the H3/H4 obstruction. Historical bounded scouts stay labeled as such; the new exact results replace only the explicitly covered branches.
