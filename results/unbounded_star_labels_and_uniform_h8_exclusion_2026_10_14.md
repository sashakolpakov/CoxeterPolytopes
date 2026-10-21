# Unbounded star labels and uniform H8 attachment exclusion

**No new finite-volume Coxeter polytope or dimension is established.** This phase closes broad last-block contact families in the compact-H9 star search, excludes the fixed raw H8 attachment in every cap direction, and builds a different fifteen-facet chamber with three disjoint caps. The [preceding checkpoint](star_contact_families_and_h8_attachment_obstructions_2026_10_14.md) and its frozen artifacts are preserved.

## Compact H9: unrestricted last-block labels

The literal eighteen-facet star case 0 remains the target: 378 actual vertices, eight missing triples and fifteen dotted pairs. With the first three blocks' own contacts fixed at pi/3, the [combined theorem](fresh_h9_star_lastblock_complete_20261014.md) excludes **every last internal label at least four and every last own-contact choice**. Other finite meeting labels and dotted distances remain unrestricted in this necessary obstruction.

The proof combines three separately reviewed families:

| Last-block family | Exact census | Independent review |
|---|---|---|
|[Internal label 4 or 5](fresh_h9_star_lastblock_20261014.md), arbitrary own contacts|290 own-profile pairs; 77,952 residual determinants. The 32 zero cases fail literal vertices.|[Complete reconstruction](fresh_h9_star_lastblock_20261014_review.md)|
|[Shared single pi/3 contacts](fresh_h9_star_shared_high_20261014.md), every internal label at least 7|7,840 affine cases: 7,392 positive throughout the continuous parameter range; 448 fail literal vertices.|[Separate Fraction arithmetic](fresh_h9_star_shared_high_20261014_review.md)|
|[Opposite single pi/3 contacts](fresh_h9_star_opposite_high_20261014.md), every internal label at least 7|47,040 quadratic cases: 43,412 coarse Bernstein certificates, 16 refined certificates and 3,612 literal vertex failures.|[Separate five-point interpolation](fresh_h9_star_opposite_high_20261014_review.md)|

For internal labels at least six, four literal strong-contact guards remove any own label at least four. Two dense-contact guards remove own profiles `(3,3)`. Exact missing-triple inequalities then leave precisely the shared and opposite single-contact cases. This supplies complete label coverage without a finite cutoff or a restriction on the field containing the large-label cosine.

The search resumes the old 51-visit position through saved batches at 103, 105 and 186 cumulative visits. The [intermediate 186-visit cursor](fresh_h9_star_lastblock_complete_20261014_search_batch_006.json) changes leaf 16's own profile on roots `(9,10)` to `(3,4)`; all other own profiles are `(3,3)`, with internal labels `(4,4,4,5)`. This lies outside the initial fixed-first-three-own family.

A subsequent [five-family dependency review](fresh_h9_star_dense_frontier_20261014_dependency_review.md) proves that ten noncore own contacts may increase to arbitrary labels at least three, while facet 14's own contacts stay exactly three and all four internal labels stay at least four. All ten changed entries lie outside the eleven-root residual Gram. Their other uses are lower-bound gates or vertex failures preserved by Perron monotonicity. This extends the proof to the 186-visit assignment without repeating the 155,680 inherited arithmetic cases.

The [widened search](fresh_h9_star_dense_frontier_20261014.md) resumes that cursor and reaches **415 cumulative visits**, with nine pending frames and every unexplored sibling retained. The [saved batch](fresh_h9_star_dense_frontier_20261014_search_batch_007.json) adds 229 visits and the new gate prunes 84 states. Its next assignment has internal labels `(4,4,4,5)` and all own profiles `(3,3)` except facet 16's `(2,5)` on roots `(9,10)`. Its right own contact lies outside the new lower-bound hypothesis. This remains a necessary label assignment; no metric realization or whole-incidence exclusion is claimed.

The separate [24 conditional cross-right guards](fresh_h9_star_last_internal_three_20261014_cross_right_review.md) prepare a possible internal-three analysis. They require the reduced dense-own hypotheses; the same cross contacts are not forced by the weaker first-three-block hypotheses alone. No complete internal-three exclusion is asserted here.

## H8 attachment: every cap direction is excluded

The [composed attachment theorem](fresh_h9_h8_raw_attachment_20261014.md) excludes the raw construction consisting of the frozen mixed twelve-wall H9 chamber, one all-dotted cap and the eleven perpendicular sidewalls of the fixed compact Bugaenko H8 source. The result covers every real cap direction and every isometric placement of that source.

Finite volume requires the twelve old-wall balls to cover the placed source. The earlier all-real source contact-gap theorem says that a projected norm below `B=sqrt(5)-1` must lie on one of three exceptional source-corner rays. The new proof exhausts every possibility:

* [At most one center can be below B](fresh_h9_h8_corner_anchors_20261014.md), by all-label pair incompatibilities, with a [separate arithmetic review](fresh_h9_h8_corner_anchors_20261014_pair_review.md).
* A low center at old walls 0–10 leaves too few balls able to reach the other exceptional corners. The [uniform radial theorem](fresh_h9_h8_global_cap_gap_20261014_radial.md) has a [separate physical-coordinate review](fresh_h9_h8_global_cap_gap_20261014_radial_review.md).
* A low center at polar wall 11 forces a small wall-4 contact. Four literal roots then force cap square greater than `2651/918>2`, contradicting its normalization. The [projection theorem](fresh_h9_h8_corner_anchors_20261014_polar_gate.md) has a [separate polynomial and interval review](fresh_h9_h8_corner_anchors_20261014_polar_review.md).
* If every center is at least B, the two source vertices at distance `acosh(17+8sqrt(5))` must be covered by an old dotted pair. Exactly three such pairs exist, and each violates a positive-support projection bound. The [all-high theorem](fresh_h9_h8_all_high_20261014.md) and its [independent algebraic-field review](fresh_h9_h8_all_high_20261014_review.md) reconstruct the source vertices and all 66 old pairs.

The [exact new cap](fresh_h9_h8_global_cap_gap_20261014.md) passed earlier scalar size tests, but its [specific radial obstruction](fresh_h9_h8_corner_anchors_20261014_radial_gate.md) showed that no orientation covers the source. The uniform theorem subsumes both that rejection and the preceding checkpoint's continuous-line exclusion. Numerical placement failures are not premises of the theorem.

A further [side-lift rigidity lemma](fresh_h9_h8_side_tilts_20261014.md), with an [independent quadratic-field review](fresh_h9_h8_side_tilts_20261014_review.md), shows that tilting the eleven source-side walls cannot evade the obstruction: every actual ambient sidewall inducing one of this source's facets must meet its containing hyperplane at a right angle. The proof uses the literal pi/3 graph and five exact angle gaps, covering all finite labels.

These are scoped results. Extra ambient caps disjoint from the source facet can supply more covering balls. Changed retained walls and different source geometry also remain open. Neither theorem excludes arbitrary larger completions or compact H9 polytopes in general.

## A changed H8-based chamber with three disjoint caps

The [new starting chamber](fresh_h9_h8_disjoint_cap_20261014.md) changes the retained core. Extend the source form by the positive coefficient `d=sqrt(5)/4`, retain its eleven sidewalls and bottom facet, and add three caps based on the exceptional source-corner rays. All three caps coexist and have squared norm one quarter. They meet their eight incident source walls at right angles and are strictly disjoint from every other retained wall.

The result has **15 actual facets, exactly 44 finite vertices and 167 finite edges**. The bottom H8 face retains its 41 vertices and 164 edges. The caps add three finite vertices and close three formerly open edges. There remain 62 finite-corner directions whose completion is unresolved. A strict spacelike escape with all wall inequalities strict proves infinite volume.

The [independent review](fresh_h9_h8_disjoint_cap_20261014_review.md) reconstructs all 105 wall pairs, all 165 source supports, the complete finite-vertex and finite-edge lists, every strict facet witness and the escape. The coefficient has nonsquare determinant ratios to the two specified earlier coordinate extensions. A two-wall cycle equals `4sqrt(5)` and has negative conjugate, so any finite-covolume completion retaining that pair would be non-quasi-arithmetic. This supplies a concrete new starting point for edge extension without entering the previously rejected arithmetic hosts.

## Verification and continuation

The exact producers, saved replays and independent reviews pass in the scopes stated above. The [archive audit](unbounded_star_labels_and_uniform_h8_exclusion_2026_10_14_audit.json) binds their source manifests and report links, and preserves the preceding frozen archive. Archive integrity alone is not a mathematical proof.

Continue the star search from the 415-visit cursor, treating the new right own contact explicitly. For constructive work, pursue the 62 unresolved finite-corner directions of the changed fifteen-wall chamber, with exact global compatibility and volume checks for every extension. Repeating raw placements or side tilts of the old mixed chamber cannot improve their established outcomes.

```sh
python3 scripts/audit_unbounded_star_uniform_h8_checkpoint_20261014.py --verify
```
