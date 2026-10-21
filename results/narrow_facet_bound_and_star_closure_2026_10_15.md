# Narrow-facet bounds, compact-star closure and nine-contact geometry

The nonexistence lane now has a complete, internally reviewed research proof candidate for an unconditional dimension bound. The construction lane has a new exact 30-wall chamber, but its volume is still infinite. No new finite-volume example or dimension is established.

The [preceding checkpoint](vertex_free_caps_and_dimension_bounds_2026_10_15.md) and its frozen certificates remain unchanged. The [archive checker](../scripts/audit_narrow_facet_star_checkpoint_20261015.py) and [manifest](narrow_facet_bound_and_star_closure_2026_10_15_audit.json) bind this checkpoint, its reports and their recursive dependencies. Archive integrity is distinct from the mathematical proof and numerical-search scopes described below.

## A complete dimension argument

The [reviewed research proof candidate](fresh_nikulin_dimension_189_candidate_20261015.md) gives **n≤189** for every finite-sided, finite-volume real hyperbolic Coxeter polytope. It imposes no arithmeticity, compactness, simple-cusp or facet-count hypothesis. This is a draft mathematical result with internal independent checks; publication priority and current-best status are not claimed.

Nikulin's narrow-facet theorem supplies a facet whose framing hyperplanes have pairwise distance below arccosh7. There are two exhaustive cases:

| Chosen narrow facet | Argument | Ambient bound |
|---|---|---:|
| Contains an ideal vertex | Affine-product cusp energy, common null-height ratio and a global contact budget | 189 |
| Contains only finite vertices | It is compact; quasi-Lannér extraction and coarse graph-distance weights on its two-faces | 143 |

The frozen [225 proof](fresh_nikulin_dimension_225_candidate_20261015.md) first established the complete dichotomy. Its [product-cusp argument](fresh_nikulin_narrow_20261015_products.md) covers every affine product, including arbitrarily many small factors, using at most four retained nonright contacts per terminating wall. The [classical formulas](fresh_nikulin_narrow_20261015.md), [independent arithmetic](fresh_nikulin_product_independent_20261015.md) and [hostile geometric review](fresh_dimension225_hostile_20261015.md) check the constants, scales, endpoints and actual-ridge bridge. Those files remain unchanged.

The [sharper contact-cost proof](fresh_nikulin_energy_charge_20261015.md) charges each non-narrow component's energy deficit to twice its actual Schur cost. The two endpoint costs sum to at most four; the narrow component retains a separate loss allowance of two. Exact exceptional calculations handle the inexpensive E7/E8 contacts that would defeat a minimum-diagonal estimate alone. Thus total energy is at least M(2N/9−10), where N=n−1 and M≥1, and the framing-wall pairing satisfies −(x,y)≥M(N/9−7). For n≥190 this is at least14, contradicting narrowness. The [separate review](fresh_nikulin_energy_charge_20261015_review.md) verifies this accounting, all exceptional identities and the A9 boundary case with independent Fraction arithmetic. Generation and saved replay pass.

The [compact-facet lemma](fresh_nikulin_narrow_20261015_compact_facet.md) deliberately permits affine subdiagrams. Compactness of a face does not imply that its ambient star is free of affine diagrams: explicit tetrahedra disprove that attempted localization. The repaired proof uses minimal indefinite diagrams of order at most ten, a path bound of17, and the simple-polytope h-vector inequalities on the facet. It does not assume the facet is itself Coxeter.

The literature distinction remains important. Nikulin's1986 ICM report already announced n<300; the prior990 refinement of Prokhorov's explicit proof was therefore not a new-best claim. The current narrow-facet reconstruction supplies explicit constants and written arguments. Its classification and narrow-facet inputs are identified in the reports with primary-source page references.

The [hypothesis ledger](dimension_bound_lane_20261015.md) continues to separate compactness, finite volume, arithmeticity, right angles, flagness and simplicity. The new general candidate does not supersede the stronger established bounds in restricted classes: compact Coxeter polytopes have n≤29; arithmetic finite-volume reflection groups also have n≤29; the noncompact arithmetic class has n≤21 and excludes n=20. These are reflection-group restrictions, not bounds on arbitrary arithmetic hyperbolic manifolds. See the [arithmetic literature audit](arithmetic_dimension_bounds_and_effective_high_dimension_search.md).

## Complete closure of the saved compact star family

The [three-large-block theorem](fresh_h9_star_three_large_blocks_20261015.md) proves that at most two of the four indicated internal labels can be at least6. Its proof uses actual spherical supports and negative index one, and places no fixed bounds on the leaves' own labels. The [independent review](fresh_h9_star_three_large_blocks_20261015_review.md) checks the literal hypotheses.

Four new placements bring the saved star0 binder to497 predicates. The [primary search and report](fresh_h9_star_large_block_resume_20261015.md) exhaust all remaining assignments after166,030 cumulative visits. A separate [empty-start traversal](fresh_h9_star_large_block_resume_20261015_independent.md), with independent Fraction/tree-Schur arithmetic and no cutoff, returns empty after163,894 visits. Both saved replays pass. The [hostile source and catalogue audit](fresh_h9_star_closure_hostile_20261015.md) additionally checks symbolic labels, bridge relaxation, symmetries and all120 recipe isomorphisms.

Together with the frozen exclusions of star representatives1–4, this closes all **120 star-tree recipes** in the240-recipe catalogue. It is a complete finite-family exclusion, not a theorem excluding all compact H9 polytopes or all eighteen-facet incidences.

The next [T-tree case6](fresh_h9_t6_frontier_20261015.md), representing two recipes, is also excluded. Its binder transfers three of the four new block placements after checking their literal supports; it does not assume that the uncovered triple023 satisfies that theorem. Both the primary search and a separate empty-start traversal exhaust all assignments in14,092 visits, with158 predicates,488 vertices and60 bridge clauses. Both saved replays and the [separate binder/source review](fresh_h9_t6_frontier_20261015_gates_review.md) pass. This leaves **21 of the28 T-tree representatives** outside the established exclusions.

## A barrier to the proposed compact-dimension29 shortcut

The [small-face constructions](fresh_compact29_small_faces_20261015.md) give79 exact local stars covering every marked pair of A4, B4, D4, F4, H4, E7 and E8. There are56 triangular and23 quadrilateral faces. A [separate implementation](fresh_compact29_small_faces_20261015_review.md) reconstructs the faces and checks all25,457 principal tests and4,864 dotted-edge shortcuts.

After orthogonal spherical padding, every positive-weight marked pair in the rank29 forest E7+E7+E7+E8 admits an individual small-face star. Consequently the proposed strategy of excluding enough marked angles independently cannot eliminate that forest using only those local conditions. This does not realize a compact H29 polytope. A stronger obstruction must constrain several incident faces simultaneously or use other global information.

## Exact construction and the renewed escape

Two [octic alternatives](fresh_h9_octic_escape_caps_20261015.md) separately cut the previous29-wall escape and give actual new facets and finite endpoints. They are mutually incompatible, and both retain one common exact escape. Their complete finite-angle alphabets have been [classified](fresh_h9_octic_angle_20261015.md). The associated [152-candidate diagnostic](fresh_h9_octic_angle_20261015_selection.md) is bounded and numerical.

The [29-wall frontiers](fresh_h9_octic_base_frontiers_20261015.md) expose a different cap determined by nine contacts and a quadratic norm equation. Its root lies in the original quartic field. This [exact nine-contact cap](fresh_h9_nine_contact_cap_20261015.md) cuts both recorded escapes. It would be missed by requiring ten linear finite contacts.

The resulting [30-wall geometry](fresh_h9_nine_contact_extension_20261015.md) has89 finite vertices, no ideal vertices,805 hyperideal extreme rays,250 finite edges and301 open finite-corner germs. Separate arithmetic and complete cone clipping independently reproduce the full census. A new strict spacelike escape proves **infinite volume**.

The [rescore of three frozen banks](fresh_h9_nine_contact_rescore_20261015.md) finds no escape cutter among their compatible rows. The [new30-wall frontiers](fresh_h9_nine_contact_frontiers_20261015.md) have eight new finite-corner germs and75 complete geodesic faces without finite endpoints. Their bounded scouts produce2,312 and71,375 numerical proposals, respectively, again with no escape cutter. None of these zero counts is an all-label or all-parameter nonexistence result.

The next construction step should test joint additions whose new vertices require multiple new walls, changed retained walls, or different sources. The existing single-endpoint proposal banks have now been tested against the latest escape; repeating them or merely increasing the facet count does not establish progress toward finite volume.

## Verification and remaining scope

The new exact cap and full cone have independent generation and saved replay. The product-cusp constants and compact-facet argument have separate arithmetic and hostile mathematical reviews. The star0 exclusion has independent complete traversal from the empty assignment and a separate source/catalogue audit. Numerical scouts remain explicitly bounded. The archive audit checks preservation of the preceding frozen sources, reports and top-level certificates.

The general construction objective remains open. The broad dimension argument is a research proof candidate, with its published inputs and new reductions exposed for further scrutiny. Stronger bounds, additional compact incidence exclusions and a finite-volume construction remain active research directions.
