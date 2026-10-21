# Compact face closures and arithmetic triage

**No new finite-volume Coxeter polytope or dimension is established.** The H20 chamber now has 133 actual facets, at least 168 verified finite vertices and three preserved bounded cusps. Both tracked two-dimensional faces are compact, but a strict rational escape still proves infinite volume. Four further compact-H9 incidences are completely excluded. Three different coupled compact-face lifts were tested; two enter a classified nonreflective arithmetic lattice, while a mixed-source lift avoids that arithmetic obstruction but exhausts its extreme-ray polar extensions without closing.

This continues the [coupled-face and T8 checkpoint](coupled_face_completion_and_t8_exclusion_2026_10_14.md), whose frozen sources and reports remain unchanged. The [literature restrictions](arithmetic_dimension_bounds_and_effective_high_dimension_search.md) and [algorithm audit](vinberg_allcock_search_strategy_20261013.md) continue to apply. Noncompact arithmetic dimension 20 and dimensions above 21 are excluded; the compact/arithmetic bound 29 is a separate statement. Allcock's published edgewalking works in a specified lattice. A general combinatorial construction must additionally solve compatible metric data and prove global finite volume.

## H20: a second face closes, but the global escape survives

The [new construction](fresh_h20_s61_multicap_20261014.md) retains the previous 107 walls and adds a chain of 26 walls perpendicular to the S61 support. Their parameter is a uniquely isolated real algebraic root, approximately 16.766732018447345454205. All normal identities and right angles follow exact algebra; outward rational intervals certify strict inequalities. The 3,107 new pairs consist of 495 right and 2,612 strictly dotted pairs.

The retained three sides and new chain form the whole S61 face as a compact 29-gon. A positive relation bounds the Klein-chart polygon, and all 406 pairs of boundary lines are checked: exactly 29 are feasible future timelike corners, while 377 violate a recorded side inequality. This is a complete face census. The two inherited angles need not be Coxeter angles intrinsically; the claim is ambient Coxeter geometry and compactness of the face.

Exactly 27 finite vertices lie on the new walls. Together with the 141 retained known vertices this gives at least 168 overall, not a complete ambient vertex count. A common strict timelike center and acute-extension preservation establish every facet's actuality. All three bounded cusp links survive, and the already compact S75 hexagon is unchanged. The [independent review](fresh_h20_s61_multicap_20261014_exact_review.md) uses separate field/interval arithmetic and a different rotation evaluation to verify the root isolation, physical pairs, polygon, vertices, center, cusps and escape. Both primary and review pass saved replay.

The old rational escape still satisfies

```text
E² = 38017800175279 / 500000000000 > 0,
```

with strictly negative pairings against all 133 walls and the new center. The global chamber therefore remains infinite-volume.

A [new arbitrary-normal gate](fresh_h20_s61_multicap_20261014_perpendicular_gate.md) explains why continuing the same wall type cannot help. The projection of E to the S61 Lorentz plane is future timelike and strictly inside the compact 29-gon. Every coherent acute added wall strictly preserves all its finite vertices: their positive-definite active Grams have nonnegative inverses. Convexity then forces every added normal perpendicular to the S61 support to remain strictly negative on E. Thus no finite collection of such additions can remove the escape, with no restriction on field, norm or angle labels. Nonright support contacts and changed retained walls are outside this theorem.

The [separate gate review](fresh_h20_s61_multicap_20261014_perpendicular_gate_review.md) checks the proof and independently reconstructs the physical projection, all 29 strict side inequalities and all 29 positive-definite active Grams with nonnegative inverses. Its saved replay passes. The proof permits the proposed extension to use its own common future strict center; it does not assume the current chosen center survives.

## H9: four incidence exclusions and one explicitly unfinished traversal

The [other-tree package](fresh_h9_other_trees_20261014.md) excludes literal cases 17, 22, 25 and 26 of the saved eighteen-facet catalogue. The complete primary searches use 48,022, 47,372, 53,274 and 62,497 nodes, respectively: **211,165 nodes in total**. Every run terminates below its cap with no survivor or frontier. Separate complete traversals with no node cutoff reproduce every node and rejection count using independent exact forest-Schur arithmetic.

Twenty-nine literal proof transfers reuse the reviewed negative-plane and sparse residual theorems only where all their incidence hypotheses hold. A [separate transfer review](fresh_h9_other_trees_20261014_transfer_review.md) verifies 2,524 actual-face premises and all transferred conditions. Category `>=7` covers every finite label at least seven; missing pairs remain arbitrary dotted distances. These are complete necessary-condition exclusions of four specified incidences, not bounded-label failures or a general H9 theorem.

The [case-5 continuation](fresh_h9_other_trees_20261014_limited.md) has a different central dotted pair. It reached the explicit 100,000-node limit without a survivor and saved its unfinished path. It remains unresolved. The four complete exclusions must not be extended to this case.

## Three coupled compact-face experiments

| Compact source and support | Seed facets | Seed finite / spacelike rays | Extended facets | Extended finite / spacelike rays |
|---|---:|---:|---:|---:|
| Makarov P2 H5, coupled D4 | 11 | 14 / 16 | 15 | 49 / 89 |
| Eleven-facet H7 source, coupled B2 | 13 | 56 / 48 | 19 | 94 / 190 |
| Mixed Makarov H5, coupled D4 | 11 | 10 / 20 | 12 | 12 / 38 |

There are no ideal rays in these censuses. Every row represents exact physical walls, a common interior and actual facets. Each final chamber has a strict spacelike escape and infinite volume. Larger finite-vertex counts are not evidence of finite volume.

The [P2/D4 experiment](fresh_h9_d4_prism_coupling_20261014.md) changes a norm-one source facet to a norm-two ambient facet by subtracting a norm-one D4 outer weight. Its π/4 source contact becomes π/3 and its dotted source end pair becomes π/5. Four of sixteen polars pass and are mutually dotted. The independent reviewer uses different scalar arithmetic and cone clipping, reproducing all 138 final rays against the producer's 5,005-support enumeration.

The [H7/B2 experiment](fresh_dimension_new_mechanisms_20261014.md) similarly couples the source's only eligible even-angle facet to B2. Six of its 48 polars pass and are mutually dotted. Exact clipping yields all 284 final rays. Its only legal first source doubling has no remaining eligible facet. This package has primary exact generation and saved replay; it does not claim a separate arithmetic implementation for the full census.

Both arithmetic constructions admit explicit odd unimodular admissible Z[φ] lattice bases. All fifteen or nineteen reflection matrices preserve the respective literal lattice. Bugaenko's classification and nonreflectivity theorem, cited precisely in the construction reports, therefore exclude finite-volume completions using only reflections preserving those hosts. Rational-space equivalence alone would not justify this conclusion. Nonintegral reflections or different arithmetic data remain outside it, and integral-root growth in these two hosts is stopped.

The [mixed-source D4 experiment](fresh_h9_mixed_d4_coupling_20261014.md) deliberately changes that arithmetic data. Its seven-facet compact H5 source glues different Makarov prism blocks. Over L=Q(√5,√((16+4√5)/11)), replacing its norm-four end by half that normal minus the D4 outer weight gives a valid eleven-wall H9 chamber. The end cycle generates the quartic extension, while an embedding fixing √5 preserves an indefinite principal block. Vinberg's criterion implies that any finite-volume completion retaining these seed facets would be non-quasi-arithmetic.

This removes the identified arithmetic obstacle but does not complete the geometry. One of twenty seed polars passes. After adding it, every one of the 38 final spacelike-ray polars fails a recorded Coxeter pair. Thus the iterated extreme-ray polar procedure is exhausted after one addition. The complete field-dependent angle gate covers all possible finite labels for these candidates. Arbitrary new normals are outside this exhaustion.

The [independent mixed-source review](fresh_h9_mixed_d4_coupling_20261014_review.md) reconstructs the quartic field with separate Fraction arithmetic and exhausts all defining supports for source, seed and extension. Every ray, active set, norm, polar trial, final pair, facet and escape matches. Both implementations pass saved replay. The external arithmetic criterion is source-checked and kept separate from the exact scalar and geometry checks.

A [complete finite-corner frontier](fresh_h9_mixed_d4_coupling_20261014_frontier.md) makes the next geometric search concrete. All twelve finite vertices are simple. Their 108 directed edge releases consist of 56 directions forming 28 finite edges, and 52 open directions reaching the light cone before their spacelike projective endpoints. Each shared eight-wall support is checked for rank and positive definiteness. Forty-four open directions release a D4 wall; eight release a source wall at the two new polar vertices. This derivative reuses the frozen independent quartic arithmetic and complete cone census. It enumerates the existing finite-corner queue, not additional mirrors, edges having no finite endpoint, or every possible endpoint label.

## Verification and next searches

The [archive audit](compact_face_closures_and_arithmetic_triage_2026_10_14_audit.json) binds the new sources, certificates and reports, checks local links and preserves the preceding checkpoint. Integrity checks, exact proofs, independent arithmetic and bounded numerical scouting remain explicitly distinguished in each report.

For H20, search walls with nonright S61 contacts or change retained walls; further perpendicular completion is excluded. For compact H9, use other incidence types or solve general spherical endpoint extensions along the mixed-source chamber's 52 open directions. Such extensions need metric compatibility with every retained wall and may require fields outside the current quartic one; the polar procedure is already exhausted. The new arithmetic host identifications should be applied before any further lattice enumeration. Case 5 has a saved unfinished traversal, not an exclusion. None of these tasks has yet produced the requested new finite-volume example.
