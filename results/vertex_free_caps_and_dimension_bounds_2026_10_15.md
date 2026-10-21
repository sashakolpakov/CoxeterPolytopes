# Geometric extensions and the parallel dimension-bound programme

No new finite-volume Coxeter polytope or dimension is established. The construction lane reaches an exactly certified 29-facet infinite-volume chamber in H9. The new dimension-bound lane gives a checked refinement of Prokhorov's published proof to **dimension at most 990**, but Nikulin's older announcement of **dimension below 300** prevents a new-best-bound claim. A separate compact argument reduces a possible exclusion of dimension 29 to explicit local targets; those exclusions remain unproved.

## Keep the hypotheses visible

All rows concern hyperbolic Coxeter polytopes or their reflection groups. Finite volume includes compactness unless the row says noncompact. These are upper bounds, not existence claims at their endpoints.

| Hypotheses | Bound / current status |
|---|---|
| Compact, arbitrary arithmeticity | Dimension at most 29, by Vinberg |
| Finite volume, arbitrary arithmeticity | Prokhorov's published 995 argument yields 990 with the refinement below; Nikulin previously announced the stronger bound below 300 |
| Arithmetic and finite volume | At most 29, combining the compact and noncompact bounds |
| Arithmetic and noncompact | At most 21, with dimension 20 also excluded |
| Compact and right-angled | At most 4 |
| Compact and flag/3-free | At most 12 |
| Simple and ideal | At most 8; neither qualifier may be dropped |
| Compact, every irreducible vertex component of rank at most r | At most 4r by the derived counting argument, also subject to the general bound 29 |

The [hypothesis matrix and primary references](dimension_bound_lane_20261015.md) give the precise statements and further facet-count restrictions. The [arithmetic literature audit](arithmetic_dimension_bounds_and_effective_high_dimension_search.md) explains why arithmetic dimension 20 is unavailable and why 21 is not an upper bound for all compact arithmetic examples. The [algorithm comparison](vinberg_allcock_search_strategy_20261013.md) separates lattice-based Vinberg enumeration from edge-walking and its completeness requirements.

The reflection/Coxeter condition is essential. Removing it leaves compact and noncompact arithmetic hyperbolic manifolds in arbitrarily high dimensions. For every n≥3, the forms `Σx_i²−2x_(n+1)²` over Q and `Σx_i²−√2 x_(n+1)²` over Q(√2) give nonuniform and uniform arithmetic lattices, respectively; torsion-free finite-index subgroups give manifolds. [Gelander–Levit, §§4.1–4.2](https://arxiv.org/pdf/1401.8003) supplies these constructions. Compactness, finite volume and arithmeticity alone do not give a dimension ceiling.

## An unconditional refinement, with its historical limitation

The [Prokhorov refinement](dimension_bound_prokhorov_20261015.md) preserves his weights, good-face inequalities and treatment of cusps. For a spherical edge diagram of rank N≥17, its weighted pair sum satisfies

`W ≤ (29N−173)/3`.

The retained additive deficit is enough: with `N=n−1`, the published lemma gives `n<996−5536/(n−1)`. At `n=991+t`, the necessary polynomial becomes `t²+985t+586<0`, impossible for `t≥0`. Thus the published lemmas imply **n≤990** without arithmetic, facet-count, label-cutoff or simple-cusp restrictions.

The [independent arithmetic review](dimension_bound_prokhorov_20261015_review.md) reconstructs distances on spherical graphs and the scalar inequalities without producer imports. The [primary-source proof audit](dimension_bound_prokhorov_20261015_proof_review.md) checks that ideal endpoints do not change the spherical edge rank, that the marked pairs are unordered, and that the unchanged good-face estimates apply.

[Nikulin's ICM 1986 report, printed pp. 654 and 659](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf), announces `dim H<300` for general finite-covolume reflection groups. It supplies a narrow-face outline and an affine-A example, but not the full constants and affine-case calculation. We do not dismiss that announcement or describe 990 as the strongest known bound. The next broad target is an auditable reconstruction of the announced stronger argument, including all affine component products and the all-finite narrow-face alternative.

Three proof obligations are explicit: derive separation estimates for every affine type and product; derive the finite-vertex functions, iteration lengths and finite-only alternative; then verify complete numerical casework contradicting the `arccosh(7)` narrow-facet bound for every n≥300, including the framing-wall hypotheses. The affine-A example alone does not close these obligations.

## Compact bounds and a concrete dimension-29 target

The [compact component inequality](dimension_bound_lane_20261015.md) counts nonright intrinsic polygon angles at vertices. If all irreducible spherical vertex components have rank at most r, the exact partition maximum and Nikulin's strict average inequality imply `d≤4r`. The [independent review](dimension_bound_vertex_components_20261015_review.md) verifies the partition and scalar arithmetic; the proof has algebraic tails, not a finite dimension cutoff.

For the additional full flag/3-free hypothesis, a quadrilateral charge to actual diagram edges gives stronger bounds. Merely having no triangular two-faces is insufficient. These conditional results do not improve the unrestricted compact bound.

The [dimension-29 reduction](fresh_compact_dimension_bound_20261015.md) identifies nine rank partitions and **2,126 typed spherical forests** that could carry the required excessive small-face weight. At least 350 relevant vertices are forced by the exact face-count identity. The [independent review](fresh_compact_dimension_bound_20261015_review.md) checks the arithmetic and the occurrence bridge obtained by zeroing weights on polygons with at least five sides. For each listed forest, ruling out triangle/quadrilateral angles of total marked weight 1, 2 or 3, according to its total weight, would suffice. None of these local exclusion programmes is claimed complete.

Vinberg's explicit E7 obstruction prevents a uniform improvement using nonnegative diagram-local weights that vanish across components. Rank-specific context, global incidence constraints or a different weighting scheme remain possible. Deleting support nodes can expose unrestricted I2 labels, so the finite list of entry types does not justify a finite label cutoff in the completion search.

## Four exact walls, but the volume obstruction persists

The [geometric construction report](vertex_free_geometric_caps_20261015.md) adds four norm-two walls to the frozen 25-facet chamber. The first two arise from complete geodesic faces without finite endpoints, a frontier absent from finite-vertex walking alone. Later walls cut the replacement escapes in succession. Every addition has exact pair checks, actual-facet witnesses and an independent reconstruction.

| Actual facets | Finite vertices | Ideal vertices | Finite edges | Open finite-corner germs | Hyperideal rays |
|---:|---:|---:|---:|---:|---:|
| 25 | 81 | 0 | 247 | 235 | 513 |
| 27 | 84 | 0 | 247 | 262 | 660 |
| 28 | 87 | 0 | 249 | 285 | 727 |
| 29 | 88 | 0 | 249 | 294 | 776 |

The [latest independent cone review](fresh_h9_fourth_escape_extension_20261015_review.md) verifies all 406 pair relations and the complete cone. A new exact strict spacelike escape proves infinite volume. More walls and more finite vertices are not a completion certificate.

The full exact frontiers at 25, 27 and 28 walls have independent reviews. The numerical scouts and 33-recipe exact contact pool have their stated bounded scopes; they are not complete wall enumerations or all-label nonexistence arguments. New coefficient fields, different retained walls and caps unrelated to endpoints remain open.

The [replayable next frontier](fresh_h9_after_fourth_rescore_20261015.md) saves 152 numerical escape cutters among 67,056 numerically compatible proposals from the bounded 68,111-row bank. Two exact nine-contact recipes require an additional quadratic extension of the coefficient field. Their full chamber compatibility and endpoint geometry remain unproved; neither is yet a certified new cap.

## Scoped obstructions and the saved compact search

The [A1×A1 direct H18 coupling theorem](fresh_h20_a1a1_h18_coupling_20261015.md), with its [independent review](fresh_h20_a1a1_h18_coupling_20261015_review.md), excludes the fixed-source direct family in H20 for all finite labels. A common strict escape works for the complete relaxed 168-normal bank. Additional noninducing walls and different sources are outside the theorem.

The [rank-three H17 coupling classification](fresh_h20_h17_rank_three_coupling_20261015.md), with its [independent review](fresh_h20_h17_rank_three_coupling_20261015_review.md), covers all spherical rank-three direct inducing attachments of the specified sources, including symbolic dihedral tails. Only the ordinary cases and twelve isometric A1-cubed variants survive the local classification, and exact escapes exclude finite volume for these direct variants. This is not an exclusion of arbitrary H20 polytopes.

The [opposite-anchor single-third theorem](fresh_h9_star_single_third_opposite_20261015.md) and its [review](fresh_h9_star_single_third_opposite_20261015_review.md) advance the compact-star search to **397 predicates and 60,015 cumulative visits**, saved in [batch 014](fresh_h9_star_single_third_opposite_20261015_search_batch_014.json) with nine pending frames. Its next necessary label assignment is not a metric realization. The [third-label-four theorem](fresh_h9_star_third_four_20261015.md) and [review](fresh_h9_star_third_four_20261015_review.md) cover that next state's family but have not yet been bound into a resumed search. Preserve the saved stack and siblings.

The [inertia ledger](fresh_h9_lorentz_inertia_20261015.md) distinguishes local configurations impossible in every Lorentz dimension from rank-only H9 obstructions. Using these local exclusions for a dimension theorem still requires a separate argument forcing their exact face and label hypotheses.

## Archive and next work

The [archive audit](vertex_free_caps_and_dimension_bounds_2026_10_15_audit.json) checks recursive source hashes, local links and preservation of the [frozen predecessor](polar_extensions_and_discriminant55_obstruction_2026_10_15_audit.json). It is an integrity check; mathematical scope belongs to the individual proofs and reviews.

The recursive H17 provenance reached an old temporary ancillary-file path. The original [Alexeev arXiv source archive](https://arxiv.org/src/2209.09110v4) was retrieved again, and the ancillary text matches the historical SHA-256 exactly. A [permanent copy](../data/references/alexeev_2209_09110v4_all_cases.txt) and [retrieval provenance](../data/references/alexeev_2209_09110v4_provenance.json) now replace that temporary dependency in this audit; no frozen predecessor was edited.

Continue three distinct tasks: reconstruct the stronger general dimension argument; test the compact dimension-29 small-face targets with exact unrestricted-label coverage; and seek a finite-volume closure or changed construction in the missing dimensions. For construction, compact H9 and nonarithmetic finite-volume H20 remain central targets. Do not restart higher-height searches in the two already excluded literal arithmetic lattices. No finite-volume discovery, global compact exclusion, or new-best dimension bound is claimed by this checkpoint.
