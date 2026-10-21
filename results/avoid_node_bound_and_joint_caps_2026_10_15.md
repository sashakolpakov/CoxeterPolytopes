# Avoid-node dimension bound and joint-cap constructions

**Checkpoint: 2026-10-15.** The unconditional finite-volume dimension argument now has an internally reviewed proof candidate **n≤181**. Three more compact H9 tree incidences are excluded. The construction lane has an exact pair of H9 walls creating a finite corner that requires both additions, and a new H20 cap retaining the full H18 source. Every chamber constructed in this phase still has an exact infinite-volume escape. No new finite-volume example or dimension is claimed.

The [previous checkpoint](narrow_facet_bound_and_star_closure_2026_10_15.md) and all its sources remain frozen. The [archive checker](../scripts/audit_avoid_node_joint_caps_checkpoint_20261015.py) and [integrity record](avoid_node_bound_and_joint_caps_2026_10_15_audit.json) preserve that predecessor, validate recursive sources and local links, and bind the current entrypoints. Integrity checks are separate from mathematical certification.

## Dimension bounds with their hypotheses

The [181 package](fresh_nikulin_dimension_181_candidate_20261015.md) contains the theorem statement, full proof chain and dependency binder. The [new local lemma](fresh_nikulin_avoid_node_20261015.md) presets omissions avoiding the narrow-facet node in every affine type. Its [independent review](fresh_nikulin_avoid_node_20261015_review.md) reconstructs exact exceptional energies and inverse diagonals and checks all analytic cases.

For the narrow affine factor of rank r, the old loss of two is replaced by

    energy/M + 2C ≥ 2(r−1)/9.

The non-narrow factors retain the stronger target `2r/9`. Across the two terminating endpoints, total Schur cost is at most four. With `N=n−1`, the Lorentz identity gives

    −(x,y) ≥ M(N−55)/9,    M≥1.

At n≥182 this reaches 14, contrary to Nikulin's narrow-facet theorem. The other case, where the selected narrow facet is compact, retains the separately reviewed bound 143. No arithmeticity, simple-cusp, label-cutoff or facet-count assumption enters the combined181 candidate.

This is an internally reviewed research proof candidate, not an externally validated theorem or a priority claim. Prokhorov's published Theorem C gives 995; Nikulin's earlier ICM report already announces the stronger unrestricted n<300. The bounded literature check found no later stronger unrestricted result, but cannot establish its absence. [Prokhorov](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&paperid=1484&what=fullteng), [Nikulin, printed pp.654–659](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf).

The [hypothesis ledger](arithmetic_dimension_bounds_and_effective_high_dimension_search.md) continues to separate the construction targets:

| Hypotheses | Constraint used in this project |
|---|---|
| Compact hyperbolic Coxeter polytope | Dimension at most 29; compact H9 remains a target |
| Arithmetic finite-volume Coxeter polytope | Dimension at most 29 |
| Noncompact arithmetic Coxeter polytope | Dimension at most 21; dimension 20 is also excluded |
| Finite-volume Coxeter polytope, arbitrary arithmeticity | New reviewed research candidate 181; historical published 995 and announced below 300 remain explicitly distinguished |
| Hyperbolic manifold without a reflection/Coxeter requirement | These reflection dimension obstructions do not apply |

Thus H20 is pursued through nonarithmetic geometry. The two previously proved nonreflective H9 arithmetic hosts remain closed; higher root heights cannot repair those hosts.

## Three more compact H9 incidences

The [compact report](fresh_h9_tnext_20261015.md) closes T-tree representatives 20, 11 and 31, with unrestricted finite labels represented symbolically. Cases 20 and 11 exhaust in 18,992 and 15,213 primary visits; independent empty-start exact searches reproduce both exclusions. A separate [literal-premise and routing review](fresh_h9_tnext_20261015_scope_review.md) checks the three initial predicate binders and their case selection.

Case 31 first produced a full necessary-label survivor. Five literal maps of an independently reviewed two-high-block theorem satisfy 400 actual-face premises and 15 missing-triple occurrences. The resulting two new predicates exclude that survivor's entire metric regime by forcing a two-dimensional negative subspace in a Lorentz Gram matrix. The [separate transfer review](fresh_h9_tnext_20261015_case31_transfers_review.json) checks the maps without importing the producer.

The resumed case 31 primary exhausts its saved cursor; a fresh independent traversal finishes in 18,615 visits with no survivor or pending domain. Primary and independent saved replays pass. The [final ledger](fresh_h9_tnext_20261015_ledger.json) records **164 of 240 recipes excluded**, including all 120 star recipes. **Eighteen T representatives, comprising 76 recipes, remain unresolved.** These are exclusions in the specific frozen construction catalogue, not a dimension bound or a classification of arbitrary eighteen-facet polytopes.

## A joint H9 corner

The [joint-wall report](fresh_h9_joint_pair_20261015.md) records bounded rank-seven/two-anchor scouts, their numerical banks and two exact constructions. This search permits a new vertex to require multiple new walls.

The selected pair remains in the existing quartic field. Each cap passes every old-wall relation, both cut the previous escape, and their exact mutual pairing is−1 for square-two normals, giving angleπ/3. The alternative 32-wall chamber has **90 finite vertices, no ideal vertices, 250 finite edges and 310 open finite-corner directions**. Its unique new finite vertex uses both new walls; its nine directions remain open. The [independent two-cut review](fresh_h9_joint_pair_20261015_two_review.md) reconstructs the complete cone with separate arithmetic. A new exact strict escape proves infinite volume.

A simpler alternative polar cut produces31 actual facets,90 finite vertices and251 finite edges, also with a complete [independent review](fresh_h9_joint_pair_20261015_review.md) and a new strict escape. The polar has positive pairing with each selected joint cap, so the two constructions are alternatives. The polar closes an old edge, demonstrating why earlier bounded zero-cutter scouts must not be interpreted as exhaustive cap exclusions.

The numerical scout banks contain 412 and 4,566 proposal records. Their40 mutual finite-angle fits in the larger bank are not all certified; rounded records can duplicate a geometric wall and near-parallel fits can be false positives. Only the separately reconstructed exact caps have the stated mathematical status.

## A new H20 cap outside the direct-attachment exclusions

The [H20 report](fresh_h20_joint_route_20261015.md) constructs a norm-two cap retaining all 40 old H18 source walls and both orthogonal support walls. Its source projection Y has square−3/2, and the complete cap is

    R = (Y, −√3/2, −1).

Because Y is future timelike, this cap is strict throughout the H18 source face. It is a noninducing wall and lies outside the earlier fixed-source direct inducing-wall obstructions.

All 43 facets and 903 pair relations are checked exactly. The cap has 18 right contacts, three π/3 contacts, one π/6 contact, one parallel relation and 19 disjoint hyperplanes. Exhausting its 1,540 possible finite-vertex supports yields exactly one finite vertex, connected to an old vertex by a finite edge. The former direct-coupling escape is cut, but a new exact strict escape proves that this 43-wall chamber still has infinite volume. The [independent review](fresh_h20_joint_route_20261015_exact_review.md) uses separate exact arithmetic. Only vertices on the new cap are exhaustively counted; no complete ambient H20 cone census is claimed.

## Next tasks

Preserve this phase's proofs, banks and certificates. Continue the 18 remaining compact T types with literal predicate transfers and independent exhaustion checks. For geometric H9, explore the newly created joint corner and its open directions, or change retained geometry; adding facets alone is not a volume argument. In H20, the exact noninducing cap establishes a useful route beyond the direct-source exclusions, but the new escape needs a different cap family or source projection.

For the dimension lane, 181 now needs scrutiny of the whole mathematical argument and a broader literature/priority check. Its local 2/9 avoidance loss is sharp for the uniform lemma because of affine A1; that observation does not rule out stronger global accounting or a different method. Compact dimension 29 remains a separate global multi-face problem: the frozen single-angle local extension barrier is unchanged.
