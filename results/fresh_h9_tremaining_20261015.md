# Remaining T-tree metric search: one exact family exclusion

T10 produced a complete necessary-label assignment, but its exact Gram constraints exclude every metric in that family. An independently reviewed all-label proof was bound into the resumable search. The next assignment is preserved as **open at 9,022 cumulative visits**, with ten pending DFS frames. No new compact H9 example or whole-incidence exclusion is claimed.

The [scope ledger](fresh_h9_tremaining_20261015_ledger.json), [ledger checker](../scripts/fresh_h9_tremaining_20261015_ledger.py) and [saved cursor](fresh_h9_tremaining_20261015_case10_metric_search_case_10_batch_002.json) record the final state. All preceding checkpoint files remain unchanged.

## Selection and exact local coverage

The preceding ledger left 18 representatives open. The three selected cases had 148 inherited analytic predicates each. Before searching, the new binders checked literal three-large-block and two-high-block theorem supports, both orderings of high blocks, and all 16 simultaneous root-pair swaps on complete vertex sets and minimal-nonface clutters. No isomorphism between the three cases is assumed.

| Representative | Actual vertices | Initial predicates | This phase |
|---|---:|---:|---|
|10|444|151|Searched; one metric family excluded; next family open|
|19|444|151|Prepared, not searched|
|29|488|153|Prepared, not searched|

The binders are [T10](fresh_h9_tremaining_20261015_gates_case_10.json), [T19](fresh_h9_tremaining_20261015_gates_case_19.json) and [T29](fresh_h9_tremaining_20261015_gates_case_29.json). They retain the uncovered three-large-block placements `023` and `123` explicitly. Category `7` always represents every finite integer at least 7; disjoint-facet distances remain continuous variables. The finite classifier and bridge clauses are necessary filters, not metric realizations.

## The first metric family

The [initial traversal](fresh_h9_tremaining_20261015_search_case_10_batch_001.json) stopped at 9,008 visits. Its root-pair labels were `(4,4,>=7,4)`. A [contact-domain screen](fresh_h9_tremaining_20261015_case10_domains.json) found 48 forced-right omitted pairs and 71 optional omitted meeting pairs. In particular the eight root normals form a positive direct sum of four dihedral blocks.

The [high-label proof](fresh_h9_tremaining_20261015_case10_high_bound.json), checked by a [separate SymPy review](fresh_h9_tremaining_20261015_case10_high_bound_review.json), excludes every high label at least 12 by a strictly negative definite residual plane. It leaves only labels 7–11, without imposing a numerical cutoff.

The [full metric proof](fresh_h9_tremaining_20261015_case10_rank_gate.md) and [certificate](fresh_h9_tremaining_20261015_case10_rank_gate.json) use 100 necessary projection/contact patterns. Each three-normal residual determinant is affine in `d=1/(4sin²(pi/m))` with nonzero slope. Exact Chebyshev equations eliminate all but four label-8 cases; each remaining case has two negative directions. This covers every omitted contact and dotted distance consistent with the saved mandatory family.

The [independent review](fresh_h9_tremaining_20261015_case10_rank_gate_review.md), [checker](../scripts/fresh_h9_tremaining_20261015_case10_rank_gate_review.py) and [data](fresh_h9_tremaining_20261015_case10_rank_gate_review.json) reconstruct all 100 determinants, 500 Chebyshev tests, four negative planes, two upper-contact witnesses and four actual face budgets. It also independently verifies the 16 literal twin maps and two distinct transported family predicates. This is an independent metric-family proof, not an independent traversal of the entire T10 search tree.

## Preserved next state

The [153-predicate composition](fresh_h9_tremaining_20261015_case10_metric_gates.json) adds the two reviewed images and hits the saved survivor. The [resumed search](../scripts/fresh_h9_tremaining_20261015_case10_metric_search.py) validates all old frame paths and sibling cursors before proceeding. It stops after 14 new visits, using about 38.26 seconds including source validation, at the next unexcluded necessary assignment. The unused 99,986-visit allowance is intentional: the configured stop is the first survivor.

The only changed mandatory label is `(12,13):4→5`, giving internal labels `(4,4,>=7,5)`. This state is not covered by the completed rank theorem. It has not been promoted to a Gram realization, and the exact symbolic high-label tail is preserved. The next productive target is a coupled low-block generalization: the two-anchor determinant has the form `(u−1)(d−2)−1`, where `u=1/[2(1−cos(pi/n))]`. That observation is a proposed next route, not a proved new exclusion.

Generation and saved replay pass for the binders, both search batches, the contact screen, the primary metric certificates and the ledger. The independent high-label and full rank reviews also pass their saved replays. The prepared empty-start independent DFS is not run because the search is nonempty. No jobs remain running in this phase.

The catalogue totals therefore stay at 15 excluded representatives covering 164 of 240 recipes, with 18 representatives and 76 recipes open. These are the verified tree-cube construction recipes; they are not all eighteen-facet compact H9 combinatorial types.
