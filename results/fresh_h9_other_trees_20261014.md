# Four other tree incidences have no compact H9 Coxeter realization

Cases **17, 22, 25 and 26** of the saved 18-facet tree catalogue are now excluded. All four complete necessary-label searches finish below their 100,000-node limits. Separately implemented traversals, using exact forest Schur elimination over `Q(sqrt(5))`, independently reproduce all 211,165 nodes and every rejection count. There is no surviving assignment or unexplored frontier in these four cases.

This result concerns four literal realizable incidence types. They are not assumed isomorphic to T8 or to one another. It does not classify the remaining catalogue, noncompact chambers, or all compact H9 polytopes; no new compact H9 example was found.

| Case | Actual vertices | Missing triples | Negative-block clauses | Complete search nodes |
|---|---:|---:|---:|---:|
| 17 | 438 | 10 | 48 | 48,022 |
| 22 | 462 | 10 | 59 | 47,372 |
| 25 | 432 | 11 | 38 | 53,274 |
| 26 | 432 | 11 | 38 | 62,497 |
| Total | | | | **211,165** |

## Literal proof transfers, not a symmetry assumption

The starting files contained one necessary mandatory-label survivor per case. Each is excluded by a locally checked transfer of the [two-high-block negative-plane theorem](fresh_h9_t8_sparse_search_20261014_two_high_blocks.md). The [new transfer checker](../scripts/fresh_h9_other_trees_20261014_transfers.py) and [certificate](fresh_h9_other_trees_20261014_transfers.json) check every proof premise in each target incidence.

There are 17 valid two-high maps: 4, 5, 4 and 4 in the four cases. These check all 45 core meeting pairs, the positive rank-eight support, all cross-block witnesses, every active-contact budget and upper-angle star, and the required literal missing triples. The same exact negative-plane calculation then applies to the renamed variables; other labels and dotted distances remain free.

The stronger [sparse all-block residual theorem](fresh_h9_t8_sparse_search_20261014_all_blocks.md) has four verified placements in case 17 and eight in case 22. Every one of its 55 core pairs, domain supports, own-contact guards and replacement stars is checked. There is **no** such transfer in the tested placements for cases 25 and 26. Those two cases therefore were searched using the negative-plane transfers alone; no missing rank-three hypothesis was assumed.

The 29 sufficient maps give 6, 10, 3 and 3 distinct theorem condition sets. These are sufficient checked maps, not a claim to enumerate all possible embeddings of the proofs. The [independent transfer review](fresh_h9_other_trees_20261014_transfer_review.md) reconstructs all four complete vertex sets and separately checks the 29 injective maps, 2,524 actual-face premises, missing triples, mapped label conditions and deduplication witnesses. Its generation and saved replay pass. The original analytic arithmetic remains bound through its frozen independent reviews.

## Complete label coverage

The [primary search](../scripts/fresh_h9_other_trees_20261014_search.py) uses the exact categories `2,3,4,5,6,>=7` on every missing triple. Category seven represents every finite integer label at least seven. The reciprocal-sum test is uniform on that category; isolated finite rank-two components retain all labels. Missing pairs remain arbitrary strict dotted distances. No finite label cutoff or choice of distance is used.

The depth-first search propagates every shared triangle pair and tests the actual finite vertices. Unassigned meeting pairs are relaxed to right angles for necessary tests. Its negative-block test forces cross pairs right only when the smallest nonright label three already fails at an actual vertex. Disjoint Lannér blocks with all cross pairs right then give two negative directions, independently of dotted distances.

The transferred theorems are used at partial assignments only when every remaining category for every hypothesis pair satisfies the theorem. Their own contacts may remain unassigned. Literal twin swaps are checked against the complete vertex set before reducing triangle orientations. The independent implementation checks all 16 simultaneous swaps and that the four ordering leaves are outside every swapped pair.

The final primary files are:

- [Case 17](fresh_h9_other_trees_20261014_search_case_17.json)
- [Case 22](fresh_h9_other_trees_20261014_search_case_22.json)
- [Case 25](fresh_h9_other_trees_20261014_search_case_25.json)
- [Case 26](fresh_h9_other_trees_20261014_search_case_26.json)

Every one has `EXACT_ALL_MANDATORY_LABELS_AND_BRIDGES_EMPTY`, no survivor and no unexplored frontier. Since a compact Coxeter realization would supply a surviving mandatory-label assignment, these are incidence exclusions, rather than merely failures of selected labelings or numerical metric fits.

## Separate full traversal and arithmetic

The [independent DFS](../scripts/fresh_h9_other_trees_20261014_independent_dfs.py) adapts the frozen independent T8 engine to these different incidences without importing producer/project code. It reconstructs each vertex set from the missing-face clutter, generates the 175 ordered Lannér-triangle categories using Fraction arithmetic, and uses immutable assignments and its own propagation/bridge checks. It has no node cutoff and compares counts with the primary only after completing the recursion.

Its finite-diagram predicate uses exact forest Schur elimination, not the primary's Dynkin lookup. A cycle has a nonpositive all-ones witness; a connected rank-at-least-three component with a label at least six has a nonpositive three-node path. The remaining tree pivots use squared cosines `1/4`, `1/2` and `(3+sqrt(5))/8`, with exact signs in `Q(sqrt(5))`. Isolated rank-two components accept every finite label.

All four independent complete runs pass and match every primary node and prune total. They check respectively 72,299, 55,865, 48,678 and 54,042 distinct finite masks. Their certificates are [case 17](fresh_h9_other_trees_20261014_independent_dfs_case_17.json), [case 22](fresh_h9_other_trees_20261014_independent_dfs_case_22.json), [case 25](fresh_h9_other_trees_20261014_independent_dfs_case_25.json) and [case 26](fresh_h9_other_trees_20261014_independent_dfs_case_26.json). Primal realization and the independently reviewed analytic theorems are explicit inherited dependencies; their arithmetic is not reimplemented by this traversal.

The [fast checkpoint](fresh_h9_other_trees_20261014_checkpoint.json) binds the transfer review and all eight full-search certificates and verifies terminal completeness, source hashes and count agreement. Its generation and saved replay pass. This is a consistency check, distinct from the already completed full independent traversals.

```sh
python scripts/fresh_h9_other_trees_20261014_independent_dfs.py --case 25 --verify
python scripts/fresh_h9_other_trees_20261014_checkpoint.py --verify
```

No shared primary traversal was repeated merely to obtain another identical count. Instead, all four were independently traversed with different finite-Gram arithmetic. Earlier frozen files are unchanged. The separate [case-5 continuation](fresh_h9_other_trees_20261014_limited.md), whose central opposite pair is dotted, reached its explicit 100,000-node cap without a survivor. It remains unresolved and is outside this four-case exclusion certificate.
