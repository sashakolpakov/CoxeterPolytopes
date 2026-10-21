# Complete compact-Coxeter exclusion of the literal T8 incidence

**The 18-facet T-case8 incidence has no compact H9 Coxeter realization.** The final necessary-label search exhausts all branches in 59,423 nodes, below its 100,000-node limit. Its saved replay passes with identical results. There is no surviving assignment or unexplored frontier. This closes one incidence from the [33-type tree catalogue](fresh_h9_tree_cube_wedges_20261014_summary.json), not that whole catalogue or the compact-H9 existence problem.

The [terminal search certificate](fresh_h9_t8_sparse_search_20261014_outside_case_08.json) records the traversal and proof dependencies. The [fast checkpoint](fresh_h9_t8_sparse_search_20261014_checkpoint.json) checks all four search stages, the successive survivor exclusions, terminal completeness and recursive provenance: 32 JSON manifests and 67 unique source references. No metric candidate remains to fit in this incidence.

## The new structural gates

The first continuation used the previously reviewed eight-condition sparse theorem at partial assignments. For an unassigned pair, it intersects the categories remaining in every triangle domain containing that pair. A theorem prunes a branch only if **every** remaining category satisfies its hypothesis. Empty intersections are left to the ordinary constraint search; they are not treated as theorem applications. Free own contacts need not have been assigned. The saved traces include example domain witnesses and a hash over all theorem-prune records.

The first outside assignment changed the internal label on `(4,5)` from four to five. The [all-block extension](fresh_h9_t8_sparse_search_20261014_all_blocks.md) proves the same sparse obstruction for every finite label on that pair at least four. It adds 384 exact polynomial cases over `Q(sqrt(2),sqrt(5))`, with 5,472 strictly positive Bernstein coefficients; their minimum is `261/2000 - 31sqrt(5)/1000`. The 432 complementary inherited occurrences and every larger internal label are accounted for. A [separate Fraction review](fresh_h9_t8_all_blocks_20261014_review.md) reconstructs the inverse/projection data, determinant coefficients, inherited partition and monotone domains.

The next assignment entered a different regime: two high internal blocks, with wall 14 touching one root in each. The [two-high-block theorem](fresh_h9_t8_sparse_search_20261014_two_high_blocks.md) then removed all fixed outside-contact assumptions. Under only

```text
m(1,2)>=7, m(4,5)>=7, m(7,8)>=4, m(12,13)>=4,
```

literal actual-vertex Schur budgets force two residual vectors to span a negative-definite plane. This contradicts Lorentz negative index one. All other meeting labels and all dotted distances are free. The [independent review](fresh_h9_t8_sparse_search_20261014_two_high_blocks_review.md) rebuilds the supports, null-vector matrices and all-label projection-cost argument without the earlier polynomial bank.

The same proof also applies to high blocks `(1,2)` and `(7,8)`, with the other two internal labels at least four. This uses wall 15 instead of 14. The [literal bridge-transfer certificate](fresh_h9_t8_sparse_search_20261014_bridge_variants.md) checks every proof premise afresh, for anchors 11 and 17: all 45 core meeting pairs, 24 cross-block guards, six budgets, four stars, the positive eight-wall support and three missing triples. No global symmetry is assumed. A [separate transfer review](fresh_h9_t8_two_high_transfers_20261014_review.md) independently checks all four placements.

One remaining apparent survivor only reversed the active root of `(1,2)` in the sparse theorem. The [exact twin certificate](fresh_h9_t8_sparse_search_20261014_twin.json) verifies that swapping roots 1 and 2 preserves all 462 actual vertices and all 21 minimal nonfaces. It transfers the entire reviewed theorem to the opposite orientation. The final search includes both orientations and both two-high placements.

## Why the search is complete for this incidence

The saved realizable simple incidence has ten missing triples and eleven missing pairs. In a compact Coxeter realization each missing triple is Lannér; its three pairs meet, so their labels are finite. Missing pairs remain arbitrary strict dotted pairs throughout. No distance is substituted or bounded.

Every missing-triple assignment is enumerated using categories `2,3,4,5,6,>=7`. The category represented by 7 denotes every finite integer label at least seven. It is sufficient exactly: a triangle with a high label and both other labels two is spherical; otherwise its reciprocal sum is at most `1/2+1/3+1/7<1`. A connected finite diagram of rank at least three cannot contain a label at least six, since an adjacent nonright edge gives a three-node determinant at most `1-3/4-1/4=0`. Isolated rank-two components allow every finite label, including the entire high category. Thus there is no finite-label cutoff in the triangle or vertex tests.

The depth-first search enumerates every consistent assignment to the ten missing triples, propagating shared pair labels and testing all affected actual vertices. Its only orientation reductions swap literal twin facet pairs, each checked to preserve the complete vertex set. The four chosen triples are `[0,1,2]`, `[3,4,5]`, `[6,7,8]` and `[11,12,13]`; their leaves 0,3,6,11 are fixed by every other twin swap, so the four orderings can be imposed independently. Arbitrary unassigned meeting contacts are temporarily set right for necessary finite-diagram tests. This is safe: Perron–Frobenius monotonicity for the normalized acute Gram matrix prevents increased off-diagonal magnitudes from repairing a nonpositive smallest eigenvalue.

There are also 60 pairs of disjoint induced missing faces. When every cross contact is forced right by an actual-vertex failure at the smallest nonright label three, the two missing-face negative directions are orthogonal. This gives two negative directions for every dotted-distance choice. These are necessary negative-index exclusions, not numerical metric failures.

Finally, the four new theorem gates prune partial domains only by the universal test described above. After these sound reductions, the traversal has no remaining branch. A compact realization would induce an enumerated mandatory-label assignment and survive every necessary test, so the empty traversal excludes it.

## Recorded runs and validation scope

| Search stage | Nodes | Outcome |
|---|---:|---|
| [Original sparse gate](fresh_h9_t8_sparse_search_20261014_case_08.json) | 1,942 | First survivor; excluded by the all-block extension |
| [All-block continuation](fresh_h9_t8_sparse_search_20261014_next_case_08.json) | 4,232 | First survivor; excluded by the two-high theorem |
| [Both high-block placements](fresh_h9_t8_sparse_search_20261014_final_case_08.json) | 21,331 | Reversed sparse orientation; excluded by the twin transfer |
| [Both sparse orientations](fresh_h9_t8_sparse_search_20261014_outside_case_08.json) | **59,423** | **Complete empty traversal** |

These are overlapping successive traversals, not disjoint counts of assignments. Every run used a 100,000-node cap; none reached it. All four generation and saved-replay runs pass.

The last traversal records 238,748 rejected spherical-vertex trials, 32,665 orthogonal-negative-block prunes, 3,420 prunes for each sparse orientation and 7,296 plus 3,648 prunes for the two high-block placements. Trial counts can exceed node counts. All intermediate survivor files remain as historical checkpoints and are explicitly superseded by the terminal result.

The full traversal replay reuses the established finite-diagram classifier and search implementation; it is not an independently implemented full DFS. The [separate coverage review](fresh_h9_t8_sparse_search_20261014_outside_coverage_review.md) independently reconstructs the 462 vertex incidences, 175 triangle categories, all 16 simultaneous twin swaps, 60 bridge clauses, 180 sampled Schur failures and 30 partial-theorem domains. Its source review finds no DFS coverage defect; its generation and saved replay pass. It explicitly does not repeat all 59,423 search nodes. The new analytic exclusions have separate arithmetic/proof reviews, and the literal second-bridge transfers have a separate incidence review. The fast checkpoint only validates provenance and terminal consistency. These verification scopes are distinct.

```sh
python scripts/fresh_h9_t8_sparse_search_20261014_outside.py --case 8 --limit 100000 --verify
python scripts/fresh_h9_t8_sparse_search_20261014_checkpoint.py --verify
```

All work in this phase uses the new `fresh_h9_t8_sparse_search_20261014` prefix; the previous frozen theorem and catalogue files are unchanged. There is no live job or remaining T8 state. A subsequent constructive step should move to a different still-open incidence, rather than another parameter family inside T8. For example, the original catalogue retained necessary survivors at cases 17, 22, 25 and 26; this result does not classify those cases.
