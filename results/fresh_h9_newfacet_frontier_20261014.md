# The changed-order Bugaenko surgery family is closed

**All 5,536 admissible incidence types in the changed-order family are excluded as compact H9 Coxeter polytopes.** The 3,244 types left open at the preceding checkpoint now have complete proofs. No full necessary label assignment survives, so there is no metric candidate or new finite-volume example in this family.

The operations are precise: truncate a codimension 2–5 face of the fixed compact 11-facet Bugaenko H8, wedge the newly introduced truncation facet, then truncate a codimension 2–5 face of the resulting H9 polytope. Every construction is polytopal. This conclusion does not cover arbitrary larger-facet surgeries or a different source polytope. All preceding archives remain unchanged.

The [final proof ledger](fresh_h9_newfacet_frontier_20261014_checkpoint.json) and its [checker](../scripts/fresh_h9_newfacet_frontier_20261014_checkpoint.py) separate the disjoint proof classes:

| Proof | Incidence types |
|---|---:|
| Duplicates of the previously closed original-carrier family | 1,040 |
| Previously certified two-Lannér-block obstruction | 1,220 |
| Previous complete 32-case pilot | 32 |
| New induced-nerve obstructions | 3,112 |
| New complete global label searches | 132 |
| **Total** | **5,536** |

The underlying operation census has 11,729 safe ordered representatives: 5,503 fail the necessary missing-face order bound, and the remaining 6,226 recipes produce these 5,536 types. The [preceding catalogue and scope](fresh_h9_bugaenko_remaining_20261014.md) remain the source for operation-orbit completeness.

## Small reusable obstructions replace a broad search

The new search first extracted induced nerve configurations on at most nine nodes. For each, every minimal nonface must be a Lannér diagram and every maximal face must be spherical. Missing faces in the full catalogue and literal containment in actual vertices are checked before applying any local obstruction.

The [local certificate](fresh_h9_newfacet_frontier_20261014_patterns.json), produced by the [pattern checker](../scripts/fresh_h9_newfacet_frontier_20261014_patterns.py), records 160 targeted tests from a 2,065-pattern inventory. Forty-three are completely excluded and already eliminate 2,676 frontier types. Thirty-seven have a local necessary labeling, while eighty reach the 10,000-state discovery limit. Those latter outcomes are not exclusions or metric constructions.

One particularly effective additional obstruction has seven nodes: an order-five and order-four Lannér block share two nodes, with no other induced minimal nonface. Its complete finite-label join has no solution after 16,237 states. It occurs in 2,902 frontier types. Together with the other local obstructions it excludes 3,112 distinct types, leaving 132. The [exact overlap certificate](fresh_h9_newfacet_frontier_20261014_relations.json) and [checker](../scripts/fresh_h9_newfacet_frontier_20261014_relations.py) record this union explicitly.

These tests do not impose a label cutoff. Category 7 represents every finite integer label at least seven. Such a label can occur in a spherical diagram only in an isolated rank-two component; its Lannér-triangle condition is constant throughout that range. The complete order-four and order-five Lannér tables supply their usual finite labels. The induced configurations used here contain no dotted pairs.

## A complete local relation makes the remaining searches small

Two order-five Lannér blocks sharing four nodes, with no other induced missing face, have exactly **96** permissible complete finite-label assignments on their six-node union. This is a necessary relation, not a nonexistence result or a realization. It is much smaller and stronger for propagation than treating the two 420-row Lannér tables independently.

The [new engine](../scripts/fresh_h9_newfacet_frontier_20261014_engine.py) inserts every applicable six-node relation, filters tables against actual spherical vertices, and propagates the shared pair labels. Values common to every remaining row are forced. Search order uses the smallest surviving table domain. Unknown meeting labels are set right only for necessary positive-definiteness pruning; every value is subsequently branched unless forced or excluded.

All 132 remaining types close in **1,465 evaluated states**, with a maximum of **67** for one type. No case reaches its configured 3,000-state limit. The actual terminal proofs use empty Lannér joins, nonspherical vertices, and two orthogonal negative Lannér blocks; no numerical signature test or dotted-distance bound is used. Although the inherited exact fixed-minor gate remains enabled, it is not needed in this run.

The [search archive](fresh_h9_newfacet_frontier_20261014_search.json) and [runner](../scripts/fresh_h9_newfacet_frontier_20261014_search.py) preserve all case records and counters. A full necessary labeling would have triggered a numerical metric fit. None appears, so **no fit was run** and no numerical failure contributes to the theorem.

## Independent verification and its limits

The [separate local checker](../scripts/fresh_h9_newfacet_frontier_20261014_local_independent.py) and [audit](fresh_h9_newfacet_frontier_20261014_local_independent.json) replay all 44 local obstructions and all 96 relation rows with a different table traversal and an independent finite-tree classifier. They reconstruct maximal faces without the producer helper and verify every applied face against a saved actual vertex. The classical Lannér catalogue remains a shared dependency.

The [separate incidence replay](../scripts/fresh_h9_newfacet_frontier_20261014_incidence_independent.py) and [certificate](fresh_h9_newfacet_frontier_20261014_incidence_independent.json) rebuild all 5,536 recipes directly from the source vertices, checking 827,634 vertex occurrences. Minimal nonfaces are recovered from the entire downward face complex; the producer's clutter surgery and nauty routines are not used. Safe recipe-orbit completeness and final isomorphism deduplication retain the preceding producer/source-review scope.

The [separate search replay](../scripts/fresh_h9_newfacet_frontier_20261014_search_independent.py) and [audit](fresh_h9_newfacet_frontier_20261014_search_independent.json) reproduce all 132 complete case records using independent finite-diagram arithmetic. They reuse the DFS, canonicalization, relation propagation and branching control flow; they are not a second implementation of those mechanisms. The relation itself has the independent traversal described above.

The final ledger checks the complete, disjoint cover of all catalogue types and the source hashes. The bounded local discovery attempts therefore leave no unresolved incidence in the final family.

```sh
python scripts/fresh_h9_newfacet_frontier_20261014_checkpoint.py --verify
python scripts/fresh_h9_newfacet_frontier_20261014_local_independent.py --verify
python scripts/fresh_h9_newfacet_frontier_20261014_incidence_independent.py --verify
```

This result identifies the obstruction to this particular construction route: repeated overlapping large Lannér blocks impose rigid simultaneous label relations. A further incidence construction must change those induced configurations; merely changing labels or dotted distances on these 14-facet types cannot work. A third carefully chosen face surgery, or a different compact source nerve, remains outside the theorem.
