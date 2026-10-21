# Dotted-central tree case5: complete necessary-label exclusion

The saved case5 incidence has **444 vertices, nine missing triples and twelve missing pairs**, including its dotted central pair. Both the resumable continuation and a fresh independent traversal from the empty assignment have exhausted their queues with no necessary-label survivor. This excludes a compact hyperbolic9 Coxeter realization of this one incidence. The final conclusion does not depend on recovering the old100,000-node prefix. It does not classify the other tree incidences or all eighteen-facet polytopes, and it produces no new compact example.

The formerly [limited search](fresh_h9_other_trees_20261014_limited.md) and all its sources remain unchanged. The advance comes from three exact local negative-index gates and a resumable search, rather than a repeated run of that100,000-node prefix.

## Structural gates

All roots in these statements are normalized to square1. The saved certificates supply actual finite faces for every argument, with no assumed permutation symmetry between incidences.

| Gate | Conditions on two blocks and their leaves | Literal placements | Distinct conditions |
|---|---|---:|---:|
| [Dense/high](fresh_h9_dotted_central_20261014_gate.md) | A≥6, B≥7; first leaf has both own contacts nonright | 44 | 24 |
| [Strong contact](fresh_h9_dotted_central_20261014_additional_gates.md) | A≥4, B≥6; first leaf dense, one second-leaf contact≥4 | 88 | 88 |
| [Two dense leaves](fresh_h9_dotted_central_20261014_additional_gates.md) | A≥6, B≥4; both leaves have both own contacts nonright | 44 | 44 |

In each case the actual finite faces force enough foreign contacts to be right. Eliminating the positive four-root block then leaves a negative definite two-dimensional Gram matrix, impossible in Lorentzian signature. The last gate also handles its allowed foreign π/3 contact using Cauchy's inequality in the positive block. These proofs allow every integer label in the stated ranges and leave every unlisted contact and dotted distance free. The [dense/high independent review](fresh_h9_dotted_central_20261014_gate_review.md) and [complementary-gates independent review](fresh_h9_dotted_central_20261014_extra_gates_review.md) check the literal supports and exact coefficient bounds. Their generation and saved replays pass.

## Resumable traversal and completed scope

The first continuation reconstructs only the nine interrupted ancestor frames. At each frame, the saved child cursor follows the previously active child; completed earlier siblings are skipped and every later sibling remains pending. Each later resume independently reconstructs every saved assignment and propagated domain from its path, and checks each parent/child cursor relation.

| Batch | New traversal visits | Result |
|---|---:|---|
| [001](fresh_h9_dotted_central_20261014_search_batch_001.json) | 20,000 | 8 pending frames |
| [002](fresh_h9_dotted_central_20261014_search_batch_002.json) | 20,000 | 9 pending frames |
| [003](fresh_h9_dotted_central_20261014_strong_search_batch_003.json) | 50,000 | 9 pending frames |
| [004](fresh_h9_dotted_central_20261014_dense_search_batch_004.json) | 6,602 | Empty queue; no survivor |

The total is **96,602 continuation visits**, in addition to the inherited100,000-node prefix. This is not a count of distinct labelings: it includes up to nine recovered ancestor rechecks, and pruning events can reject partial families. No pending branch was discarded because of a node limit. The terminal batch stopped by exhausting its stack before its30,000-node budget, and its saved replay passes.

The solver uses175 exact hyperbolic-triangle categories, verified literal twin swaps, actual-vertex finite-diagram tests, orthogonal negative-block clauses, the five inherited analytic conditions and the156 new ones. Category7 denotes **every finite integer label at least7**. It is never substituted as a numerical angle in a metric computation: an isolated rank-two component permits the whole category, while a connected higher-rank finite component cannot contain it. Dotted distances remain unconstrained real variables. An empty necessary-label search therefore excludes all metric completions of this incidence.

The continuation reuses the frozen producer's finite-diagram classifier. It does not independently replay the old100,000-node prefix. The [independent traversal](fresh_h9_dotted_central_20261014_independent.md) starts from the empty assignment and finishes in **89,819 nodes**, with no cutoff, survivor or unexplored frontier. It independently reconstructs all444 vertices from the clutter and all175 triangle categories, checks the16 literal twin permutations, and uses its own immutable recursion and exact Fraction arithmetic in Q(√5). It checks278,137 finite masks through176 distinct rooted-tree Schur pivots. Generation and full saved replay pass. It has no producer imports and does not assume the resumed search is complete. The geometric realization and reviewed analytic theorems remain inherited dependencies.

```sh
python scripts/fresh_h9_dotted_central_20261014_dense_search.py --batch 4 --verify
python scripts/fresh_h9_dotted_central_20261014_independent.py --case 5 --verify
```

The next constructive opportunity is a different still-open tree incidence whose low-label blocks evade these conditions, followed by full omitted-contact and rank/signature completion of any necessary survivor. The three gates can be transferred only after checking their literal face premises. No further case or label search was started in this bounded phase.
