# Direct rational Gale constructions for compact H9

No compact hyperbolic metric survives this bounded construction. Two walks through explicit integral rank-four Gale configurations produced **770 distinct realizable simple nine-dimensional polytope types with 14 facets**. Exact necessary Coxeter gates exclude every saved type. **749 are absent from the union of both previously closed Bugaenko surgery catalogues**; this comparison does not assert that they are new relative to every earlier repository experiment.

This is a sampled construction and exclusion result. It does **not** exhaust rank-four Gale diagrams, any open missing-edge graph shadow, or all compact H9 polytopes with 14 facets. The older surgery theorems are not used as a classification of that larger problem.

The starting points come from the eight still-open graph shadows `0,1,3,4,5,7,9,31` in the [H9/14 ledger](h9_14_facet_completion.md). Their designated old nerves were already excluded; the purpose here is to change the actual nerve by moving Gale vectors. The walk is allowed to change the missing-edge graph as well. Burcroff's [Lemma 10.10](https://arxiv.org/pdf/2201.03437) motivates rejecting missing triangles; every retained proposal has only missing faces of orders 2, 4 and 5 and at least two missing pairs.

## Construction and scope

Each proposal is an integral `14 × 4` matrix. All 1,001 four-row determinants must be nonzero. Positive five-element circuits and their complements must each cover all fourteen labels. Their complements are the proposed nine-vertex facets; complete minimal-nonface reconstruction supplies the incidence filter.

The first walk starts from the rational realizations of the eight old selected nerves. Numerical SVD produces a starting Gale proposal, followed by integer rounding and exact acceptance. It makes 10,000 coordinate-change proposals per starting point. The second walk starts from a matrix with few missing four-faces in each first-walk group and makes 12,500 proposals per group, favoring fewer missing four-faces. The exponential acceptance rule guides proposals only; it is not part of any exclusion proof.

| Bounded batch | Proposed moves | Accepted moves, including repeated nerves | Distinct new nerves | Duplicates in the two surgery catalogues | Outside their union |
|---|---:|---:|---:|---:|---:|
| Initial walk | 80,000 | 69,718 | 351 | 15 | 336 |
| Few-four-face walk | 100,000 | 80,883 | 419 | 6 | 413 |
| Total | 180,000 | 150,601 | **770** | **21** | **749** |

The second batch excludes the first batch's canonical keys before saving, so the 419 records are additional. It includes 71 new nerves with **no** missing four-face, compared with four in the first batch. Thus it changes the relevant overlap structure rather than merely extending the first walk.

All coordinates satisfy `abs(entry) ≤ 8192`; the bound `24·8192^4 < 2^63` justifies the primary determinant expansion in signed 64-bit integers. The separate replay uses arbitrary-size FLINT integers and rational matrices.

## Exact polytopality

For every positive circuit, choose its positive integral cofactor dependence. Their sum gives a vector `w > 0` with `Gᵀw = 0`. Rescaling row `i` by `w_i` makes the row sum zero. The transpose kernel has dimension ten and contains the constant vector, so it defines a rational affine configuration of fourteen points in dimension nine.

Gale duality identifies its facets exactly with the complements of the positive circuits. Uniformity makes each such facet an **8-simplex with nine vertices**. Every label lies on a facet and is consequently an actual vertex. This is an actual rational polytope realization, stronger than feasibility of an abstract oriented-matroid sign table. The polar gives the simple nine-dimensional polytope whose Coxeter metric is sought.

The independent checker goes further than this argument: it constructs rational homogeneous primal coordinates, verifies every supporting functional and every facet's affine rank eight, reconstructs the full minimal-nonface clutter, and checks degree two at every ridge. Across the 770 cases this covers **770,770 integer determinants and 121,500 supporting facets**.

## Coxeter gates

| Exact terminal gate | Saved nerves |
|---|---:|
| Literal occurrence of one of the 44 previously verified induced-nerve obstructions | 715 |
| Two disjoint Lannér blocks with all proper unions spherical | 32 |
| Complete propagated full-nerve label search | 23 |
| **Total** | **770** |

The [local obstruction and relation certificates](fresh_h9_newfacet_frontier_20261014.md) and [two-block theorem](fresh_h9_bugaenko_remaining_20261014.md) are inherited, source-bound proofs. Every new occurrence is checked against the actual Gale-derived faces. The independent replay rechecks these literal embeddings without importing the producer's clutter or canonicalization routines.

The remaining 23 searches exhaust **189 states**, with a configured limit of 5,000 per nerve; none reaches that limit. They enforce shared Lannér labels, every actual spherical vertex, the complete 96-row overlapping-block relation, and early negative-index/rank tests. They cover all finite integer labels through the previously verified `2,…,6,≥7` categories. No complete necessary labeling remains, so there is no metric fitting or numerical Gram failure being promoted to an exact obstruction.

The saved search replay uses a separate finite-tree classifier and finite-field determinant implementation. It reuses the DFS, symmetry reduction and constraint propagation. The existing small-obstruction proofs are inherited rather than reproved in this new checker. The canonical relabelings are checked explicitly; the uniqueness convention itself uses the existing nauty implementation.

## Files and replay

- Proposal producers and literal matrices: [initial walk](../scripts/fresh_h9_gale_frontier_20261014_walk.py), [data](fresh_h9_gale_frontier_20261014_walk.json), [second walk](../scripts/fresh_h9_gale_frontier_20261014_low4.py), [data](fresh_h9_gale_frontier_20261014_low4.json).
- Exact realization and catalogue comparison: [checker](../scripts/fresh_h9_gale_frontier_20261014_geometry.py), [certificate](fresh_h9_gale_frontier_20261014_geometry.json).
- Exact label gates: [checker](../scripts/fresh_h9_gale_frontier_20261014_labels.py), [certificate](fresh_h9_gale_frontier_20261014_labels.json).
- Separate primal/arithmetic replay: [checker](../scripts/fresh_h9_gale_frontier_20261014_independent.py), [certificate](fresh_h9_gale_frontier_20261014_independent.json).

```sh
PYTHONPATH=src:. python scripts/fresh_h9_gale_frontier_20261014_geometry.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_gale_frontier_20261014_labels.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_gale_frontier_20261014_independent.py --verify
```

The first two commands regenerate their exact saved results. The last reconstructs every saved rational polytope and all literal gate applications, replays the 23 complete searches with separate arithmetic, and rejects forged Gale, positive-kernel and facet records. It does not replay all 180,000 proposal moves. No third batch or broader classification is included.

Both primary certificates passed generation. The separate checker passed generation and its final saved `--verify` replay, including all three forged-record rejection checks. All five scripts compile. No search or replay remains running.

A separate [source and theorem review](fresh_h9_gale_frontier_20261014_review.md) also passes. It checks the Gale realization, explicit primal construction, literal gate scope, counts and provenance; it does not claim a second full 770-matrix replay.
