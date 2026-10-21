# Accelerated exact finite attachment closure through dimension 29

The dimension-at-most-29 finite attachment queue is exhausted: **5,172 marked states, none pending, and 2,070 singular weight obligations**. The [terminal certificate](fresh_sqrt5_nonsingular_rank29_fast_20261015.json) and its [SQLite archive](../data/enumerations/fresh_sqrt5_nonsingular_rank29_fast_20261015_running.sqlite) are immutable. This is complete necessary attachment coverage from the specified finite seeds; weight payment and an improved dimension bound require separate certificates.

The [source](../scripts/fresh_sqrt5_nonsingular_rank29_fast_20261015.py) starts with all 1,649 stored states, including all 216 pending siblings, and all 311 literal seed maps. The starting batch has a separate [independent arithmetic and coverage review](fresh_sqrt5_nonsingular_rank29_exhaustion_20261015_review.md). The new registry uses exact marked full-label isomorphisms. Complete contact vectors are sorted lexicographically before insertion, preserving reproducible traversal. Each processed core, all retained children and the next-core cursor commit in one SQLite transaction. Time and work guards pause only between complete cores.

The continuation processed 3,739 further states in two atomic batches: 1,433→3,473 in 600.590 seconds, then 3,473→5,172 in 543.379 seconds. The first pause occurred between complete cores; the second batch stopped because the queue was empty. The original producer and cursor were preserved. The resulting status partition is:

| Status | States |
|---|---:|
| Expanded nonsingular | 2,760 |
| Sealed nonsingular, excluded by rank | 132 |
| Nonsingular at the 31-root boundary, excluded by rank | 210 |
| Sealed singular weight obligations | 1,871 |
| Expanded singular weight obligations | 43 |
| Singular weight obligations at the 31-root boundary | 156 |
| Total | 5,172 |

The acceleration changes the exact arithmetic path. A positive diagonal radical gauge puts the Gram matrix in Q(√5). The white support is a forest. Eliminating one white leaf at a time gives strictly positive pivots in both embeddings and leaves a four-by-four residual on the black nodes. One field-valued congruence therefore supplies both physical and conjugate inertia; its conjugate is the corresponding congruence in the other embedding. Kernel vectors are recovered by reversing the same elimination steps and are checked against the full original matrix. This avoids repeated dense elimination without introducing approximation.

Every finite-corner/no-short contact profile is still enumerated in the complete finite alphabet {2,3,4,5,6,10}. Nonzero products with an inherited kernel, two physical negative directions, conjugate indefiniteness and nonnegative root dependencies are rejected exactly. A short marked path or an entire H4 endpoint corner gives a persistent inherited payment. A final weight rule must preserve these baseline payments. Finite-edge cycle products are automatically integral when they lie in the arithmetic field, since each factor 2cos(π/m) is an algebraic integer.

The saved complete domains contain 51,599 contact vectors: 9,774 accepted necessary children, 21,308 two-negative rejections, 8,301 finite-cycle-field rejections, 11,629 nonzero inherited-relation rejections, 427 conjugate-signature rejections, and 160 entire-H4 payments. Accepted vectors may lead to an already stored state; all 9,774 child maps are literal full-label bijections.

The 31-root boundary follows from the inherited compact dimension bound29. A full quadrilateral star in dimension d has d+2 roots and rank d+1. A nonsingular 31-root core exceeds Lorentz rank30; a singular one is retained as a boundary payment obligation. Below the boundary, every touching external white component can be added from its first contact and along its connected edges. Orthogonal white spectators are positive definite and cannot turn a sealed nonsingular matrix into the required corank-one full star.

Every singular state is exported as a possible stopping case, even if it admits further touching extensions. Its four finite corner forests, marked-pair orbits and spectator budgets are saved. The explicit dictionaries list target dimensions28 and29. For any lower target d, the same catalogue applies: keep rows with at most d+2 nodes, and the spectator rank is d+2 minus the row size. In particular dimension27 uses rows with at most29 nodes and spectator rank29 minus the row size. The full queue nevertheless continues through31 nodes.

This is necessary attachment coverage from the 311 literal nonsingular first extensions, reduced to 199 initial marked roots. Separate certificates cover the other singular-seed and dotted families. The catalogue does not assert actual realization of every retained graph, complete affine-subdiagram absence, payment by a weight rule, all-rank sealing, or an improved dimension theorem.

Generation passed. A terminal structural audit separately checked database integrity and equality with the exported records, all source manifests and counters, root preservation, reachability of every state, absence of pending siblings, and all 9,774 child maps (2,600,547 label comparisons).

The [full independent review](fresh_sqrt5_nonsingular_rank29_fast_20261015_review.md) passed generation and a fresh saved replay without prefetched arithmetic. It reconstructs all 5,172 states, 51,599 complete contact profiles, 10,085 seed and child maps, and 8,280 corner patterns belonging to the 2,070 singular obligations. Its finite-corner inversion and bordered Schur complements use a separate field implementation from this producer's white-leaf elimination. It also checks the terminal database and literal agreement with the unchanged original producer's immutable 2,007-state checkpoint. A second full same-engine arithmetic replay was unnecessary and has not been run.

The frozen primary hashes are:

| Artifact | SHA-256 |
|---|---|
| Source | `21f07de631a9946519b5219100c387f8ad0c3c18fbed7280c4b0469ba23118d3` |
| Terminal JSON | `f6f14c2b226221a0a507c3c2df6179c0aa1d3ec5bba91136b2db46eeab94d22f` |
| Terminal SQLite | `d2ecf4b036032d87d6ccebd940a896c2eacbb6cb079bb671d3112d5c17ce2d53` |

The source retains a full saved-arithmetic replay command; it reads the terminal archive without resuming the exhausted queue:

```sh
python3 scripts/fresh_sqrt5_nonsingular_rank29_fast_20261015.py --verify
```
