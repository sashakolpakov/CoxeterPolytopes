The two new caps pass a separate exact arithmetic replay. They form a coherent 29-wall Coxeter chamber with **66 finite vertices and 455 geodesic edges**, but it still has infinite volume. Both caps cut the previous exact common-bank escape; the saved replacement escape is strictly inside all 29 halfspaces and spacelike.

The [independent checker](../scripts/fresh_h9_joint_frontier_20261013_coupled_independent.py) and [certificate](fresh_h9_joint_frontier_20261013_coupled_independent.json) bind the frozen [two-cap source](../scripts/fresh_h9_makarov_joint_completion_20261013_coupled_exact.py), its [data](fresh_h9_makarov_joint_completion_20261013_coupled_exact.json), and the corrected helper-source manifest. They import no cap-producer arithmetic. They reuse our separate standard-library Fraction polynomial implementation and 160-bit outward rational square-root intervals, and freshly replay the complete [old frontier](fresh_h9_joint_frontier_20261013.md) and [all-edge census](fresh_h9_joint_frontier_20261013_all_edges.md).

The cap normals are positive multiples of vectors in \(K=\mathbb Q(\sqrt5)\). The replay uses those projective vectors with their actual positive norms, checking the normalization identity \(q d=2\). This eliminates the two normalization square roots from the geometry without assuming they are independent of the old coordinate field. Each nonzero new pairing is certified negative, and its square is strictly greater than the product of the two wall norms. Thus the 55 new pairs consist of exactly 15 right angles and 40 strictly dotted pairs, with no finite-angle approximation. All original 351 pair checks are also replayed.

Both new walls are orthogonal to

\[
S=\{0,2,3,5,7,8,9\},
\]

and to each other; every other old wall is strictly dotted to them. Their common center pairings are negative. The new corner is reconstructed independently by projecting that center off the nine incident normals. Its Gram is positive definite, its point is timelike and future-directed, and its only zero wall pairings are at \(S\cup\{27,28\}\). All 65 old finite vertices remain strictly inside both new halfspaces. Acuteness and the strict center also ensure every listed wall is a genuine facet: projecting the center to any wall gives a timelike point strict in all the other inequalities.

The full finite-pair graph has 352 nine-cliques, exactly 66 of which are positive definite. The only new one is the corner above. It has 1,114 eight-cliques, exactly 455 of which are positive definite. The nine added edge supports are obtained by deleting one entry from \(S\cup\{27,28\}\). For every PD8 support the checker tests every existing finite-pair ninth extension by its exact Schur complement; none has zero complement. The resulting complete edge classification is:

| Existing endpoints | Edges |
| --- | ---: |
| Two finite endpoints | 144 |
| One finite endpoint, one open end | 306 |
| No finite endpoints | 5 |

The old 446 supports retain their endpoint types. The new corner forms a separate, isolated component of the finite-edge graph; its other component contains the original 65 vertices. Each of its nine outgoing edges has an independently constructed exact future spacelike continuation, with exactly its eight supporting wall pairings zero and all other pairings negative. These points are saved in the certificate. The lack of ideal extensions is a statement about genuine ideal vertices, not an assertion that the cone lacks an open ideal-boundary region.

The strict escape for the whole chamber is

\[
(15272584,19332841,13954294,8575997,3198152,706462,
-27462556,-33326515,-22217679,-11108841),
\]

with exact positive norm

\[
2712511767010814-1210493722603484\sqrt5.
\]

All 29 inequalities and its future orientation are checked. Therefore this is positive local construction progress, not a new finite-volume or compact H9 example.

There is also a conditional arithmetic obstruction. For old wall 1 and new wall 27 the scale-invariant Cartan two-cycle is exactly

\[
\frac{4\langle w_1,v\rangle^2}{w_1^2v^2}
=\frac{1349500093}{480261426}
+\frac{94041207}{160087142}\sqrt5.
\]

Its trace over \(\mathbb Q\) is \(1349500093/240130713\), which is not an integer, so this number is not an algebraic integer. Vinberg's necessary cycle-integrality condition therefore implies that **any finite-volume Coxeter completion retaining this pair of facets is nonarithmetic**. The condition is recalled in [Bogachev–Kolpakov, Theorem 2.8](https://arxiv.org/pdf/2002.11445). This does not assert finite volume or classify the current infinite-volume reflection group as an arithmetic lattice.

Generation and the saved-result replay both passed:

```sh
python scripts/fresh_h9_joint_frontier_20261013_coupled_independent.py
python scripts/fresh_h9_joint_frontier_20261013_coupled_independent.py --verify
```

This audit independently recomputes the exact arithmetic and geometry of the saved chamber. It does not claim completeness of a search for further walls or first-mirror completeness in a full arithmetic reflection group.
