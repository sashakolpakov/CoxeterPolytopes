# Independently implemented complete T8 traversal

**PASS:** a second complete implementation exhausts all 59,423 mandatory-label search nodes and finds no survivor. This independently confirms the compact-Coxeter exclusion of the literal 18-facet T8 incidence in the [primary report](fresh_h9_t8_sparse_search_20261014.md). It does not exclude the other tree incidences or all compact H9 polytopes.

The [new script](../scripts/fresh_h9_t8_sparse_search_20261014_independent_dfs.py) imports only the Python standard library. Its [certificate](fresh_h9_t8_sparse_search_20261014_independent_dfs.json) records the complete result and source hashes. It uses new immutable-assignment recursion, new propagation and bridge checks, independently generated triangle categories, and separate exact finite-Gram arithmetic. It has **no node cutoff**. The final node and prune totals are compared with the original only after the independent recursion has returned with no survivor.

## Independent finite-Gram arithmetic

The original diagram-classification function is not imported. For each connected nonright-contact component, the new implementation applies the following exact rules to a normalized Gram matrix with diagonal one:

1. A component of order at most two is positive definite for every finite integer label. Thus labels six and the entire category `>=7` remain allowed in isolated rank-two components.
2. A graph with a cycle is not positive definite: the vector equal to one on a shortest cycle has quadratic value at most zero, since every cycle bond has cosine at least one half. Additional acute contacts only decrease the value.
3. A connected component of order at least three containing a label at least six is not positive definite. An adjacent nonright bond gives a three-node path with determinant at most `1-3/4-1/4=0`.
4. Every remaining component is a tree with labels three, four or five. Exact leaf elimination subtracts `cos(pi/m)^2 / pivot` from the parent pivot. Every pivot must be strictly positive. The squared cosines are `1/4`, `1/2` and `(3+sqrt(5))/8`, so all calculations use separate Fraction arithmetic in `Q(sqrt(5))`. Signs are decided by rational comparisons and squares, with no floating-point tolerance.

This is a complete positive-definiteness test for every diagram encountered. The run checks 85,018 distinct finite masks and caches 145 distinct rooted-tree Schur computations. It uses neither a finite Dynkin catalogue nor arithmetic functions from the producer.

## Independent coverage and traversal

The script reconstructs all 462 nine-facet vertex incidences directly from the literal minimal-nonface clutter and compares them with the frozen geometry. It finds ten missing triples and eleven missing pairs. It independently generates all 175 ordered triangle categories using exact Fraction reciprocal sums. Category seven means every finite label at least seven; it is not a fixed angle substitution.

All 16 combinations of the four twin swaps are checked on both the complete vertex set and the minimal-nonface set. The selected ordering leaves 0,3,6,11 lie outside every twin pair, so the four order reductions are simultaneously valid. The 60 disjoint-negative-block clauses are also rebuilt from the clutter, and every cross pair is checked to occur at an actual vertex.

Assignments are immutable tuples. Each recursion creates a new assignment and new filtered triangle domains. It tests newly changed finite contacts on every affected actual vertex; unassigned finite contacts are relaxed to right angles. Bridge pruning uses independently recomputed label-three failures, and the four analytic theorem predicates are applied only when every remaining category satisfies every bound or equality.

The literal incidence's primal realizability and the separately reviewed analytic family theorems are explicit inherited dependencies. Their recursive source hashes are checked. Their polynomial or negative-plane arithmetic is not repeated by this script. The new independent part is the entire discrete traversal, including its finite-Gram arithmetic, geometry extraction, category generation and propagation. The [separate coverage review](fresh_h9_t8_sparse_search_20261014_outside_coverage_review.md) also audits the original DFS logic and its symmetry/category premises.

## Complete result

| Check | Independent result |
|---|---:|
| Search nodes | 59,423 |
| Rejected finite-vertex trials | 238,748 |
| Orthogonal-negative-block prunes | 32,665 |
| Sparse orientation prunes | 3,420 + 3,420 |
| Two-high-placement prunes | 3,648 + 7,296 |
| Surviving assignments | 0 |
| Unexplored frontier | None |
| Cutoff | None |

All counts agree with the primary traversal. The matching search order is a deliberate control; arithmetic and implementation are separate. A compact realization would induce a surviving mandatory-label assignment, so this empty necessary search excludes that incidence without choosing any dotted distances or attempting a numerical metric fit.

```sh
python scripts/fresh_h9_t8_sparse_search_20261014_independent_dfs.py --verify
```

Generation and the final complete saved replay both pass. This new supplement is frozen; no previously frozen source or report is changed.
