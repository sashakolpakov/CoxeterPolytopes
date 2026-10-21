# Coupled five-cube cuts and four wedges

The explicit nonproduct simple nine-polytope constructed here has **18 facets and 512 vertices**, but admits no compact hyperbolic Coxeter realization. The final necessary-label search is exhaustive for this one incidence: it terminates after 44,035 nodes, below its 100,000-node ceiling, with no surviving assignment. The nine disjoint-facet distances remain arbitrary continuous variables throughout. This is not an exclusion of other cube cuts, wedge choices, or eighteen-facet polytopes.

## An explicit incidence

Start with the five-cube `−1 ≤ x_i ≤ 1`, and successively add

\[
x_0+x_1\le1,\quad x_1+x_2\le3/2,\quad
x_2+x_3\le7/4,\quad x_3+x_4\le15/8.
\]

Each cut truncates the current codimension-two face supported by its two positive coordinate facets. Its depth is half the smallest positive old-vertex slack outside that face. Exact edge intersections give a simple five-polytope with 14 facets and 70 vertices. Replace each of the four negative inequalities with index `i=0,1,2,3` by

\[
-x_i+z_i\le1,\qquad -x_i-z_i\le1.
\]

These are ordinary wedges at four mutually meeting facets. The [primary geometry](fresh_h9_coupled_cube_wedges_20261014.json) stores every rational inequality and vertex. The resulting nine-polytope has 11 minimal nonfaces of order three and nine of order two, with no higher minimal nonfaces. Its missing-face hypergraph is connected, so the incidence is not a nontrivial Cartesian product. Unlike the earlier edge-cut three-cube multiwedges, it has no mandatory Lannér block of order four or five.

The [separate primal checker](../scripts/fresh_h9_coupled_cube_wedges_20261014_independent.py) enumerates all 2,002 five-row supports of the base, reconstructs the wedge vertices, checks their nine independent active equations and all 9,216 inequalities, and recomputes the entire downward face complex. This verifies realizability, simplicity, the 512 vertices, and the missing-face list without importing the producer.

## From a local survivor to the global signature obstruction

The initial breadth join reached its explicit 20,001-state stopping threshold and was left as an unresolved intermediate archive. An MRV depth-first search then found a complete mandatory assignment in 12 nodes: every missing triangle has labels `(3,3,4)`. With all other meeting pairs relaxed to right angles, all 512 vertex Grams are positive definite. A separate exact `Q(√2)` LDL replay verifies these triangles and vertices. This was a necessary-label survivor, not a Lorentzian metric.

For that assignment, 111 of the 118 omitted meeting contacts are forced right: introducing even the smallest nonright label three makes a literal actual-vertex Gram non-positive-definite. The [separate contact replay](fresh_h9_coupled_cube_wedges_20261014_contact_independent.json) checks all 111 witnesses by exact LDL. The [short metric obstruction](fresh_h9_coupled_cube_wedges_20261014_first_metric_gate.md) needs only nine rational witnesses: triangles `[0,1,2]` and `[9,10,11]` are forced orthogonal. Each has determinant `−√2/4`, so their union has two negative directions, impossible for a hyperbolic Gram matrix.

The [final DFS](../scripts/fresh_h9_coupled_cube_wedges_20261014_bridge_csp.py) applies this reasoning during the search, rather than collecting further local survivors. The incidence contains 73 induced pairs of disjoint missing faces: 12 of sizes `(2,2)`, 38 of `(2,3)`, and 23 of `(3,3)`. Each block must have a negative direction. Thus some cross contact must be nonright. If every such contact is either already right or forced right by an actual-vertex test, the partial assignment is excluded. No value for a dotted distance is needed.

The [terminal certificate](fresh_h9_coupled_cube_wedges_20261014_bridge_csp.json) records:

| Check | Count |
|---|---:|
| Complete DFS nodes | 44,035 |
| Partial assignments excluded by two orthogonal negative blocks | 38,812 |
| Triangle extensions excluded by spherical vertex tests | 134,053 |
| Unexplored frontier | None |
| Surviving complete assignments | 0 |

All four swaps of the two facets in a wedge block are explicitly checked incidence automorphisms before using them to order one pair of labels. The search includes every finite integer label: `7` is a category for all labels at least seven. It is valid for Lannér triangles because the only spherical triple involving such a label has its other two labels equal to two. In a larger spherical connected component, a label at least seven is impossible: an adjacent edge has label at least three, and the corresponding path determinant is strictly less than `1−3/4−1/4=0`. Parallel pairs are not included in this compact-polytope category; the nine missing pairs are strictly dotted.

## Independent arithmetic and scope

The [arithmetic replay](fresh_h9_coupled_cube_wedges_20261014_bridge_independent.json) reproduces the whole DFS and its branch digest while replacing the inherited finite-diagram classifier. It uses a nonpositive cycle witness or exact Schur elimination of weighted trees over `Q(√5)`. Squared edge cosines for labels three through six are `1/4`, `1/2`, `(3+√5)/8`, and `3/4`; 102 distinct rooted Schur subtrees occur. Labels at least seven use the preceding strict path bound. The independent predicate is called 81,909 times.

The traversal itself is reused and source-reviewed, not independently reimplemented. The primal geometry has a separate exhaustive implementation. No numerical metric fit or numerical failure is part of the exclusion. The intermediate breadth-limit archive and the initial mandatory survivor are superseded by the complete bridge certificate, not silently relabeled as complete searches.

Replays:

```sh
python scripts/fresh_h9_coupled_cube_wedges_20261014_independent.py --verify
python scripts/fresh_h9_coupled_cube_wedges_20261014_contact_independent.py --verify
python scripts/fresh_h9_coupled_cube_wedges_20261014_bridge_independent.py --verify
```

The next distinct proposal is to replace the positive-cut path by another tree on the five cube directions and vary the negative facet left unwedged. Those incidences require their own rational realization and complete label tests; this path result does not exclude them.
