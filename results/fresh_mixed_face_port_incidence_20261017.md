# One shared incidence budget for simplex and quadrilateral faces

A vertex has a useful additional constraint linking the triangle and
quadrilateral calculations: **for any three incident edge directions, at most
one face aggregate can occupy two or more of their three two-face positions**.
This supplies a shared charge that can pay both simplex and quadrilateral
aggregate rows. The earlier simplex-only charges cannot pay a quadrilateral.

The [exact checker](../scripts/fresh_mixed_face_port_incidence_20261017.py)
and [certificate](fresh_mixed_face_port_incidence_20261017.json) reconstruct
these indicators on the saved simplex and quadrilateral stars. No optimization
or new dimension exclusion is asserted here.

## Geometric statement

Let P be a simple compact acute hyperbolic polytope, and v a vertex. Label
its n incident edges by their omitted facet normals, as in the
[simplex port theorem](fresh_simplex_port_incidence_20261017.md). Each unordered
pair of ports determines exactly one actual two-face through v.

Partition all triangular two-faces by their
[unique maximal simplex extension](fresh_sqrt5_rank13_triangle_simplex_flags_20261017.md).
Partition all quadrilateral two-faces by their
[complete incident-wall star](fresh_quad_product_aggregation_20261017.md).
These are partitions of actual faces; overlapping maximal product faces are
not separate groups. At v, associate to each aggregate A the set E_A of port
pairs corresponding to its incident triangular or quadrilateral two-faces.
An aggregate with no incident face contributes an empty set.

For distinct aggregates A and B, E_A and E_B are disjoint. Otherwise the
unique two-face determined by their common port pair would belong to two
groups. A triangle and a quadrilateral cannot be that same two-face either.

Fix a triple U of ports. Define

\[
 c_U(A,v)=\mathbf 1\left\{
       |E_A\cap\tbinom U2|\ge2\right\}.
\]

Two subsets of the three-element set of pairs, each of size at least two,
intersect. Therefore

\[
 \boxed{\sum_{A\ni v} c_U(A,v)\le1.}
\]

This is a necessary incidence condition, without any sufficiency or global
realization assertion. It applies to simplex and quad aggregates together.
More generally, the same proof works for any fixed set D of port pairs and
any threshold k with 2k>|D|, using the indicator |E_A intersect D|>=k.

At a fixed spherical vertex forest F, average c_U over the full orbit of U
under the labelled diagram automorphism group. The average still has total
at most one. A nonnegative coefficient beta for this averaged indicator costs
at most beta in F's vertex budget. Its coefficient in an aggregate row is
the sum of the averaged indicator over that aggregate's actual vertices.
This defines a valid additional LP column, in the same way as the previous
simplex packing columns. All geometric hypotheses used by the two aggregation
theorems remain in force.

## It adds information beyond pair capacities

On three ports a,b,c, consider three hypothetical aggregate footprints

`{ab,ac}`, `{ab,bc}`, `{ac,bc}`.

Give each mass 1/2. Every pair then has load one, so all three individual
pair capacities pass. The new indicator sum is 3/2, violating its bound one.
The checker exhausts all disjoint families of nonempty subsets of the three
pairs and verifies that their indicator sum is at most one. This is a
combinatorial relaxation witness. It does not assert that these three
footprints occur in the arithmetic catalogue at a common forest, or that
they are simultaneously geometrically realizable.

## Exact features on the saved face catalogues

For a maximal d-simplex, the corner at a deleted root o has port set
S minus {o}, where S is the full spherical-deletion set. Its footprint is the
complete graph on those d ports. Thus its selected triples are exactly the
three-element subsets of that port set. Across all corners it contributes
(d+1) binomial(d,3) literal indicators.

For a quadrilateral full star, the positive-definite two-root deletions form
the inherited bipartite graph. A K2,2 with sides {i,a} and {j,b} contributes
the port pair {a,b} at vertex omitted {i,j}. The checker builds these sets
directly and rejects any duplicated pair at a vertex. The full aggregate
footprint is bipartite; a selected triple therefore contains exactly two
occupied pairs.

For a dotted star whose unique prism has type Delta1 times Delta(k-1), the
footprint at each of its 2k vertices is a star with k-1 edges. It has
binomial(k-1,2) selected triples. The exact census gives:

| Saved domain | Aggregate vertex occurrences | Literal selected triples |
| --- | ---: | ---: |
| 12,420 simplex stars | 40,683 | 22,668 |
| 1,399 marked dotted stars | 6,640 | 4,224 |

The dotted rows contain 10,128 quadrilateral-corner occurrences. Their
source records are deliberately not quotiented by full-star isomorphism;
these totals are feature-census counts, not global polytope incidence counts.
Orthogonal finite spectators do not change the local occupied pairs. The
coefficient at a full forest still requires the correct automorphism average,
which can change when spectators are added.

The all-finite overlap example from the aggregation theorem is also checked:
its two overlapping triangular prisms contain five distinct quadrilaterals
and exactly 20 quad-corner pairs. We use the whole aggregate, so the shared
quadrilateral is counted once. Its literal footprint features are saved.

Simple triangular faces and individual quadrilaterals occupy only one pair
at each vertex and receive zero from this charge. In particular, this new
lemma does not by itself repair the zero-payment single quadrilateral found
by the separate coverage test. Pair weights and complete face constraints
are still necessary.

**FINAL/FROZEN:** generation and a fresh saved replay passed. A separate
agent's narrow geometric review passed, checking aggregate uniqueness,
one footprint per aggregate at each vertex, and the shared automorphism-averaged
budget. This is internal review. No separate arithmetic reconstruction, complete face
coverage, feasible augmented rule or dimension improvement is claimed.
