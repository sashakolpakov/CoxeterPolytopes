# Complete disjoint-base packing inequalities on at most five vertices

The abstract packing domain for prism bases inside a simplex of dimension
at most4 now has a complete exact description. There are respectively
1,1,8,294 nontrivial inequalities for simplices with2,3,4,5 vertices,
together with nonnegativity. On five vertices this is a26-dimensional
polytope with52 integer packings. All294 inequalities are verified as
facets. This class strictly strengthens every subset-floor capacity.

The [standalone checker](../scripts/fresh_simplex_base_packing_cone_20261017.py)
and [data](fresh_simplex_base_packing_cone_20261017.json) use exact
integers/Fractions, no project imports and no optimizer. This is the full
abstract disjointness domain, not a claim that every packing is geometrically
realizable or that it improves an ambient dimension bound.

## The packing problem and its exact dual

Let V have m vertices, 2<=m<=5. A base is a subset B of size at least2;
a packing is a family of pairwise disjoint bases. Single vertices are
unused. Nonnegative demands d_B are feasible at capacity c if they are
the marginals of a nonnegative mixture of packings of total mass at most c.
Equivalently one may require coverage at least d_B: deleting individual
bases from packings decreases any excess marginal without increasing mass.

Every packing has at most two bases. Its complete list is therefore empty,
all singleton-base packings, and all pairs of disjoint bases. The lists
have2,5,15,52 members for m=2,3,4,5. The minimum covering mass has the dual

    maximize sum_B d_B y_B,
    0 <= y_B <= 1,
    y_B+y_C <= 1 whenever B and C are disjoint.

Finite-dimensional linear-programming duality applies: the primal is
feasible using singleton-base packings and has a finite optimum. Thus
all the dual inequalities are necessary and sufficient, including at
zero capacity. This statement concerns arbitrary rational demands and
does not rely on a numerical optimization.

For m<=3 all bases intersect, so the sole inequality is sum d_B<=c.
For m=4 the compatibility graph consists of three disjoint edges joining
complementary pairs, and five isolated larger bases. Its covering mass is

    sum_{|B|>=3} d_B + sum_{complementary pairs {B,C}} max(d_B,d_C).

Choosing one term from each maximum gives the complete eight inequalities.

## Five vertices: a Petersen-graph reduction

There are ten pairs, ten triples, five four-sets and one five-set. The
compatibility graph on these26 bases has the Petersen graph on the ten
pairs (two pairs adjacent precisely when disjoint), one pendant triple
at each pair (its complement), and six isolated larger bases.

For nonnegative demands a dual optimum can always be chosen with

    y_{V\p}=1-y_p for every pair p,
    y_B=1 for every |B|>=4.

Consequently its value is

    sum_{|B|>=3} d_B
      + max sum_{pairs p} (d_p-d_{V\p}) x_p,

where 0<=x_p<=1 and x_p+x_q<=1 for disjoint pairs p,q. The objective
coefficients may have either sign; no monotonicity assumption is made.

Every vertex of this ten-variable polytope is half-integral. To prove it,
consider the tight-edge graph on coordinates strictly between0 and1. A
bipartite connected component permits a small alternating positive/negative
perturbation, preserving its tight equalities, all strict inequalities and
the bounds. That contradicts extremality. Each component therefore has an
odd cycle; the equalities x_p+x_q=1 then force every coordinate in that
component to be1/2. Conversely, coordinates fixed at0 or1 and fractional
components with odd cycles have a full-rank active system and are vertices.

The checker enumerates all3^10 assignments from{0,1/2,1}. Exactly1,959 are
feasible, and the tight-component test retains exactly294 vertices:

| Number of half coordinates | Number of vertices |
|---|---:|
|0|76|
|5|12|
|6|60|
|7|90|
|8|45|
|9|10|
|10|1|

Every retained active system independently has rank10. Lifting these
vertices back to the26 base coordinates gives the complete inequalities.
The checker verifies every inequality on all52 packings and verifies
affine rank25 among its tight packings. Nonnegativity plus these294 facets
is therefore a complete description, not a sampled or partial cut family.
Coefficients and capacity are stored doubled, so all saved entries are
integers0,1,2 and the saved capacity is2.

## A strict separation from all subset-floor cuts

On V={1,2,3,4,5}, assign mass1/3 to each of the four bases

    {1,2}, {1,3}, {2,3,4}, {1,4,5}.

They are pairwise intersecting, so their total mass must be at most1 in
any packing mixture of capacity1. Their total is4/3. Nevertheless, for
every nonempty T subset V and every 1<=r<=|T|,

    sum_B floor(|B intersect T|/r) d_B <= floor(|T|/r).

All80 inequalities are checked exactly, including r=1 and thus every
individual vertex capacity. A full-cone facet detects the violation4/3>1.
This supplies a new local obstruction beyond the previous labelled-floor
family. It is an abstract fractional witness, not an actual polytope.

## Application scope

For a fixed maximal simplex S and immediate coface U, the reviewed
geometric argument makes the selected prism base vertex sets disjoint.
The new inequalities therefore apply to the simultaneously averaged
coface slices, with their original capacities and provenance coefficients.
They add no interaction between different cofaces and impose no further
geometric compatibility between disjoint bases. Those restrictions remain
available if this entire disjointness cone passes at the arithmetic witness.

Generation and a fresh saved replay pass. The
[independent checker](review_simplex_base_packing_cone_20261017.md) rebuilds
all packings from set partitions and finds the294 vertices by exact
active-matrix ranks, without using the producer's odd-cycle predicate.
It independently checks all304 facets and the80 subset-floor witness
inequalities. The [successor diagnostic](fresh_sqrt5_complete_coface_base_cones_20261017.md)
applies the complete domain to the arithmetic witness. No dimension
conclusion is asserted by this abstract package.
