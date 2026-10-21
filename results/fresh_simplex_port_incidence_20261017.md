# Maximal-simplex port packing excludes the saved dual incidence distribution

The sparse dual of the [enlarged rank-13 flexible-star model](fresh_sqrt5_rank13_full_flexible_simplex_lp_20261017.md) satisfies every old root-pair capacity, but violates an additional necessary incidence condition at **22 of its 203 positive-mass vertex forests**. The condition concerns which maximal simplex faces can coexist through one vertex. It does not reject any individual arithmetic local star.

The [exact checker](../scripts/fresh_simplex_port_incidence_20261017.py) and [saved certificate](fresh_simplex_port_incidence_20261017.json) reconstruct the literal port blocks of all 346 positive dual rows. They also export 22 valid new charge columns over all 10,185 old simplex constraints. There is no optimizer in this package, no augmented feasibility result, and no new dimension bound.

Let P be a simple compact acute-angled hyperbolic polytope. At a vertex v, label its n incident edges by the n incident facet normals: an edge has the label of the one normal omitted from its defining set. A d-dimensional face through v selects exactly d of these **ports**. If the face is a simplex, its two-port subsets are triangular two-faces through v.

The frozen [maximal-simplex extension theorem](fresh_sqrt5_rank13_triangle_simplex_flags_20261017.md) says that every triangular two-face of P belongs to a unique maximal simplex face. Its proof uses the invertible full triangle-wall star and all its spherical one-root deletions. Consequently, if two distinct maximal simplex faces through v shared two ports, the unique two-face on those two edges would be a triangle in both, contradicting that theorem. Their port sets therefore intersect in at most one element. This assertion is about actual faces; it does not assume that all individually admissible local stars coexist.

Here is a useful consequence. Fix a set U of u ports and an integer k with 2k−u≥2. Any two subsets of U of size at least k intersect in at least two ports. Thus at most one maximal-simplex port set B can satisfy

\[
 B\subseteq U,\qquad |B|\ge k.
\]

This remains true if the indicator is averaged over all images of U under the labelled finite-diagram automorphism group. Its sum over all maximal simplex faces at a vertex is still at most one. The statement has no bound on simplex dimension and no Coxeter-label cutoff. The bounded diagnostic below tests only (u,k)=(4,3),(5,4),(6,4).

For comparison, root-pair capacities alone do not imply this condition. On four ports, assign mass 1/2 to each of the four three-port blocks. Every pair has total load one, but the total block mass is two, whereas the packing condition permits only one. This is why the new condition is not a restatement of a pair-weight budget or an edge-endpoint balance.

To interpret the saved dual, let α_f be the hypothetical count of a maximal simplex star f and μ_F the hypothetical count of vertices of forest type F. The coefficient of a root-pair orbit is the number of its triangular corner occurrences. The old dual constraints bound their average load by μ_F times the orbit multiplicity. For a simplex whose spherical deletion indices form S, its vertex at deletion o has port set S\{o}. The checker reconstructs this entire block, instead of only its pairs.

At each vertex type, average identifications with a canonical labelled forest uniformly over its diagram automorphisms. If a block orbit has total mass m and b literal members, each member has mass m/b. This is a consequence of symmetry averaging, not an additional geometric symmetry assumption. For a fixed U, summing these literal masses gives the expected packing load. Any actual global incidence distribution with the specified star and vertex counts must obey the bound μ_F. The calculation also independently reproduces all 15,834 literal pair loads and checks that the old μ_F is exactly the largest pair load at each of the 203 forests.

The smallest convenient witness uses forest index 390,

\[
 F=2A_1+A_2+3H_3.
\]

Let a,b be the two A1 ports and c,d the two A2 ports, so U={a,b,c,d} is invariant under all automorphisms of this forest. The two positive dual contributions are:

| Old simplex row / literal star | Ports at its relevant vertex or vertices | Total mass |
| --- | --- | ---: |
| 1008 / 1208 | {a,b,c} or {a,b,d} | 48532913 / 10000000000 |
| 8907 / 10619 | {a,c,d} or {b,c,d}, with two vertex occurrences | 46911 / 8000000 |

Every pair of these three-port sets intersects in two ports. The valid inequality is therefore

\[
 \alpha_{1008}+2\alpha_{8907}\le\mu_{390}.
\]

The saved dual has μ390=46911/8000000, while its left side is 107171663/10000000000. The ratio is 107171663/58638750, approximately 1.82766. The second row alone uses the entire allowed vertex mass; the first row cannot be added there. All six pair capacities within U, and all remaining pair capacities of the full forest, nevertheless pass.

The new condition also meets the dominant exceptional part of the old dual. At forest 951, A1+A5+H3+H4, let U be the invariant six-port A1+A5 subset. Any two four-port subsets of U intersect in at least two elements, giving

\[
 \alpha_{599}+\alpha_{957}\le\mu_{951}.
\]

Row 599 is the dominant star 681 identified in the old report. The left side is 300792247/10000000000 and μ951=141551883/5000000000, a ratio of 300792247/283103766>1. The other three highlighted forests—2A1+A4+H3+H4, 2A1+H3+2H4, and A1+A3+H3+H4+I2(5)—pass these particular small-clique tests. No general packing sufficiency is asserted for them.

The checker tests 754,754 literal clique inequalities and retains one strongest tested inequality per violating forest. Holding every saved face multiplier fixed, pair capacities together with these 22 inequalities require at least

\[
 \frac{7706660539}{300000000000}
\]

additional vertex mass beyond the old normalizer. The same fixed face demand divided by this necessary vertex mass is at most

\[
 \frac{2104176880435}{179495552776}\approx11.7227243121<12.
\]

This rejects the old dual masses as a candidate global incidence distribution. It does not construct a realizable distribution after increasing the masses, and it does not rule out a different obstruction to target 12.

The successor columns have a direct geometric meaning. Assign a nonnegative coefficient to one averaged clique indicator at a specified vertex forest. Pay that coefficient at each maximal-simplex corner selected by the indicator. The total contribution to that vertex budget is at most the coefficient, because the clique packing bound is one. The coefficient in a simplex row is the sum of the indicator values at its vertices. For a block orbit, this value is

\[
 \frac{\#\{\text{literal blocks in the orbit contained in }U\}}
      {\#\{\text{literal blocks in the orbit}\}},
\]

provided its size is at least k. Double-counting block–U incidences verifies that this equals averaging over the U orbit. The export gives each column's destination vertex, cost one, exact rational face coefficients, U orbit, and literal dual provenance. It covers every old face, including rows with zero dual multiplier. There are 45 nonzero coefficients across the 22 columns. Summing these new charges with the old pair charges is valid even when the different clique charges overlap: their vertex bounds are added separately.

The unique-maximal-simplex theorem is the inherited geometric premise. The new local packing proof and the exact incidence test are supplied here; no independent second arithmetic implementation is claimed. Generation and fresh saved replay both passed. The remaining task is a separate augmented partial-model solve, retaining every old face and vertex budget. Even a feasible result there would leave the old model's omitted face classes and vertex forests unresolved.
