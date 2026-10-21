# Independent review of the 815 cap capacities and their strengthenings

The [independent checker](../scripts/fresh_sqrt5_orbit_and_orientation_caps_20261017_review.py)
and [exact result](fresh_sqrt5_orbit_and_orientation_caps_20261017_review.json)
accept the geometric capacity argument and independently reconstruct the
selected arithmetic and incidence data in all three packages:

- [815 orbit and integrality capacities](fresh_sqrt5_orbit_integral_cap_capacity_20261017.md);
- [215 cross-orientation tightenings](fresh_sqrt5_cross_orientation_cap_capacity_20261017.md);
- [36 distinct-profile clique capacities](fresh_sqrt5_cap_profile_cliques_20261017.md).

No correction was identified. The final family has 851 valid selected
columns, of which 46 violate the unchanged diagnostic dual. This review
does not repair that dual, solve a model, or establish a dimension bound.

## What is independently reconstructed

The checker imports only the Python standard library. Its field is
`u+v*phi`, with `phi^2=phi+1`, rather than the producers' square-root basis.
It uses a separate reverse-order Gaussian solver and exact substitution
checks. Where physical core entries require an additional square root,
it constructs a positive tree gauge rooted at the opposite end from the
older helper. Cross-product squares are recovered from that gauge;
the producer's special two-node closed formula is not called.

It independently reconstructs all 47 selected core automorphism groups
by adjacency-preserving backtracking, with no producer colour partition.
It also checks every spherical one-root deletion, a negative determinant,
and an existing positive principal submatrix. These establish the needed
Lorentz signature and that the graph symmetries preserve the finite
vertex labels.

Among the 1,827 isolated-A1 incidence records, it reconstructs the 1,100
with exactly one internal A1 outside the core, then recovers exactly the
815 selected occurrences. For every selected occurrence it checks the
raw core isomorphism, the cap and coface roots, spectator orthogonality,
the complete contact vector, and the full-rank dimension counts.
The base vertices are rebuilt from the literal quadrilateral partitions:
for a deleted cap root b, a vertex is labelled by the other omitted root
p in a positive-definite deletion pair {b,p}. This independently verifies
the `canonical_vertices` convention rather than using it as an untested
list of roots incident to the base.

All simultaneous contact-vector/base-subset orbits, projection norms,
positive-side incidence coefficients, and signed column loads are
recomputed. The 241 old columns are retained coefficient by coefficient.
The counts are 572 capacities one and 243 capacities two, 167 positive
loads, and 19 violated columns before cross-orientation tightening.

The 243 cross-orientation products are then recomputed. Exactly 215
have intersecting base subsets and a nonintegral cross-product square;
the resulting 42 violated columns include 23 new violations. Finally,
all 80 pairs of capacity-one profiles sharing a base type are checked.
The 70 excluded pairs reproduce the incompatibility graph. An independent
recursive maximal-clique algorithm gives 22 pairs, 12 triples and two
quadruples. Their incidence sums and 36 new columns are checked exactly,
including the four additional violated columns.

## Capacity semantics

Fix one actual maximal simplex S, including its full n+1-wall star,
and choose one identification of its Lorentzian core with the canonical
core. This identification is fixed simultaneously for all incidences.
The remaining walls of the star are finite orthogonal spectators. In
the selected prism occurrence the cap b has nonzero contact with exactly
one spectator, an isolated norm-2 A1 root q; appended spectators are
orthogonal to b. Consequently the complete contacts determine

    b = v - sqrt(tau) q/2,    v^2 = 2 - tau/2.

The core is nonsingular, so its contact vector fixes v uniquely. At the
same q, the spanning full star fixes b uniquely. Fixing its base-vertex
subset also fixes the base face. The inherited canonical dotted-prism
incidence theorem then prevents this one incidence from being counted
twice. At two distinct A1 copies q,q', the spectator vectors are
orthogonal, so two caps of the same orientation have product v^2.
The saved finite-angle or arithmetic obstruction forbids that repetition.
Thus each orientation has capacity one. The number of simultaneous
orientations, rather than the order of the automorphism group, is the
correct capacity for their union. No extra spectator-copy multiplicity
is available: the exclusion already compares different copies.

For two orientations with overlapping base-vertex subsets, occurrences
at the same q use the same immediate coface U of S and share a base
vertex. The inherited disjoint-base theorem forbids this. It also covers
the possibility of a repeated prism aggregate: its distinct bases in S
cannot share that vertex. At different q,q', the cap product is just
`v0.v1`; its nonintegral square excludes coexistence in an arithmetic
chamber. Both tests are necessary for the stated capacity-one union,
as is the preceding within-orientation bound.

The clique extension applies exactly the same argument to every pair of
orientations of two different profiles. Their tau values may differ;
the spectator terms still disappear at distinct copies. Each individual
profile already has capacity one, and every pair of profiles in the
selected clique is incompatible. Their total incidence count is therefore
at most one. Source-base lists are checked disjoint before coefficients
are added. Separate columns need not describe disjoint resources: each
is its own universally valid signed inequality.

These arguments use the inherited whole-star and shared-base theorems
for actual maximal simplices and canonical dotted-prism aggregates.
They are not claims about arbitrary all-finite overlapping products.
Omitted positive occurrence types have coefficient zero; omitting them
cannot invalidate the capacity for the selected subfamily.

## Angle and arithmetic normalization

With wall squares two, the normalized length-two cycle is `(b0.b1)^2`.
It is not the unsquared product, nor a quantity missing a factor four.
Algebraic integrality of this cycle is a necessary arithmetic condition;
see the normalization, cyclic products and Vinberg criterion in
[Antolin-Camarena, Maloney and Roeder](https://arxiv.org/pdf/0708.2109).
The cited paper states its arithmeticity application in dimension three;
the necessary cycle condition here is the inherited general arithmetic
reflection-group criterion, not a new dimension-specific sufficiency claim.
Nonintegrality cannot disappear on adjoining a larger number field.

The independent field test uses the ring of integers Z[phi]. For
`z=u+v*phi`, it compares integral u,v with integral trace `2u+v` and
norm `u^2+uv-v^2`. These agree exactly on every tested cycle.

The finite-angle test is also correctly scoped. If a normalized negative
product h lies strictly between -2 and 0 and belongs to Q(sqrt(5)), then
`h=-2cos(pi/m)` has degree at most two. The degree formula gives
`phi(2m)<=4`. A prime divisor p of such an order has p-1 dividing the
totient, hence only primes 2,3,5 occur; their exponents are at most3,1,1.
This leaves m=3,4,5,6. The square roots for m=4,6 are outside Q(sqrt(5)),
so only -1 and -phi remain. Parallel and dotted cases are not rejected
by this finite-angle argument. The arithmetic test may separately reject
a dotted product. No finite alphabet is assumed for an unspecified pair.

## Verification limits

Final generation63104 and fresh saved replay73810 both terminated with
exit0/PASS. Syntax, links, eleven direct manifests and 115 hash references
also pass. A separate conceptual read by the frontier agent accepts the
same-copy and distinct-copy capacity proof under the stated fixed-core
and literal-base hypotheses; that read is not an arithmetic census.
The reconstruction does not import either cap
producer or the older arithmetic review helper. It does inherit the
frozen raw arithmetic/local-realization catalogues and the established
maximal-simplex/prism incidence theorems; those large catalogues are not
re-enumerated here. This is an independent review of the complete selected
family, not a census of all possible cap profiles or of omitted compatible
pairs. No optimizer was invoked.
