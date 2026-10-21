# Cross-orientation arithmetic cycles tighten 215 cap capacities

The [exact certificate](fresh_sqrt5_cross_orientation_cap_capacity_20261017.json)
and [producer](../scripts/fresh_sqrt5_cross_orientation_cap_capacity_20261017.py)
retain all [815 orbit capacities](fresh_sqrt5_orbit_integral_cap_capacity_20261017.md)
and lower 215 of their capacities from2 to1. There are now787 capacity-one
columns and28 capacity-two columns. At the unchanged full-model dual,
42 capacities are violated, including23 new violations. This excludes
that witness more strongly; it is not a new optimum or dimension theorem.

## Why two orientations cannot coexist

Keep the inherited geometry: the full simplex star spans the ambient
space, its nonsingular Lorentzian core is orthogonal to its finite
spectators, and each selected cap has exactly one isolated A1 spectator
q with nonzero product. Its norm-2 wall vector is

    b_i = v_i - sqrt(tau) q/2.

The previous theorem bounds each individual core-contact/base-subset
orientation by one. For an orbit with two orientations, write their
canonical base-vertex subsets as B0,B1. In every strengthened column
these subsets intersect. At the same q, both incidences would use the
same simplex-vertex/coface flag, excluded by the inherited disjoint-base
capacity. They therefore cannot coexist at the same q.

At distinct isolated A1 copies q,q', the spectator vectors are orthogonal
and both are orthogonal to the core. The mutual wall product is thus

    (b_0,b_1) = (v_0,v_1) = p_0^T G_core^(-1) p_1,

where p_i is the complete core-contact vector. Arithmeticity requires
the square of this product to be an algebraic integer. In each of the
215 strengthened columns it is not. Hence different orientations also
cannot coexist at distinct cofaces, and their total capacity is one.
No finite-angle restriction is needed for this cross-orientation step.

This argument uses the literal simultaneous orbit of contacts and base
subsets, with a fixed identification of the actual core. It does not
choose incompatible canonical identifications independently for the
two incidences. The old within-orientation theorem remains necessary;
the cross-orientation check alone would not prove capacity one.

## Exact projections and checks

All243 two-orientation columns have only labels2 and3 in their cap/core
contact vectors. For cores with labels3 and5, the checker inverts the
physical norm-2 Gram matrix over Q(sqrt(5)), checks the matrix product
against the identity, and checks both projection norms against the
inherited value h=2−tau/2. It computes their cross product exactly.

The other cores have two nodes with dotted squared product d>4. Their
Gram inverse is

    [[2,sqrt(d)],[sqrt(d),2]] / (4-d).

The two orientations contact opposite nodes by a bond3. Their projection
norms are2/(4-d), while their cross-product square is d/(d−4)^2.
This gives an exact calculation in Q(sqrt(5)) without adjoining sqrt(d).

For each square z=a+b sqrt(5), the checker independently compares the
integral-basis test (a−b,2b integral) with the trace/norm test
(2a,a^2−5b^2 integral). The data record both rational coefficients,
trace, norm, the intersecting base vertices, and the projection method.
The remaining28 capacity-two columns are left unchanged. No assertion
is made that this is the strongest possible capacity for them.

Generation and fresh saved replay pass. All815 inherited positive-side
coefficients and source-base lists are retained. Of the215 tightened
columns,51 have positive load at the diagnostic witness and24 are now
violated; one was already violated at capacity two. Together with the
19 inherited violations this gives42 distinct violated columns.

The arithmetic helper comes from the prior independent cap review and
is explicitly source-bound. This producer is not itself an independent
implementation of that helper or of the inherited815-cap census.
Repairing the witness requires updating its face masses and covers;
this package performs neither that repair nor an optimization.
