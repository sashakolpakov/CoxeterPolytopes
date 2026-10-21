# Finite-vertex caps, mixed arithmetic roots, and constrained completion

**No new compact or finite-volume Coxeter polytope in a new dimension has
been verified.** This round changes how candidate facets are constructed
and selected. In particular, it creates new finite vertices, permits
previously excluded dotted relations, and searches compatible completions
with explicit escape constraints.

## H9 caps that create finite vertices

The [changed-core continuation](fresh_h9_gap_completion_20261013.md) first
uses a different retained quartic wall as an anchor, adds two compatible
cuts, replaces one older wall by its seed-symmetric image, and reaches a
48-wall chamber. Exact escapes remain.

The next construction changes the geometric requirement. Choose eight
retained walls with positive-definite Gram matrix and construct a norm-two
normal orthogonal to all eight. The pairing with polar wall 348 is now
allowed to be dotted. Thirty exact caps satisfy eight right angles and
strictly dotted relations to every other retained wall. Each has an
explicit actual finite vertex. Compatible selections reach 70 facets;
the numerical finite-vertex count rises from 15 to 37, while exact
spacelike escapes still survive.

A further coupling imposes a finite angle between two of these new cap
families. Fifteen exact additional caps give two coherent 73-wall
configurations. Seven old walls and two new caps have a positive-definite
rank-nine Gram matrix and an explicit actual joint finite vertex. Five
such joint vertices are certified across the two configurations.

The final producer-free replay checks 5,256 selected pairs, 1,514 affine
frontier/inertia entries, five joint finite vertices, and 24 strict
spacelike escapes. Thus the new facets and finite vertices are real, but
both chambers still have infinite volume. Their leading remaining escape
has seed support `[2,3,6,8,9]`. No general obstruction to further
completion is claimed.

The [separate geometry audit](fresh_h9_gap_completion_independent_geometry_20261013.md)
reconstructs the first 73-wall case with independent quadratic-tower
arithmetic. It verifies the common interior, an exact strict escape, and
the actual joint vertex involving a new pi/3 pair, including every one of
the 73 inequalities. This complements the full pair and affine replay.

## Arithmetic mirrors that close an actual edge

The [conductor-sixteen construction](fresh_h9_quartic_mixed_20261013.md)
works over `Q(T)`, where `T=2cos(pi/8)`, in the lattice
`<-(T-1)> + F4 + F4 + <1>`. A complete two-block root census through
height 1,000 extends the earlier 22-wall cone to 62 walls. Targeted mixed
and neighboring-weight slices through height 100,000 give a mutually
compatible 170-wall H9 cone and a 68-wall H5 face. Separate exact
enumerations verify the declared slice censuses.

The next step solves an endpoint equation directly. Fixing the finite
couplings at an edge of type `A2+A1+A1` reduces the norm-two equation to
`6t(t+1)=(T-1)x^2`. The substitutions `x=3s`, `t=(y-1)/2` give
`y^2-6(T-1)s^2=1`. An explicit integral solution produces four compatible
symmetry-related mirrors, reaching **174 H9 walls and 70 H5 walls**.
Their height is approximately `2.0269e22`; this is an algebraically
constructed root, not a uniform search to that height.

The new wall closes the previously escaping H5 edge at a point with
square `-(T-1)/y^2 < 0`. All 174 inequalities hold there, and the nine
incident independent walls have positive definite diagram
`F4+A2+A2+A1`. Thus it is an actual finite vertex. The independent
Fraction replay checks the Pell identity, all 17,466 H9/H5 wall pairs,
the endpoint inequalities, and the complete 30-profile finite-label
endpoint census.

Both final cones nevertheless retain exact strict spacelike escapes.
The construction does not yet prove that the full reflection group has
a compact chamber. A separate all-height local obstruction rules out
roots confined to one specified binary plane; it does not rule out the
full host.

## Exact reconstruction of the larger split-seed bank

The previous larger scout had 3,040 numerical cap profiles. Their
[exact reconstruction](../data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json)
and [independent audit](fresh_h9_split_leaf_3040_independent_20261013.md)
now verify every normal, all 48,640 seed pairings, and the common strict
timelike interior.

The independent audit corrects the census to **3,034 distinct mirrors**.
Six pairs of numerical profiles had split the repeated root `d=phi^2`
into two nearby approximations. Exact reconstruction identifies the pairs.
The new selection search evaluates the exact reconstructed normals, not
those displaced numerical roots.

The [constraint-guided search](../scripts/fresh_h9_split_leaf_z3_20261013.py)
uses Boolean variables for candidate facets. Each geometric escape
requires a selected facet that cuts it; every discovered incompatible
pair is forbidden. A compatible selection is extended and tested
geometrically, and its surviving escapes become further constraints.
This replaces repeated independent greedy choices by one accumulating
completion problem.

The [exact escape-clause checker](../scripts/fresh_h9_split_leaf_escape_constraints_20261013.py)
replaces numerical points by rational seed coefficients. It proves positive
Lorentz norm and strict seed inequalities. Exact integer interval bounds
provide conservative cutter sets: every omitted cap has nonpositive
pairing, and an uncertain sign stays allowed. A finite-volume completion
must satisfy every resulting clause. These necessary clauses alone do
not prove that completion is impossible.

The [independent escape audit](fresh_h9_split_leaf_escape_constraints_independent_20261013.md)
checks all 3,164,640 cap pairings for the initial 1,041 clauses. In this
archive the conservative possible-cutter sets equal the exact cutter
sets, with no omitted or extra caps.

The [rigorous graph](fresh_h9_split_leaf_3034_degree_graph_20261013.md)
certifies 1,427,166 incompatible cap pairs. Every other pair remains
allowed, including any unresolved interval. A field-degree argument
limits possible finite labels to fourteen values, without a numerical
label cutoff. The separate replay verifies every rejected pair with
general interval arithmetic. Thirty rounds of the
[certified-constraint search](../scripts/fresh_h9_split_leaf_certified_sat_20261013.py)
accumulate 1,761 exact escape clauses; these conditions alone remain
satisfiable. Satisfiability of necessary conditions is not a finite-volume
certificate.

### Four mandatory edge endpoints close the whole bank

A stronger geometric condition makes the problem small. Each orthogonal
eight-wall seed group has a Lorentzian two-dimensional orthogonal
complement, defining an actual edge in every acute-compatible completion.
At either ideal endpoint, some added cap must have nonnegative pairing:
otherwise a slight outward movement, followed by an interior perturbation,
gives a strict spacelike escape. Zero pairings are retained, since an ideal
endpoint may contribute to finite volume.

The [exact endpoint certificate](fresh_h9_split_leaf_edge_endpoints_20261013.md)
checks all 12,160 pairings. Each of the four endpoint lists has fourteen
caps. Exhausting their compatible choices leaves exactly **36 pairs**:
one entry from `0,...,5` and one from `6,...,11`. Each pair admits only 28
possibly compatible cap entries; their union has 332 entries.

All these entries fail to cut the same small integer-coefficient vector,
written in the ordered sixteen seed normals:

```
w = (8,8,8,1,-3,-3,-3,-3, -3,-3,-3,5,5,5,5,5).
w^2 = (7242-2800*phi)/11 > 0,
(w,S) = -116-232*phi < 0.
```

Every pairing with the sixteen seed walls and all 332 possible remaining
caps is strictly negative. The [exact exhaustive cover check](fresh_h9_split_leaf_endpoint_cover_20261013.md)
and this single escape therefore exclude **every finite-volume completion from
the entire frozen bank of 3,034 distinct mirrors**. The final proof needs
neither a SAT solver nor the 1,761-clause search history. It does not
exclude additional mirrors outside the bank or arbitrary completions of
the seed.

The [independent endpoint-and-cover replay](fresh_h9_split_leaf_endpoint_independent_20261013.md)
reverses the order of the quadratic extensions in the sign calculation,
checks all 12,160 endpoint pairings, independently exhausts all 38,416
endpoint choices, and verifies the same strict universal escape.

## Other seed geometries tested numerically

Two additional discovery searches change the initial Gram matrix rather
than only its cap bank. Both retain exact symbolic parameter formulas,
but their parameter-root selection, compatibility checks, and geometry
are numerical. They are not exact exclusion certificates.

The [multi-block scout](../scripts/fresh_h9_multiblock_seed_20261013.py)
uses two orthogonal sets of `m` roots with cross Gram `-2I-K`, where
`K=Z(aJ-bI)Z^T` and `Z` records blocks of sizes `(2,2,1,1)` or `(3,2,2)`.
The sum/difference spaces have total Lorentz rank ten. The norm-two cap
condition is a quadratic in `a` for fixed `b`. Over the declared six
choices of `b`, the numerical equations give 2,055 parameter points;
132 points prioritized by cap count receive 2,315 greedy compatible-set
tests. Every tested set has a numerical spacelike vertex. The
[archive](../data/enumerations/fresh_h9_multiblock_seed_20261013_scout.json)
records the precise restrictions.

The [ten-wall exchange scout](../scripts/fresh_h9_exchange_seed_20261013.py)
uses two orthogonal sets of five independent roots, with cross Gram
`-(2-b)I-aJ`, where `a=2cos(pi/m)`. It solves quartic norm equations in
`b` for all 1,287 permutation orbits of profiles with magnitudes
`0,1,phi`. Nine declared labels give 5,429 numerical parameter points;
144 points prioritized by cap count receive 4,315 compatible-set tests.
All tested sets still have numerical spacelike vertices. See the
[archive](../data/enumerations/fresh_h9_exchange_seed_20261013_scout.json).
These counts do not exhaust all parameters, all compatible subsets, or
other leaf-angle choices.

## Simultaneous H20 section repairs

The [adapted-section construction](fresh_h20_adapted_sections_20261013.md)
tests four cusp pencils and all eight subsets of three mutually
orthogonal additional doubling seams. The ambient convex unions reach
64 Borcherds chambers and 5,577 true ambient boundary facets.

An exact actual-ridge argument covers **every real parameter in all 32
specified families**. Retaining wall 23 forces one possible Coxeter angle
parameter, yielding four explicit normals over `Q(sqrt(6))`; each has a
separate actual obtuse ridge. Doubling wall 23 instead produces an actual
obtuse ridge throughout the admissible interval. The other allowed seams
preserve these witnessed ridges. An independent standard-library replay
checks all ambient boundaries and the continuous argument.

This rules out the declared joint-repair families. It does not exclude
other cusp directions, further seams, other gluings, or a Coxeter polytope
in H20.
