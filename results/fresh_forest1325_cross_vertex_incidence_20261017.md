# A forced vertex strengthens the forest1325 corner budget

This package concerns the current compact arithmetic model in hyperbolic
13-space with fixed ground field Q(sqrt5); it is not a complete condition
for unrestricted finite-volume Coxeter polytopes.

The two represented corner types at forest1325 force a second finite
vertex. Combining that incidence with the complete current corner graph
gives the new necessary inequality

\[
3\alpha_{45769}+3\alpha_{26918}+8\alpha_{26362}+4\alpha_{26448}
\ \le\ 4\mu_{1325}+3\mu_{6537}.
\tag{1}
\]

Here1325 is `3 A1 + 2 H3 + H4`, while6537 is `3 H3 + H4`.
The coefficients count actual whole-face corner occurrences, not diagram
orientations. Inequality(1) rejects both current exact fixed witnesses.
Increasing only their1325 vertex allowance gives ratios
**11.999087132529974** for the full witness and
**12.054542455466887** for the named248-Q restriction. The latter remains
above12. This is one additional cross-vertex constraint and an exact
monotone repair, not a global optimization or a dimension theorem.

The [source](../scripts/fresh_forest1325_cross_vertex_incidence_20261017.py)
and [exact data](fresh_forest1325_cross_vertex_incidence_20261017.json)
bind the [current actual baseline](fresh_2599_raw110_raw24_pair_geometry_20261017.md),
the current ownership-preserving graph, and the
[native cyclic-cover reconstruction](fresh_restricted_290_cyclic_cover_20261017.md).
No pair outcome, face mass, or frozen predecessor changes.

## The complete physical configuration

At the representative source vertex, the active incident normals are
`qB,qA,w,h5,hcenter,h3`, consisting of three A1 components and one H3
component with edges `h5--5--hcenter--3--h3`. The other seven incident
normals form the passive `H3+H4` spectator forest. The three external
normals are B,A,C, in that order. The only additional nonzero products,
with all native normals of square2, are

\[
\begin{aligned}
q_B\!\cdot B&=-(1+\sqrt5), &q_A\!\cdot A&=-\sqrt{(7+\sqrt5)/2},\\
w\!\cdot A&=-1, &w\!\cdot C&=-1, &h_3\!\cdot C&=-1,\\
A\!\cdot B&=-(1+\sqrt5)/2, &B\!\cdot C&=-1.
\end{aligned}
\]

In particular A and C are orthogonal. These entries are reconstructed
from the complete native contact covectors at one finite vertex, using
the unique positive timelike height of each external outward normal.
The source repeats this frame's native omitted-product, cross-product,
and cycle checks; it does not select relative signs independently.

The raw25 star12 has external walls{B,C}, face26918 and provenance541;
the raw23 star46 has external walls{A,C}, face26448 and provenance490.
The shared C is an actual identical wall. Scaling only qA by the positive
square root of `(7+sqrt5)/2` puts the entire active Gram over Q(sqrt5).
Its inertia is(6,1,2). Including the seven passive roots gives(13,1,2).

There are84 six-element supports among these nine active normals. Exact
congruence checks find15 positive definite supports, with the inverse
of each positive block checked entrywise nonnegative. Each support plus
the seven passive normals is an actual finite vertex, by the argument
below. The certificate saves every positive and nonpositive support.

The raw25 prism fixes qA,w and the passive roots, and has exactly eight
vertices. The raw23 prism fixes qB and the passive roots, and has exactly
ten. Their union contains14 of the15 vertices. The remaining vertex is
the support `{w,h5,hcenter,h3,B,A}` plus passive roots, of type
`3H3+H4`. The ninth positive support in the raw25's entire wall set is
not a ninth prism corner: it omits one of that prism's fixed supporting
walls. It is already a vertex of the raw23 prism. The source explicitly
separates fixed prism supports from arbitrary positive principal subsets.

## Why the completions are actual faces

Let h be a strict future timelike interior point of the actual acute
chamber, and let n_i be any selected13 distinct walls with positive
definite Gram G. Set b_i=h·n_i<0 and c=G^{-1}b. A positive definite
matrix with nonpositive off-diagonal entries has nonnegative inverse,
so c≤0. For

\[
x=h-\sum_i c_i n_i
\]

one has x·n_i=0, and every other actual wall m satisfies

\[
x\cdot m=h\cdot m-\sum_i c_i(n_i\cdot m)\le h\cdot m<0.
\]

Also x²=h²−bᵀG^{-1}b<0 and x·h=x²<0. Thus x is future timelike,
and the selected13 walls determine its unique projective finite vertex.
Every additional ambient facet is strictly inactive there. This proves
the actual completion, rather than assuming that a positive principal
Gram is automatically realized in an abstract diagram.

Fix the11 walls `{w,h5,hcenter,h3}` plus the passive roots. The four
positive supports adjoining the pairs

`{qB,qA}, {qB,A}, {B,A}, {qA,B}`

give a compact quadrilateral. Its minimal five-root core is
`{qB,qA,w,B,A}`; it has inertia(3,1,1), and its unique relation has
alternating signs on the four side roots. The same cone argument as in
the [frozen full-star aggregation proof](fresh_quad_product_aggregation_20261017.md)
therefore gives exactly these four extreme rays. They are actual finite
vertices, strict against all other walls. Convexity makes every other
ambient inequality redundant on this compact quadrilateral.

At its new corner the marked native5 pair is A,B. The H3 roles are
`w--3--A--5--B`, so A is the center and B the5-leaf. The opposite corner
is the original source vertex.

One additional triangle is essential for counting. Fix the11 walls
`{qB,qA,h5,hcenter}` plus the passive roots. The three remaining side
walls are w,h3,C. All three pair completions are among the saved15
positive supports. Their residual3-by-3 Gram has inertia(2,1,0), and its
inverse is strictly negative entrywise. Its dual rays are the three
finite vertices on a common future cone, proving that these halfspaces
bound the compact triangle. At the source vertex its two edge ports
are w,h3 and its unique nonincident side wall is C. The source saves
the common support, all three corners, and the source-corner identity.

## An injection into three pointed target corners

Let X count simultaneous occurrences of the specified raw25 and raw23
corner orbits at source vertices of type1325. The current complete
compatibility graph has exactly12 allowed literal pairs between these
two orbits. The full finite forest has12 weighted automorphisms, and the
source checks that their images of stars12/46 are exactly those12 pairs,
with native contacts and occupied port pairs preserved. Thus the preceding
geometry applies to every possible selected coexistence. It does not
multiply the whole-face count by12.

Map such an occurrence to its forced quadrilateral, pointed at the new
native5 corner. A marked H3 component there fixes A,B,w in their distinct
roles. The marked pair A,B determines the unique actual two-face through
that vertex, hence the entire quadrilateral and its opposite source
vertex. It also fixes the w edge port at that source vertex.

The literal triangle proved above contains this w port and has C as its
unique nonincident wall. By
[disjoint maximal simplex ports](fresh_simplex_disjoint_ports_20261017.md),
there is only one maximal simplex through that port. The triangle lies
in it and shares its unique nonincident facet C. Consequently C is
recovered uniquely from the pointed quadrilateral and source vertex.
This step uses an actual triangle; a square port alone would not suffice.
The unordered external wall sets{A,C} and{B,C}, together with the source
vertex's13 incident walls, now recover both complete quadrilateral-star
aggregates uniquely. Each dotted star has its unique maximal prism by
the frozen aggregation theorem. Thus two distinct originating corner
pairs cannot map to the same pointed target.

Every vertex of type6537 has three H3 components, each with exactly one
eligible native5 pair. Its H4 component does not have the required whole
H3 role. It follows that

\[
X\le 3v_{6537}.
\tag{2}
\]

This geometric incidence already underlies row50 of the
[52-context simplex family](fresh_sqrt5_core23_25_spectator_incidence_20261017.md),
whose inequality is `2a26448+a26918−2a311≤3v6537`.
We explicitly recognize that older row; the new information is its
coupling to the complete source-vertex corner graph below.

## The source-vertex charge

The eight current corner orbits at1325 have whole-face indices

`45769,26918,26918,26362,26362,26448,45573,26683`.

Each is one distinct corner occurrence of its Q class. Give these
orbits respective integer weights

`3,1,2,4,4,4,0,0`.

For an actual set of simultaneous stars at one source vertex, let x be
the number of selected raw25/raw23 pairs. The source re-enumerates all522
maximal cliques of the unchanged conservative compatibility graph and
checks, for every one,

\[
\text{corner charge}\le4+x,\qquad x\le1.
\tag{3}
\]

All unknown pairs remain allowed. Every actual set is a clique. Removing
a star from a maximal clique cannot invalidate(3): the12 marked pairs
form a matching, so at most one x is lost, and its endpoint weight is
at least1. All other weights are nonnegative. This proves(3) for every
actual subset, without supposing that an attaining relaxation mixture
describes a physically realized chamber.

Summing(3) over source vertices counts the Q faces with coefficients
`3,3,8,4` as in(1). Combining with(2) proves that inequality. Unlike the
older simplex row, it uses the other corner classes competing for the
same source-vertex allowance. The old row is strictly satisfied by both
current witnesses, while(1) is violated.

## Exact repair and scope

Both face-mass vectors give left side `353962329/5000000000` in(1).
The full witness violates it by `1558119/5000000000`; the named248-Q
restriction violates it by `9223961/27500000000`. Increasing only μ1325
by one quarter of the respective gap makes(1) an equality. The exact
repaired ratios are

\[
\frac{55714268997452065}{4643208969323058}<12,
\qquad
\frac{18670750401222348}{1548855999321229}>12.
\]

All other vertex masses and all face masses are unchanged. Increasing
one nonnegative vertex allowance preserves every inherited column,
including the52 columns with signed face coefficients: their variables
are nonnegative, their vertex costs are nonnegative, and their fixed
face loads do not change. The same monotonicity preserves the historical
weighted capacities and graph-cover witnesses. No old-mass implicit
packing or joint-family coverage is claimed at these newer face masses.

No pair is excluded or has its ownership changed. In particular the two
owned unknowns at1763 remain. There is no global optimizer or inference
of target12 feasibility. This package completes one finite geometry and
one pointed incidence inequality. Global compatibility, additional actual
classes and complete face/field coverage remain separate.

Generation99156 and fresh saved replay8426 both terminated exit0/PASS,
in9.32 and7.27 seconds respectively. All35 direct source bindings, eight
local links and Python syntax checks passed.

Construction_audit completed a separate full source/report mathematical
read. It passed the actual strict-center completion, the source-incident
triangle and its compact cone, the maximal-simplex recovery of C,
native5 role uniqueness, the3-per-target count, literal fixed-support
distinction for the8/10 prism corners, the12-map and522-clique interfaces,
matching-based downward closure, and both monotone repairs. Root also
completed a full source/report mathematical read, including the inherited
simplex-port and full-star aggregation proofs, and passed the same
geometric and counting composition. These are source/math reviews using
inherited interfaces, not independent replays of the84 supports,12 maps
or522 cliques.

A separate [standard-library saved-data checker](../scripts/check_forest1325_cross_vertex_saved_20261017.py)
imports no producer or project arithmetic. It independently enumerates
all522 maximal cliques from the saved graph by a different recursion,
checks all eight count patterns and12 matching pairs, recomputes every
charge and the four whole-face coefficients, and verifies both exact
repairs using both vertex-budget coefficients. Its first final-location
run terminated exit0/PASS (tool chunk0eced9), and a final report-link
check rerun also terminated exit0/PASS (chunk70cfd2). It also checked the35
source bindings and local links. This checks the saved graph and exact
accounting; it does not independently reconstruct the native walls,
84 supports or12 orientation maps, or prove the pointed injection.
Source, data and report are **FINAL/FROZEN**. No external refereeing or
proof-assistant verification is claimed.
