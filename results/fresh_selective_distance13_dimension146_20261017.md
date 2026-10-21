# Selective distance13 and the unrestricted dimension146 bound

Every finite-sided, finite-volume real hyperbolic Coxeter polytope has
dimension **at most146**. There is no arithmetic, compactness, simplicity,
cusp-rank, factor-count or facet-count restriction. This improves the
[internally checked156 proof](fresh_distance14_unrestricted_dimension156_20261017.md).
A compact-facet application also improves the ambient bound127 to111.

The argument uses a selective corner payment: weight1 if the two marked
nodes have graph distance at most13, **or if they lie in one connected
component of order at most17**; weight0 otherwise. Distance13 alone is
false. Its exceptional configuration is the known finite-volume H17
pyramid. The crucial new step is proving that the relevant small corner
component cannot acquire additional supporting white attachments without
creating a short paid path. No new construction dimension is claimed.

The [producer](../scripts/fresh_selective_distance13_dimension146_20261017.py)
and [exact certificate](fresh_selective_distance13_dimension146_20261017.json)
check the geometric matrices and1255 reduced cusp profiles. The
[separate review](fresh_selective_distance13_dimension146_20261017_review.md)
records the independent checking and its scope. External refereeing and
proof-assistant formalization are not asserted.

## Actual quadrilateral stars

Use the actual-face incidence and closed-opposite-side arguments from the
[distance theorem](fresh_simple_cusp_distance_bound_20261017.md).
For a quadrilateral, W is the positive-definite set of supporting whites,
and U,V are its opposite black pairs. W plus one black is positive
definite; W plus adjacent blacks is positive semidefinite. Each W+U and
W+V is indefinite. Choose inclusion-minimal indefinite cores L,M in them.
Both black deletions of each core are positive definite. Minimal cores
have order at most10, and the complete high-rank census supplies384
marked contexts, all of radius at most7. Low-order cores have radius at
most3, independently of their labels.

Suppose no path of length at most13 joins U to V with only white internal
nodes. If the cores overlap, every shared node must have distance7 from
both marked sets: any smaller sum gives the forbidden short path. The
only radius7 contexts are L23 with blacks{0,1}, far node8, and L25 with
blacks{1,2},{1,3},{2,3}, far node9. Each has only one farthest node.

For a positive-definite black deletion A, its inverse is entrywise
nonnegative. A new white contacting its node i nontrivially has contact
magnitude at least1 and squared norm2. Its Schur residual is at most
2−(A inverse)ii. If that inverse entry is at least2, the enlarged block
cannot be positive definite, contradicting its containment in W plus
one black. This excludes equality as well as strict negativity.

The producer recomputes every marked distance and all21 distance-at-least6
ports in16 contexts. The only ports surviving both black-deletion tests
are L25 node9 for the six pairs

    {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3},

and L25 node8 for{1,2}. In particular L23's far node8 cannot contact an
external white. If L23 shared that node with M, a neighbor of the shared
node in M would be an external white: another shared node would already
give a shorter marked path. The same inverse test excludes this overlap.

## Disjoint cores cannot be exceptional

Disjoint cores cannot be orthogonal, since the ambient Gram has negative
index at most one. Any cross bond must have endpoint distances summing
to at least13. Both endpoints are white. The preceding inverse tests
leave only two underlying possibilities, with both cores of type L25.

For a far9-to-far9 bond, the white A3 tail on each side has residual4/3.
Positive definiteness of their union bounds the squared bond magnitude
by16/9, below2. Thus the nonright bond is simple, including when arbitrary
finite labels, parallel or dotted pairs were initially allowed. Removing
the two fork centers from the resulting20-root Gram leaves the positive
block A12+2A2+2A1. Its two-center Schur complement is

    [[−7/78, −1/13], [−1/13, −7/78]],

with two negative eigenvalues. This is the obstruction already used in
the156 proof. If there were additional far-to-near cross bonds, the next
obstruction applies instead.

For a far9-to-near8 bond, both marked contexts have radius7. There is a
wholly white T-shaped subgraph with arm lengths1,2,5, namely affine E8.
For the bond9--18 in two copies numbered0..9 and10..19, its nodes are

    [18,19,17,16,9,8,7,6,5].

The certificate supplies its positive kernel. Every native nonright bond
has magnitude at least1. Increasing that bond or adding other nonpositive
contacts makes the quadratic form on the positive kernel no larger.
Thus the actual white block cannot be positive definite. Inducedness of
the displayed graph is unnecessary. This excludes every far-to-near
case, including simultaneous cross bonds.

Consequently the only no-short-path case is two L25 cores sharing node9.
They share no other node and have no additional cross bond, since either
would supply a path of length at most13.

## The shared core and its isolation

Number its first copy0..9 and its second copy10..18, using9 for the shared
far node. Each L25 has edges

    0--1, 0--2, 2--3, 0--4--5--6--7--8--9,

with all bonds simple; the second copy uses the analogous labeling.
The19-root Gram has signature(17,1,1) and radical vector

    z=(−6,−3,−4,−2,−5,−4,−3,−2,−1,0,
        6, 3, 4, 2, 5, 4, 3, 2, 1).

This is an actual linear relation among the roots, not just a relation
in an abstract Gram matrix. Their span contains a timelike vector. The
Gram-radical combination is null and orthogonal to that vector, hence
lies in a positive-definite space and is zero.

For the left marked pair choose a preferred side as follows:

| Marked pair | Preferred side | End type |
|---|---:|---|
|{1,2}|2|Unbranched|
|{1,3}|1|Fork|
|{2,3}|2|Fork|

Make the analogous choice on the right. Deleting the two unchosen black
nodes gives a connected17-node corner: A17 for two unbranched ends,
D17 for one fork, or affine D16 for two forks. There are respectively
one, four and four marked contexts. Its preferred sides are distance14.

Let e be any further SUPPORTING white of the actual quadrilateral. In
an A17 or D17 preferred corner, positive semidefiniteness after adding e
forces e to be orthogonal to every node whose inverse diagonal exceeds2.
For A17 the possible contacts are only the two nodes nearest either end;
for D17 only the long-tail endpoint and its neighbor. The inverse is
entrywise nonnegative, so multiple contacts cannot evade these tests.
For affine D16, its positive kernel forces e to be orthogonal to all17
corner roots.

Including the two omitted black roots, all remaining possible contacts
are at distance at most1 from U on the left or V on the right. Every
left candidate has z-coordinate strictly negative and every right
candidate strictly positive; the central zero-coordinate node is
excluded. Pairing e with the actual radical relation gives

    sum_i z_i (e,e_i)=0.

All products are nonpositive. A nonzero contact therefore requires at
least one contact on each side. These give a path from U to V through e
of length at most4, with only white internal nodes. This contradicts
the assumption. Every additional supporting white is consequently
orthogonal to the whole core.

At a finite preferred corner, all incident roots are W plus its two
sides, so its marked connected component has exactly17 nodes. At an
ideal corner, use the actual [cusp-factor incidence lemma](fresh_cusp_factor_payment_framework_20261017.md):
the varying affine factor containing both side marks lies WHOLLY in
W plus those two sides. Thus an isolated finite A17 or D17 component
cannot be the varying component of an ideal corner. The affine D16
case is already an entire irreducible cusp factor. Missing spectator
walls at a nonsimple cusp cannot enlarge the paid marked component.

This proves the selective payment for every quadrilateral. Triangles
already have two paths of length at most8 by the inherited minimal-core
argument, so receive at least2. Larger polygons require only nonnegative
payment. Every k-gonal two-face therefore receives at least5−k.

## Vertex costs and the dimension contradiction

A component of order m<=17 has at most m(m−1)/2<=8m pairs in total.
For a larger component the payment counts only distances at most13.
The inherited all-rank path, cycle and single-/double-fork formulas give
at most m pairs at each distance; exceptional finite/affine diagrams
have order at most9. Thus finite vertex cost is at most13n.

At an ideal vertex let the irreducible factor ranks be r_i and let
V=product(r_i+1). The same pair bound and the trivial total-pair bound,
multiplied by the exact factor-fiber multiplicity V/(r_i+1), give

    C(v)/V <= sum_i min(13,r_i/2).

The [shape-curvature theorem](fresh_cusp_link_shape_curvature_20261017.md)
gives normalized ideal curvature

    kappa/V = 3d/4 − (5/2) E|sum_i X_i|,

where d=n−1 and X_i is centered uniform on r_i+1 consecutive integers.
Use the all-rank convex-order transfer already proved in the156 package,
now with cap26. Concentrating two ranks below26 until one is0 or26, or
two above26 until one is26, preserves the cost and cannot increase the
curvature. The strictly increasing sum of squared ranks proves
termination. Every partition is therefore dominated by copies of26,
at most one part1..25, and at most one part>=27.

For d=146,...,155 there are121,...,130 such profiles,1255 in total.
Integer convolution and a separate bounded-composition moment formula
agree exactly. All curvature-minus-cost margins are positive; the
minimum is830407/9565938 at ranks(26,26,26,26,26,16).
The finite slack beta_n−13n is positive throughout n=147,...,156.
The global shape lower bound also has the positive5(n−1) term. These
inequalities contradict the corner-payment upper bound. The inherited
156 theorem excludes n>=157, completing the unrestricted bound146.

At n=146, the valid local link profile(26,26,26,26,26,15) instead has
margin−515185/19131876. Thus this scalar comparison alone gives no145
bound; that profile is not an ambient finite-volume example.

For a compact facet of dimension d=n−1, all vertices are finite and
its strict shape inequality gives beta_d<13(d+1). For d=111+2t and
112+2t, eight times the difference has coefficients respectively
(12,224,4) and(120,228,4) in increasing powers of t. Both are positive
for t>=0. Hence d<=110 and ambient n<=111 in this restricted class.

## Known geometric sharpness and verification scope

The shared core is a finite-volume H17 pyramid with19 facets, already
within [Tumarkin's classification, Theorem2 and Table13](https://arxiv.org/pdf/math/0301133).
The producer independently reconstructs its entire chamber cone in
Gram-dual coordinates p_i=(x,e_i). Its image is z-perpendicular and its
inequalities are p_i<=0. The extreme rays are the single zero-coordinate
ray and the81 two-coordinate rays pairing opposite signs of z. Exact
inverse-Gram calculations give80 finite vertices and2 ideal vertices,
all in one time cone. Their sum is a strict timelike interior point;
the tight rays of each wall span a facet. Thus the complete projective
polytope lies in the closed Klein ball and has finite volume. This is
an exact reconstruction of a known example, not a novelty claim.

It realizes the distance14 obstruction. The unrestricted improvement
comes from its forced small component, not from asserting distance13
for every quadrilateral. The producer verifies384 marked contexts,
21 inverse-budget ports, nine preferred-corner contact domains, the
physical Gram obstructions,82 chamber rays and1255 profile margins.
Published classifications, inherited actual incidence arguments and
the written all-rank reductions remain explicit mathematical inputs;
finite matrix checks do not replace those arguments.

Producer generation and a fresh saved replay both terminate with PASS.
Source syntax and all19 direct source bindings pass. A separate internal
mathematical/source review and an independent Fraction-only implementation
also pass, including the nine preferred-corner domains and all1255
profiles. The separate review binds this final producer package and
records its own saved replay. No mathematical gap was identified by
these internal reviews; external refereeing remains unperformed.
