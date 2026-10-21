# Selective distance12 and the unrestricted dimension136 bound

Every finite-sided, finite-volume real hyperbolic Coxeter polytope has
dimension at most136. The argument has no arithmeticity, compactness,
simple-cusp, factor-rank or facet-count restriction. For a chamber having
a compact facet, the ambient dimension is at most103.

This strengthens the [dimension146 proof](fresh_selective_distance13_dimension146_20261017.md)
by replacing distance13 with distance12. As before, pairs in the same
connected corner component of order at most17 are also paid. The new
exception consists of two minimal cores sharing a terminal edge. A
Lorentzian Schur-complement argument isolates that core, replacing the
radical argument used for the previously known shared-vertex exception.

The [producer](../scripts/fresh_selective_distance12_dimension136_20261017.py)
and [certificate](fresh_selective_distance12_dimension136_20261017.json)
give the exact finite case split and the remaining scalar inequalities.
The [separate review](fresh_selective_distance12_dimension136_20261017_review.md)
independently reconstructs the critical arithmetic and reviews the written
geometry. No new finite-volume construction or publication-priority claim
is made.

## Geometric setup and inverse tests

Use the inherited actual-face and cusp-factor incidence statements from
the146 proof. For a quadrilateral, let W be its supporting white roots
and U,V the opposite pairs of side roots. W and W plus one black root
are positive definite. W plus an adjacent black pair is positive
semidefinite. Each of W+U and W+V is indefinite, since the corresponding
closed opposite sides do not meet. Minimal indefinite cores L,M have
order at most10 and both black deletions positive definite.

The complete high-rank catalogue has384 marked contexts and maximum
distance7 from the marked pair. Cores of order at most5 have radius at
most3, independently of their edge labels. Assume there is no path of
length at most12 between U and V with white internal nodes.

For a black-deletion Gram A, every external supporting white contacting
node i must have Schur residual at most2−(A inverse)ii. The normalized
contact magnitude is at least1 and A inverse is entrywise nonnegative.
If either black deletion has inverse diagonal at least2, this contradicts
positive definiteness of W plus the remaining black. Thus equality2 is
excluded here too. The producer recomputes all marked distances and all
inverse tests at distance at least5. Exactly18 marked ports survive.

## Disjoint cores

The cores cannot be orthogonal in a space of negative index one. Every
cross bond must join white ports whose distances from their black pairs
sum to at least12; both distances are therefore at least5. The18-port
list leaves115 ordered marked joins, with six underlying possibilities:

| Ports, ignoring marked-context multiplicity | Ordered occurrences |
|---|---:|
|L25 node9 to L25 node9|60|
|L25 node9 to L25 node8|24|
|L25 node9 to L22 node8|18|
|L25 node9 to L18 node7|6|
|L25 node9 to L22 node4|6|
|L25 node8 to L25 node8|1|

L25 has simple edges0–1,0–2,2–3 and0–4–5–6–7–8–9.
The node9-to-node8 join contains the same wholly white affine E8 block
as in the146 proof. The node8-to-node8 join contains wholly white affine
D5: its two trivalent nodes are the joined nodes8,18 and their four
leaves are7,9,17,19. Positive kernels certify both obstructions. Larger
bond magnitudes or additional nonpositive contacts only decrease the
quadratic form on those kernels, so neither actual white block can be
positive definite.

After excluding those types, each marked context pair has only one
possible cross bond. Each endpoint has a wholly white simple A3 tail.
Their Schur endpoint residuals are4/3, so positive definiteness forces
the squared cross-bond magnitude below16/9. The only nonright Coxeter
bond in this range is simple. This reasoning covers arbitrary finite
labels as well as parallel or dotted candidates; there is no cutoff.

Delete the two fork centers, numbered0 and10, from each resulting
joined Gram. The retained block is positive definite. Its two-center
Schur matrices are as follows:

| Right core/port, with left L25/node9 | Schur matrix | Determinant |
|---|---|---:|
|L18/node7|[[-1/15,-1/10],[-1/10,-7/30]]|1/180|
|L22/node8|[[-5/66,-1/11],[-1/11,-7/44]]|1/264|
|L22/node4|[[-1/15,-1/10],[-1/10,-1/5]]|1/300|
|L25/node9|[[-7/78,-1/13],[-1/13,-7/78]]|1/468|

Every matrix is negative definite. Thus every disjoint case has two
negative directions and is impossible.

## Complete overlap split

Every shared root must have its two distances sum to at least13.
Consequently both distances are at least6. Exactly16 marked contexts
have such roots, and each has at most two of them. Enumerating injective
partial identifications of these roots therefore exhausts the overlap;
no intersection of three or more nodes is possible.

There are208 candidate identifications before injectivity and common
Gram-entry consistency, and162 after those tests. In120, a shared node
has an external neighbor from the other core and its black-deletion
inverse diagonal is at least2. That neighbor is white, since the shared
root's distance from the other black pair is at least6. The inverse
test excludes these identifications.

The42 remaining identifications all use two L25 cores:

* 27 share only node9 with node9. Their marked pairs belong to
  {0,1},{0,2},{0,3},{1,2},{1,3},{2,3}, with at least one pair from the
  last three (radius7).
* Six share just node9 with node8, or node8 with node9. The shared node
  has wholly white arms of lengths1,2,5, so the positive affine-E8
  kernel excludes these cases, including any extra contacts.
* Nine share the terminal edge in reverse: left8=right9 and left9=right8.
  Both marked pairs are from {1,2},{1,3},{2,3}.

For the shared-vertex cases, the only additional cross bond not already
creating a path of length at most12 would connect the two neighbors of
the shared node. These three whites would be an affine triangle; its
positive kernel excludes the bond. In the shared-edge cases every
additional cross bond immediately gives a short marked path. Therefore
the two surviving core Grams are completely fixed.

## Isolation of the shared-vertex core

This is the19-root Gram of the146 proof. Number the first copy0..9 and
map the second's0..8 to10..18, with9 shared. Its Schur matrix after
deleting centers0,10 is the constant2×2 matrix with entries−1/12. Its
radical is the actual root relation

    z=(-6,-3,-4,-2,-5,-4,-3,-2,-1,0,
        6,3,4,2,5,4,3,2,1).

The relation is actual because its Gram-null linear combination is
orthogonal to a timelike vector in the core span, and hence must vanish.

Choose side0 for a marked pair containing0. For the remaining pairs
choose2 from {1,2},1 from {1,3}, and2 from {2,3}; make the analogous
choice in the second copy. In all27 contexts the two preferred sides
lie in one corner component of order at most17. Their distance in that
component is13 or14.

The producer decomposes each full preferred corner into connected
finite or affine blocks. An additional supporting white e has no
contact with an affine block, by its positive kernel. In a finite block,
positive semidefiniteness permits contact only at inverse diagonals at
most2. Together with the two omitted blacks, every allowed contact has
nonzero z-coordinate. Every negative-z contact is within distance2 of U;
every positive-z contact is within distance2 of V.

Pairing e with z gives zero. Acuteness therefore forces contacts of both
signs whenever there is any contact. This gives a white-interior marked
path of length at most6, contrary to the assumption. Every additional
supporting white is orthogonal to the core.

## Isolation of the shared-edge core

The new core has18 roots. Number the first copy0..9; map the second's
0..7 to10..17, its8 to9, and its9 to8. Deleting centers0,10 leaves a
positive-definite block. The two-center Schur matrix is

    [[-5/66,-1/11],[-1/11,-5/66]],

with negative determinant. Thus the full core is nonsingular Lorentzian;
its determinant is−1. Use the same preferred side choices for the nine
radius7 marked contexts. Their16-node preferred corners are one A16,
four D16 and four affine D15, and the preferred marked distance is13.

The preferred-corner inverse or positive-kernel tests allow additional
supporting-white contacts only in the end sets

    E_left={1,2,3}, E_right={11,12,13}.

Each available end node is within distance1 of its marked pair. If e
contacts both ends, it gives a marked path of length at most4. Otherwise
all contacts lie at one end. In the inverse of the FULL Lorentzian Gram,
the principal block on either end set is

    [[8,10,5],
     [10,14,7],
     [5,7,4]].

All entries are positive and the diagonal minimum is4. For a nonzero
contact vector b supported at that end, acuteness and the minimum
nonzero magnitude1 imply b^T G^{-1} b>=4. The core span already has the
unique negative direction, so its orthogonal complement in the ambient
space is positive definite. The additional root's orthogonal residual
must therefore satisfy

    0 <= 2-b^T G^{-1}b <= -2,

a contradiction. Every additional supporting white is again orthogonal
to the whole core. This step applies to arbitrarily large allowed
contact magnitudes and needs no radical relation.

At finite preferred corners the isolated marked component is the actual
ambient vertex component. At ideal preferred corners, the inherited
factor-fiber incidence theorem puts the entire varying affine factor
inside W plus the two sides. An isolated finite component cannot be
such a factor; an isolated affine component is already the full factor.
Additional omitted spectator cusp walls do not enlarge it. This covers
nonsimple cusps and completes the selective quadrilateral payment.

## Dimension inequalities

Pay weight1 on a marked pair if its graph distance is at most12 OR its
connected component has order at most17; otherwise pay0. Every triangle
has two paths of length at most8 by the inherited minimal-core argument.
The preceding proof pays every quadrilateral at least1. Nonnegative
weights suffice for polygons with five or more sides. Thus every
k-gonal two-face receives total weight at least5−k.

Components of order at most17 have at most8 times their order in pairs.
For larger finite or affine components, the all-rank path/cycle/fork
counts give at most the component order at each distance. The finite
vertex cost is consequently at most12n. At an ideal vertex of factor
ranks r_i and link vertex count V=product(r_i+1), the factor-fiber
multiplicity gives

    C(v)/V <= sum_i min(12,r_i/2).

The inherited shape inequality gives ideal normalized curvature

    kappa/V = 3d/4 - (5/2) E|sum_i X_i|,   d=n-1,

where each X_i is centered uniform on r_i+1 consecutive integers.
The already proved all-rank convex-order transfer applies with cap24:
concentrating two ranks below24, or two above24, preserves this cost
and cannot increase curvature. Increasing the sum of squared ranks
proves termination. The reduced partitions are copies of24, at most
one part1..23, and at most one part at least25.

For d=136,...,145, the1175 reduced profiles all have positive
curvature-minus-cost. Integer convolution and the independent moment
formula agree. The minimum is172096/1953125, at ranks(24,24,24,24,24,16).
The finite slack beta_n−12n is positive throughout n=137,...,146; the
global shape inequality also has its strict positive5(n−1) term.
Hence these dimensions are impossible. The completed146 theorem
excludes all larger dimensions, proving the bound136.

At n=136 the profile(24,24,24,24,24,15) instead has margin−59157/3906250.
This same scalar comparison does not prove135. It is a local link
profile, not a finite-volume chamber construction.

For a compact facet of dimension d=n−1 the finite payment gives
beta_d<12(d+1). For d=103+2t and104+2t, eight times the opposite
difference has coefficient lists(12,208,4) and(112,212,4), respectively.
Both are positive for t>=0. Thus d<=102 and ambient n<=103.

The finite/quasi-Lanner classification, actual incidence results,
all-rank distance counts, convex-order theorem and shape inequality
remain explicit inherited inputs. The finite arithmetic does not
formally verify these written arguments or replace their hypotheses.

Primary generation and a fresh saved replay both terminate with PASS.
The separate Fraction-only implementation also passes all384 marked
contexts,75 ports,162 consistent overlaps,115 disjoint joins,36 corner
domains and1175 profiles. Its written mathematical/source review finds
no gap. These are internal checks, not external refereeing or formal
proof-assistant verification. The primary source, data and report are
final/frozen; the separate review records its own replay and final
bindings.
