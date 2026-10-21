# Perpendicular-edge facet capacity in compact right-angled H4

18 October 2026. Own proof attempt before further literature searching.
Written geometric and combinatorial lemmas; independent review pending.
No computation, optimizer, compilation, Lean execution, minimum-volume
theorem, or second-volume theorem is claimed. Frozen files are unchanged.

## Statement and consequence

Let P be any compact right-angled hyperbolic four-polytope, with f facets.
If one facet F has V vertices, then

    f >= 2V+5.                                    (1)

Write W for the total ridge excess above pentagons. Since W includes all
the ridge excess belonging to F, the elementary three-dimensional Euler
identity gives

    W >= d(F)=(V-20)/2.

Consequently

    f+W >= (5/2)V-5.                              (2)

The exact compact volume identity is Vol(P)=pi^2(f+W+16)/12. A hypothetical
polytope of volume strictly below34pi^2/3 therefore has V<=48 in every
facet. Here V is even. If Vol(P)<=34pi^2/3, all its facets instead have
V<=50. Likewise a hypothetical polytope of volume strictly below the
120-cell double has V<=104 in every facet, since that volume corresponds
to f+W=256. These are finite local-size reductions; they do not classify
the surviving polytopes or settle either volume question.

## Compact acute-normal facts used below

Use future unit vertex vectors and outward unit wall normals. Distinct
normals have product0 when their facets meet and otherwise product<-1.
The weak gap is the standard right-angled acute-normal gap. Equality-1
would give a future-null feasible sum of the two outward normals by
acuteness, contradicting compactness, so the strict version holds here.

Pairwise orthogonal wall normals have a common actual face: starting from
an interior future timelike point, orthogonally remove its components in
those normals. Its norm stays negative; adding nonnegative multiples of
outward normals keeps all other wall inequalities strict by acuteness.
The projected point is future and belongs to the indicated face.

Fix a facet A with normal a. Any two wall normals b,c orthogonal to a
and to the normal r of a facet disjoint from A must be orthogonal to one
another. Indeed span(a,r) is Lorentzian, so its orthogonal complement is
positive definite. The acute right-angled gap excludes b.c<=-1 there:
product-1 would make b=-c, incompatible with an interior point satisfying
both inequalities strictly; a product below-1 is impossible in a positive
space. Thus b.c=0.

In particular, a facet disjoint from A meets a clique of the neighbors of
A. These facts prove the flag and induced-square exclusions needed here;
no arbitrary combinatorial realization is assumed.

## Perpendicular-edge opposites are distinct

This argument applies in both dimensions three and four. Let A be a
compact right-angled facet and x one of its vertices. The one incident
edge leaving A is perpendicular to A. Its other endpoint is finite.
The facet R_x meeting that endpoint and not containing the edge has normal

    r_x=sinh(L_x)x-cosh(L_x)a,   L_x>0.             (3)

Here x is regarded in a-perp and L_x is the length of that edge. Formula
(3) follows because r_x is orthogonal to all the wall normals defining x
inside A, has norm1, and faces the correct half-space. It shows R_x is
disjoint from A. Its intersection with any neighbor B of A is nonempty
exactly when x belongs to B: its product with B's normal is sinh(L_x)
times the corresponding product with x.

The R_x are pairwise distinct. If r_x=r_y, their projections to a-perp
give sinh(L_x)x=sinh(L_y)y, forcing x=y. Moreover an actual facet
disjoint from A meeting all its neighbors at x must equal R_x. Its normal
has the same form sinh(L)x-cosh(L)a; two distinct such normals have
positive product cosh(L-L_x), contrary to acuteness. This also avoids
assuming a more distant supporting wall creates another endpoint on the
same edge.

## A polygon in a compact right-angled three-polyhedron

Let Q be a compact right-angled three-polyhedron containing a k-gonal
face A. Label its neighboring faces A_1,...,A_k cyclically, and label by
B_i the perpendicular-edge opposite at the vertex A_i cap A_(i+1) of A.
Indices are cyclic. The previous paragraph makes all these faces distinct.
The acute-normal facts show that the only neighbors of B_i among the
A_j are A_i,A_(i+1), and that the A_j induce precisely their cyclic k-gon.

Call any further face extra. An extra face meets at most one A_i. Two
nonconsecutive A_i cannot both meet it by the clique property. Two
consecutive A_i would force it to be the already counted B_i.

**Three-dimensional capacity lemma.** Q has at least2k+2 faces. If it has
exactly2k+2, its unique extra face meets none of the A_i.

Here is the required closing argument. Work in the triangulated dual
two-sphere of Q. Every vertex in it has degree at least5, because every
face of Q is a compact right-angled polygon.

If there are no extra faces, each A_i has exactly the five known neighbors
A,A_(i-1),A_(i+1),B_(i-1),B_i. Its cyclic link therefore forces the edge
B_(i-1)B_i. The outer triangles form a disc whose complementary disc has
boundary B_1,...,B_k and no interior vertices. In a triangulation of a
polygon, an ear has degree2 within the disc. That B_i has only its two
known outside neighbors A_i,A_(i+1), so its total degree is4, a contradiction.

If there is exactly one extra face C and it meets A_i, it meets no other
A_j. For j!=i the same five-neighbor argument forces B_(j-1)B_j. At A_i,
the cyclic link instead contains B_(i-1),C,B_i in this order. The remaining
disc is a polygon with boundary the B-cycle with C inserted between
B_(i-1) and B_i; it has no interior vertices. Every B_j needs degree at
least3 within this disc, since it has exactly two outside A-neighbors.
C needs degree at least4 there, since it has only one outside neighbor
A_i. An ear again contradicts these requirements. Thus C meets no A_j.

For completeness, in this last case C is the sole interior vertex of the
remaining B-disc. A chord between two B_j would cut off a triangulated
polygon not containing C and produce an ear of total degree4. Thus there
are no such chords, and C joins every B_j. The face lattice is the usual
Loebell combinatorics with two k-gons and2k pentagons. Only the already
proved disjointness of C from all A_j is needed in (1).

## Four-dimensional counting with singleton footprints

Return to F in P. It has M=V/2+2 two-faces, and hence M neighboring
facets F_1,...,F_M. Together with F and the V distinct R_x from (3), these
give1+M+V known facets.

Let the ridge F cap F_i have k_i sides. Inside the three-polyhedron F_i
we have already counted exactly the following faces:

* F cap F_i itself;
* its k_i neighboring faces coming from the other F_j;
* the k_i perpendicular-edge opposites coming from the R_x at its vertices.

The perpendicular edges used in this three-polyhedron are the same actual
edges of P used in (3). The orthogonal-projection fact ensures that a
different F_j meeting F_i must already be one of those k_i neighbors.
Thus every further face of F_i comes from a facet of P outside the known
1+M+V facets. Let e_i be their number. The three-dimensional capacity
lemma proves e_i>=1.

Every such extra facet has a footprint of at most two among the F_i.
Its footprint is a clique by the acute-normal facts. Three F_i in that
clique meet F at a vertex x, since their normals in a-perp are orthogonal.
An actual facet meeting all three is R_x by (3), already counted. A
larger footprint is impossible because a-perp has signature(3,1), so it
contains at most three mutually orthogonal unit spacelike normals.

If e_i=1, the equality part of the capacity lemma says its unique extra
face in F_i is disjoint from every one of the k_i neighboring faces of
F cap F_i. Its parent extra facet cannot meet a different F_j: such a
meeting would give three pairwise orthogonal walls and hence the actual
common edge inside F_i with F_j. Thus this extra facet has a singleton
footprint {i}.

Let s count the i with e_i=1. Their extra facets are distinct and meet no
other F_j. All other M-s neighbors require at least two incidences with
extras, each of footprint at most two. There are therefore at least

    s + (1/2)sum_(e_i>=2)e_i >= s+(M-s)=M

distinct extra facets. Adding the known facets gives

    f >= 1+M+V+M = 2V+5,

as asserted. The absence of ideal vertices is essential to the finite
perpendicular-edge endpoints and the strict normal gap in this argument.

## What remains unproved

There is no propagation here from one dodecahedral facet to an entire
120-cell, and no proof that all compact right-angled four-polytopes are
unions of regular120-cells. The capacity lower bound starts at45 facets
for V=20; it does not reach120. Equations (1)--(2) reduce possible
smaller-volume counterexamples to facets with at most48 vertices, but do
not remove those possibilities. A proof of the minimum or the double gap
must control compatibility among those remaining facets, beyond this
single-facet neighborhood count.
