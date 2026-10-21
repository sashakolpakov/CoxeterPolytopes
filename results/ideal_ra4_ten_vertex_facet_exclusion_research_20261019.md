# The ten-vertex octahedrite cannot enlarge the octahedral chain class

19 October2026. New direct cusp-geometry argument, awaiting independent
review. No checker, classification search, numerical solver, Lean job or
TeX build ran. No dimension exclusion is claimed. The current verified
16-page ideal4 manuscript is unchanged.

## Statement

Let O be the regular ideal right-angled octahedron, D its triangular-face
double, A_k the ideal right-angled k-antiprism, and B the ten-vertex
octahedrite in Kolpakov's Figure7. Suppose every facet of an entirely
ideal finite-volume right-angled hyperbolic four-polytope P is one of

    O, D, B, or A_k for any k≥4.

Then neither B nor any A_k occurs. Hence P satisfies the previously proved
O/D chain theorem. In particular the complete published small-facet
classification implies:

**If every facet of P has at most ten vertices, then P is a chain of
regular ideal24-cells with the original opposite-cut restriction.**

For an entirely ideal right-angled five-polytope, some four-facet must
therefore have a three-face with at least eleven vertices. More precisely
each positive-kappa four-facet must contain a three-face outside the
displayed spectrum, and every five-face in a higher-dimensional entirely
ideal polytope inherits this condition.

These are restricted structural conclusions and necessary conditions,
not nonexistence in dimensions5 or6.

## 1. Inputs and the intrinsic B model

The only three-dimensional classification input is the face lattice of B
(and, for the at-most-ten corollary, completeness of the small list).
Standard uniqueness of ideal right-angled three-polyhedra realizes its
combinatorial automorphisms as isometries. The cusp shapes below are
derived explicitly, rather than inferred from a picture or numerical fit.
The inherited O/D and antiprism cusp descriptions are as in
[the all-antiprism argument](review_and_extension_ideal_ra4_antiprism_spectrum_20261018.md).

The vertices of B are A,...,J. Its four quadrilaterals and eight triangles
are

    ABCD, ABFE, EHIJ, CGHI;
    AEJ, ADJ, BFG, BCG, CDI, DIJ, EFH, FGH.

The permutations

    S1=(A I)(B H)(C E)(D J)(F G),
    S2=(A B)(C E)(D F)(G J)(H I)

preserve these faces. Their product is the half-turn fixing C,E and
interchanging the opposite pairs A/H,B/I,D/G,F/J. Its action on the
rectangle at E exchanges both opposite face pairs, so it is a genuine
half-turn rather than a reflection. Send E to infinity and C to0 on the
boundary plane of the upper-half-space model of H3. That half-turn is
(x,y)↦(−x,−y). Scale and rotate the rectangle at E so its four vertical
walls are

    ABFE: x=a,    EHIJ: x=−a,
    AEJ:  y=1,    EFH:  y=−1,                       (1)

with a>0. Set

    R=(1+a²)/2,    t=(1−a²)/2.

The hemisphere for ABCD is orthogonal to x=a and contains C=0 and
A=(a,1), so its center is(a,t), radius R. Its half-turn image CGHI
has center(−a,−t), radius R.

The BCG hemisphere passes through0 and is orthogonal to both these
hemispheres. If its center is(q,v), this gives aq+tv=0. It is tangent
to x=a at B, since those two nonadjacent facet walls are opposite in
the rectangular link at B. Thus B=(a,v) and q=(a²−v²)/(2a).
Solving gives v=1 or−a²; v=1 would identify B with A. Consequently

    B=(a,−a²),    center(BCG)=(at,−a²), radius aR.

The four vertices of ABFE appear on its boundary line in the order
A,B,F,E=∞. Hence −1<−a²<1, so0<a<1. In particular t>0.
The half-turn gives I=(−a,a²) and the CDI hemisphere with center
(−at,a²), radius aR.

The ADJ hemisphere is orthogonal to y=1 and tangent to x=±a at A,J,
so its center is(0,1), radius a. The DIJ hemisphere is orthogonal to
x=−a and passes through I=(−a,a²), J=(−a,1). Its center is(−a,R),
radius t. These two hemispheres are orthogonal, as required along DJ.

At D the ABCD and DIJ hemispheres are opposite in the cusp rectangle,
and hence externally tangent. Their centers and radii give

    4a²+(R−t)²=(R+t)²,
    a^4+4a²=1.

Write phi=(1+sqrt5)/2 and s=sqrt(phi). The positive solution is

    a²=sqrt5−2=phi^−3,   a=s^−3,
    R=phi^−1,           t=phi^−2.                   (2)

The remaining hemispheres are half-turn images: FGH has center(0,−1),
radius a, and BFG has center(a,−R), radius t. The remaining boundary
vertices are

    F=(a,−1), H=(−a,−1), J=(−a,1),
    D=(−a³,1−2a²), G=−D.

This derivation takes place inside the existing B, so it assumes no
unverified claim that a candidate hemisphere arrangement is complete.

## 2. All B cusp signatures, including which sides are long

There are three intrinsic classes of cusps:

|Cusps|Rectangle aspect|Quadrilateral sides|
|---|---|---|
|A,B,H,I|s|Two adjacent sides, one long and one short|
|D,F,G,J|s|One long side|
|C,E|s³|Both opposite long sides|

In particular **every B cusp is nonsquare**. The two displayed
automorphisms are transitive on each row, so it suffices to compute A,J,E.

At E the rectangle is2a by2, proving aspect1/a=s³. Its x=±a sides
belong to the quadrilaterals and are long.

At A translate A to0 and invert in a unit Euclidean sphere. The two
vertical walls become x=0,y=0. The ABCD and ADJ hemispheres become
planes at distances1/(2R) in the y-direction and1/(2a) in the
x-direction. Thus the rectangle has width1/(2a), height1/(2R),
and aspect R/a=s. The quadrilateral ABFE is a short side, ABCD a
long side, and the triangle AEJ is a long side.

At J the analogous rectangle has width1/(2a), height1/(2t),
so aspect a/t=s. Its unique quadrilateral EHIJ is long; AEJ is
short. At E that same triangle AEJ is short. We have established
the particular cross-cusp signature needed later:

    AEJ: long at A, short at E and J.               (3)

Every B cusp requires at least one long quadrilateral side. The
opposite-quadrilateral cusp requires two. No cusp of B is square.

## 3. Complete multiplicative list of possible box shapes

In a cusp box of P, each nonsquare face ratio belongs to

    Rset={s,s³,2} union {rho_k=2cos(pi/k): k≥4}.

The minimum is s, the maximum s³; rho4=sqrt2, rho5=phi=s²,
rho6=sqrt3, and rho_k increases to2. For three distinct side lengths
x>y>z, the two consecutive ratios r=x/y and q=y/z are in Rset,
as is their product. The only unordered possible pairs are

    {s,s}, {s,rho5}, {rho4,rho4}.                   (4)

Here is a complete check, without a cutoff on k. A factor s³ times
any unequal ratio already exceeds s³. The same holds for2, since
2s>s³. Two antiprism ratios have product at least2; unless both
are rho4, their product is at least sqrt2*phi>s³ because phi<2.
For s*rho_k, all k≥6 give s*sqrt3>s³ since3>phi². The k=5
product is s³. The k=4 product s*sqrt2 lies below s³ but is
neither s,s³ nor2. It cannot equal an antiprism ratio either:
if s*sqrt2=2cos(pi/j), then cos²(pi/j)=phi/2 generates Q(sqrt5),
so phi(j)=4 and j is5,8,10 or12. The latter two quadratic fields
Q(sqrt2),Q(sqrt3) exclude8,12; the explicit cosine-squares at5
and10 are(3+sqrt5)/8 and(5+sqrt5)/8, neither phi/2. This excludes
the last possibility. The pair s*s gives rho5, completing(4).

Up to scaling the three-distinct shapes are therefore exactly

    (2,sqrt2,1),
    (s²,s,1), (s³,s²,1), (s³,s,1).                 (5)

The first shape uses D and A4 faces and is excluded by the previous
long-side argument: the D rectangles force both long sides of an A4
rectangle to be nontriangular, while an A4 cusp has only one.

Each of the three new shapes has exactly one opposite pair of A5
rectangles (ratio s²) and four B rectangles (ratios s or s³).
But an A5 cusp is incident to its pentagonal base. At least one
edge of its rectangle must therefore correspond to a pentagonal
ridge of P, shared with an adjacent facet. Every adjacent facet
in this box is B, which has no pentagonal faces. This is impossible.
Thus **no three-distinct box occurs**.

## 4. Equal-long boxes fail, and a B triangle finishes the proof

A noncubic box now has proportions r:r:1 or r:1:1 for r in Rset.
The first shape has four nonsquare belt rectangles and two square
end rectangles. An O cusp has no nontriangular side and a square
D cusp has exactly one, so the two ends together permit at most
two nontriangular edge incidences. Every belt rectangle requires
at least one long nontriangular side, for B by Section2 and for
D/A_k by the inherited signatures. In an r:r:1 box all those long
sides lie on the square ends. They demand at least four such edge
incidences, a contradiction.

Consequently every noncubic box has shape r:1:1. If r=s or s³,
all four nonsquare belt faces belong to B and both square ends
belong to O or D. Within a B rectangle a long side meets another
belt face, so a triangular ridge on a long side has B as its
neighbor. A short side meets an end face, so its neighbor is O or D.

Apply this to one actual AEJ triangle of a B facet. At A its long
side forces the neighboring facet to be B. At E and J its short
sides force that same neighboring facet to be O or D. It is a
single ridge and has one neighboring facet globally, so this is
impossible. No B facet occurs.

We are now in the already proved O/D/all-antiprism class. Its
all-antiprism exclusion removes every A_k, proving the theorem.
This order avoids extending a statement which assumed B absent.

## 5. Higher-dimensional consequence and verification boundary

For an ideal right-angled four-polytope Q put

    kappa(Q)=3v(Q)−4f(Q)+24.

The inherited ideal5 incidence identity is sum_Q kappa(Q)=48, so
some four-facet has positive kappa. An m-cell O/D chain has
kappa=2−2m≤0. Therefore each positive-kappa four-facet has a
three-face outside O,D,B and all antiprisms. By the complete
published classification through ten vertices, that face has
at least eleven vertices. Every ideal5-face has this property;
in particular every five-facet of an ideal6-polytope does.

The external inputs here are the stated small3D classification,
three-dimensional right-angled uniqueness, and the previously
reviewed O/D chain and all-antiprism theorems. The new work is the
exact Figure7 cusp geometry, full box-ratio classification and
same-triangle neighbor contradiction. None is claimed kernel-
verified or proved by a finite model census. The dimension ceiling
for entirely ideal right-angled polytopes remains6.
