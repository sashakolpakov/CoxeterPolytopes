# Triangle consumers of the whole D19/D20 fork units

19 October 2026. **New written consumer theorem, pending independent
review.** No computation, census, proof script or build was run. This note
does not classify quadrilaterals, construct a compact chamber, give a valid
full refund, or improve a dimension bound.

## 1. Statement and exact scope

Let P be a compact arithmetic simplest-type hyperbolic Coxeter polytope of
dimension n=19 or20. At a finite vertex assume its ENTIRE diagram is D_n.
Write its fork as a, its short leaves as b,e, and its long arm as
a-l1-...-lL, where L=n-3. All old roots are unit roots and all old bonds
are simple. Consider the actual two-face at the fork/whisker angle (a,b).

Start from the distance-at-most-two baseline. In the Q(sqrt(3)) branch use
the actual inherited high-label rule, which also gives a unit to an entire
I2(m) marked component for m>=6. Suppose the unit is deleted at EVERY
whole-D_n fork/whisker angle, including other corners of this same face.

**Triangle consumer theorem.** If this two-face is a triangle, its total
remaining charge is at least two, except for the following exact profile:

    the third wall c is perpendicular to every old root except b,
    b.c=-cos(pi/5).

That exceptional profile forces actual ground field Q(sqrt(5)). Retaining
the original unit at that angle pays it. The new theorem concerns the
whole D_n units only; it does not alter any other basic or bonus role.

The required background is finite Coxeter classification, the actual
finite-corner correspondence and acuteness, absence of affine principal
diagrams in a compact chamber, arithmetic cyclic integrality and
nonphysical positive definiteness, and rational isotropy in dimension at
least five. No high-dimensional field LIST or attachment census is used.

## 2. Complete finite continuation at one triangle corner

Put W={e,l1,...,lL}. A triangular two-face has third side c, and its three
whole vertex diagrams are W+a+b, W+a+c and W+b+c. In particular W+a+c
is finite. Its old path is

    e-a-l1-...-lL = A_(n-1).

A finite diagram is a forest. Thus c contacts at most one root of this
path. At ranks19,20, finite irreducible classification permits exactly:

1. no contact;
2. a simple or4 contact to an endpoint e or lL;
3. a simple contact to a next-to-end node a or l_(L-1).

Indeed the connected possibilities of rank n are A_n,B_n,D_n; the E,F,H
diagrams have smaller rank. There are no suppressed contacts or spectator
components: this is the entire old finite vertex and the entire triangle
support W. The remaining b.c contact is finite or right, since W+b+c is
finite.

## 3. Contact at the near endpoint e

Suppose c.e has label3 or4. If b.c is nonright, finite classification of
the three-node path e-c-b permits ordered labels (3,3),(3,4),(3,5),(4,3).
At the two other triangle corners, (a,c) is a distance-two pair and (b,c)
is an adjacent pair. Neither is a deleted whole-D_n fork/whisker unit.
They alone pay two, even after deleting the old (a,b) unit.

If b.c=0, label3 at e produces an induced affine E8: its arms about a
have lengths1 (b),2 (e,c),5 (l1,...,l5). Label4 produces the positive-null
five-root diagram on {a,b,l1,e,c}; positive kernel coefficients in this
order are

    (1, 1/2, 1/2, 1, 1/sqrt(2)).

Both contradict compactness. These are principal obstructions and do not
depend on unlisted roots.

## 4. Contact at the far endpoint lL

Let its label be k=3 or4. If b.c=0, label3 makes the complete old-D_n-plus-c
diagram finite D_(n+1), violating the positive index n. Label4 gives a
positive-null diagram: put coefficient1 on a and every old long-arm root,
coefficient1/2 on b,e, and coefficient1/sqrt(2) on c. This is incompatible
with compactness.

If b.c is nonright, W+b+c contains the long path b-c-lL-...-l1 of rank
n-1. It can only be A_(n-1) or B_(n-1), so k=3 and b.c has label3 or4.
Label3 gives the induced all-simple affine cycle

    a-b-c-lL-...-l1-a.

For label4 put beta=sqrt(2)/2. The independent old-D_n-plus-c block has
projection cost

    q=n beta^2/2+1/2+beta.

The cycle above now has one4 bond; its cyclic product puts sqrt(2) in the
ground field. At an embedding changing its sign, the projection cost is

    q'=n/4+1/2-sqrt(2)/2 > 1.

This contradicts arithmetic nonphysical positive definiteness for both
n=19,20. Therefore the far-endpoint branch is empty.

## 5. Contact at the far next-to-end node l_(L-1)

The contact is simple. If b.c=0, the complete diagram is affine D_n: put
coefficient1 on a,l1,...,l_(L-1) and coefficient1/2 on b,e,lL,c.
If b.c is nonright, the component in W+b+c has a fork at l_(L-1) with
arm lengths1,2,L-2. Its rank is n-1>=18; it is neither a path nor a D
diagram, and cannot be a finite E,F,H diagram. Hence W+b+c is not finite.
This eliminates the whole branch.

## 6. Contact at the near next-to-end node a

Here a.c=-1/2 and c is perpendicular to W. If b.c=0, the principal
diagram on {a,b,e,c,l1} is affine D4. Otherwise write b.c=-beta with
beta=cos(pi/m), m>=3. The two corners W+a+b and W+a+c are both whole
D_n fork/whisker corners, so BOTH units are being deleted in this audit.

The old-D_n projection cost of c is

    q=[n beta^2+2(n-2)beta+n-2]/2,
    q-1=(beta+1)(n beta+n-4)/2.

The black three-cycle makes 2beta an algebraic integer in the ground field.
Set A=2(beta+1). Its physical value lies in [3,4). At every nonphysical
embedding, finite two-root positivity gives |beta'|<1, and q'<1 then gives

    0<A'<8/n <= 8/19.

Any embedding of the full ground field fixing beta would preserve the
negative physical determinant ratio, so this ground field is Q(beta).
For degree at least three, the positive nonzero integer norm of A would
be smaller than 4(8/19)^2<1. For degree two its norm is1 and trace4,
so A=2+sqrt(3). Thus beta=sqrt(3)/2 and m=6. In degree one A=3, giving
m=3, whose black triangle is affine and violates compactness.

The only remaining case has ground field Q(sqrt(3)). Its third corner
is W+I2(6), and the actual high-label rule gives its marked I2(6) pair
TWO units. Both D_n fork units may therefore be deleted on this triangle.
The local arithmetic feasibility of m6 is not a deficient payment.

## 7. No contact with the long A_(n-1) path

Here c is perpendicular to a and W. If b.c=0, the old-D_n-plus-c block is
positive definite of rank n+1, impossible. Otherwise put beta=cos(pi/m),
m>=3, and z=4beta^2. The physical projection cost is q=nz/8>1. Cyclic
integrality makes z an algebraic integer, with

    1<=z<4,           0<z'<8/n<=8/19

at all nonphysical embeddings. The same independent-block argument makes
the ground field exactly Q(z).

Degree at least three has impossible norm less than4(8/19)^2. In degree
two the positive norm is1 and the integer trace is2,3 or4. Trace2 gives
no degree-two value. The other two possibilities are

    z=(3+sqrt(5))/2,  m=5,  K=Q(sqrt(5));
    z=2+sqrt(3),      m=12, K=Q(sqrt(3)).

In degree one, z=1,2 or3 and K=Q; rational isotropy excludes a compact
arithmetic simplest-type chamber in these dimensions. This uses that
global arithmetic input, rather than incorrectly rejecting the local
rational Lorentz root block by conjugate tests alone.

The m12 triangle has a whole I2(12) marked component at its third corner.
The inherited Q(sqrt(3)) high-label rule pays two there. The m5 triangle
has payments (1,0,1) before deletion, and is precisely the exception in
Section1. Its local model and absent saved second-node A1 cross bonus are
spelled out in the [independent triangle/quad review]
(review_dn_triangle_quad_capacity_finite57_20261019.md) and the [shared
capacity theorem](arithmetic_dn_golden_triangle_quad_capacity_20261019.md).

## 8. What is proved and what remains

Sections2–7 exhaust actual finite continuations of a third triangle wall
at this whole-D_n angle. With the actual field-specific charge, the sole
newly deficient triangle is the golden singleton profile. The shared
capacity theorem shows that it and the golden opposed-fork quad cannot
occupy both whisker angles at one old vertex. The silver quad forces a
different ground field and cannot mix with the golden triangle.

This settles the TRIANGLE part of a prospective activated refund. It
does not prove the corresponding complete quadrilateral interface.
Every affected quad with inadequate remaining charge, all simultaneous
deletions, and any new compensating coupons must still be covered before
there is a valid refund theorem. Even a full one-unit saving at every
whole D20 vertex leaves cost37, above the target35 for excluding dimension20.
