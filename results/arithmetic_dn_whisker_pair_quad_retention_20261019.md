# A complete candidate whisker-pair unit refund in D19/D20

19 October 2026. **New written proof, awaiting independent adversarial
review.** No enumeration, exact-checker run, numerical solve, Lean
process or build was run. All earlier source-bound reports are unchanged.

This note completes the quadrilateral consumers of the third natural
D_n unit, the pair of short whiskers. Together with the separate
triangle-retention lemma it gives a candidate actual one-unit refund.
The numerical dimension ceiling is not improved by this statement.

## 1. Precise local and weight statements

Let n be 19 or 20 and let P be a compact arithmetic simplest-type
hyperbolic Coxeter polytope. At a whole D_n finite vertex, write a
for the fork, b,e for its short whiskers and a-l1-...-lL for its
long arm. Put N=n-2 and L=N-1. Unit outward roots are used.

Let B be the distance-at-most-two finite-angle baseline, with the
usual nonnegative H4 extra if applicable. At every vertex whose
ENTIRE rank-n diagram is D_n, consider its three distinguished
marked pairs:

    (a,b), (a,e), (b,e).

Each has baseline 1. Define B* by deleting that unit on all three
pairs at every such vertex and leaving every other baseline intact.
This is only a comparison function: it is not claimed to pay every
face, since the golden/silver fork-pair exceptions are known.

**Quad-retention claim.** Any actual quadrilateral containing a
modified WHISKER-pair angle (b,e) has total B* at least 1.

The [separate triangle lemma](arithmetic_dn_whisker_pair_triangle_retention_20261019.md)
proves that any actual triangle containing such an angle has total
B* at least 2. Therefore the whisker-pair deletion has all its
triangle and quadrilateral consumers paid even if both fork units
are also deleted. It composes, as far as these consumers are
concerned, with any separately justified subset or activation of
the fork deletions.

In particular, from any correct angle rule at least B, subtracting
one only at the full-D_n whisker pairs preserves all face payments.
All remaining faces are unchanged, and k-gons with k>=5 require
only nonnegativity. The retained units proved here are baseline
units, not adjustable bonuses. This conclusion requires the proof
below and its explicitly listed inherited local lemmas to pass
independent review before promotion.

## 2. Geometric preliminaries

The white support for a quadrilateral at (b,e) is the connected path

    W={a,l1,...,lL}=A_N.

Write its side roots in cyclic order b,e,c,d; the opposite pairs
are (b,c),(e,d). Put

    b.c=-x,     e.d=-y.

Every distinct actual wall pair is acute. A product -1 is impossible
in a compact chamber: the positive null sum of the two unit normals
is feasible against every other wall by acuteness, giving an ideal
point. Thus a nonzero native magnitude is a finite cosine or is
strictly greater than 1.

At a full finite vertex the old roots span a positive n-space.
Every new wall has a nonzero timelike component, because the simple
old vertex already has n incident facets and lies strictly inside
the new wall. Thus old-plus-new determinant ratios are negative and
nonzero physically. At every nonphysical embedding of the totally
real arithmetic ground field the corresponding independent block
is positive definite. We use cyclic integrality and those strict
conjugate determinant inequalities. Unit entries are not assumed
to lie in the ground field without a cycle argument.

Any positive-definite root subset of size n determines an actual
finite vertex of P. To see this directly, project an interior future
timelike point perpendicular to those normals. Their positive-definite
acute Gram matrix has an entrywise nonnegative inverse; the correction
is a nonnegative combination of outward normals. It preserves all
other strict wall inequalities and keeps the point timelike.

Consequently four cyclic candidate corners that are positive definite
and share a white support really close an actual quadrilateral.
For each chosen side the two finite candidate vertices are the two
endpoints of its edge in the white face; every other wall is strict
there. A convex polygon cannot have an additional vertex between
them on the same side. The four edges form its entire boundary.

## 3. The two complete adjacent-corner lists

Assume for contradiction that all other corners of this quad have
B*=0. The corner W+e+c before c is added is the long path

    e--a--l1-- ... --lL=A_(n-1).

If its ORIGINAL baseline is zero, finite classification and the
marked distance condition give exactly these c contacts with W:

|Profile|Contact with W|
|---|---|
|I|None|
|F|Simple bond to lL|
|B|A 4-bond to lL|
|T|Simple bond to l_(L-1)|

Indeed finite one-root extensions of a path of this rank are
isolated, endpoint 3/4 extensions, or next-to-end simple extensions.
The near endpoint e is excluded by adjacency positivity; the near
next-to-end node a is excluded by its marked distance 2 to e.
Two contacts create a cycle and exceptional finite diagrams are
too small. Also e.c=0.

The only other way B*=0 at (e,c) is that it is one of the deleted
distinguished pairs at a whole D_n vertex. It cannot be a fork/
whisker pair: e already has its white neighbor a, and adding c gives
it at most degree two, while a fork needs three neighbors; if e
were the leaf of fork c, its two neighbors a,c would violate that
role. It can only be a WHISKER pair. That forces

    c.a=-1/2,     c perpendicular to every other white root,
    c.e=0.                                      (H)

Its common fork with e is a. The same classification holds for d
at (d,b), exchanging b,e and c,d.

Thus each adjacent new zero is either original-zero I/F/B/T, or
the new near-whisker profile H. These lists are complete, and none
of the adjacent new zeros is a fork/whisker deletion.

## 4. Three inherited algebraic wall tests, with their exact scopes

We will use the following local facts.

**Far-wall finite-product test.** If a new wall has profile F,B or T
against the old white path and contacts one old whisker with a finite
native cosine, while being perpendicular to the other whisker, it
is impossible at n=19,20. This is the algebraic content of
[the nonisolated finite-product note](arithmetic_dn_nonisolated_quad_finite_opposite_obstruction_20261019.md),
independent of its c contacts. A zero whisker contact gives an
old-plus-new positive-definite block (F) or a positive-null block
(B,T), and is impossible too. Parallel products -1 are forbidden by
Section 2. Therefore every such actual far wall has magnitude >1
against that old whisker.

**Isolated-whisker finite labels.** If a new wall contacts only one
old whisker by -cos(pi/m), the only possible finite m at n=19,20
are 5 and 12. Here is a short complete proof. Put u=cos^2(pi/m).
The old-plus-new ratio gives 0<u'<2/n<=2/19 at every nonphysical
embedding. For d=[Q(u):Q],

    Tr(u)=d/2+mu(m)/4>=d/2-1/4.

If d>=3 this exceeds 1+(d-1)*2/19, an upper bound from the
conjugates; the difference is (30*d-87)/76>0. Thus d<=2 and
m=3,4,5,6,8,10,12. Rational u forces K=Q, excluded by rational
isotropy for compact simplest-type chambers in these dimensions.
The other conjugate for m=8 or 10 exceeds 2/19. Only m=5,12 remain.
This fact concerns squared contacts in K, not an unsupported claim
that cos(pi/m) itself is in K.

**Near-whisker finite labels.** For a wall of profile H, let its
opposite old-whisker contact have magnitude x. Its old projection
cost is

    q_c=n*x^2/2+N*x+N/2.

The simple old/new triangle cycle puts 2*x in O_K. At a nonphysical
embedding q_c'<1 gives

    -1<x'<-1+4/n.

If x is finite, A=2*(1+x) is integral, lies in [3,4) physically
and in (0,8/n) otherwise. Degree at least 3 gives norm less than
4*(8/19)^2<1. The resulting finite cosine list of degree at most 2
is m=3,4,5,6. Label 3 makes an induced affine triangle at a and the
chosen old whisker. Labels 4,5 violate the displayed conjugate
interval. Only m=6 is possible, with x=sqrt(3)/2. Its actual ground
field must be Q(sqrt(3)): any nonidentity extension embedding fixing
x would preserve the physically negative old-plus-new determinant.

These tests use no published high-dimensional ground-field list.
They use the arithmetic criterion, ordinary finite classification,
elementary cosine traces and rational isotropy.

## 5. Neither adjacent corner has profile H

Suppose c,d each have original-zero profile I/F/B/T. The finite
corner W+c+d allows at most one nonisolated profile, except F+F:

* F+F produces a D_n fork at lL, with c,d its short whiskers.
* Every other two-nonisolated combination has a 4-bond in a fork,
  two forks, degree four, or arms of lengths 1,2,and a long arm;
  none is finite in these ranks. Equivalently, add the second root
  to the large A/B/D component made by the first.

The F+F case makes the opposite angle (c,d) a deleted full-D_n
whisker pair. Each of x,y must be >1 by the far-wall test. But the
principal Gram on {b,c,e,d} is the direct sum of their two blocks

    [1,-x;-x,1],     [1,-y;-y,1].

Both would have a negative eigenvalue. A Lorentz space has only one
negative direction. Thus F+F is impossible in an actual polytope.

If both profiles are I, c contacts only b and d only e. Projecting
off W gives b_bar^2=e_bar^2=t=n/[2*(n-1)] and
-b_bar.e_bar=h=(n-2)/[2*(n-1)], with t+h=1. The unchanged
c,d are orthogonal unit roots. Opposite projected walls are strictly
ultraparallel, so x^2>t and y^2>t. The rank condition is

    (x^2-t)*(y^2-t)=h^2.

Both x,y finite would make both positive factors less than h;
both >1 would give two negative black blocks. Neither can equal 1.
Thus one is finite and the other >1. This is precisely the
[reviewed double-whisker obstruction](arithmetic_dn_double_whisker_obstruction_20261018.md),
which excludes the pair in these dimensions without a field list.

It remains one I and one far profile. By symmetry let c be I and
d far. The far-wall test gives y>1. The two-negative-block argument
then forces x<1, so the b.c contact is finite and nonzero.
These same actual walls form a second quadrilateral, now at the
old fork/whisker angle (a,e), with white support

    W'={b,l1,...,lL}=A1+A_L.

Its four n-root corner diagrams are

    W'+a+e = old D_n,
    W'+e+c = I2(m) on b,c + A1(e) + A_L,
    W'+c+d = I2(m) on b,c + (A_L+d),
    W'+d+a = the finite far extension of the path b-a-l1-...-lL.

All are positive definite. The geometric observation in Section 2
makes this a genuine second actual quadrilateral. Its other marked
pairs (e,c),(c,d),(d,a) all have baseline zero: the first two lie
in different components and the last pair is far apart in its
large classical component. The new far-wall contact is F,B or T.
The corresponding original-unpaid fork quadrilateral is excluded
by the [far-simple](arithmetic_dn_far_simple_unpaid_quad_obstruction_20261019.md),
[far-double](arithmetic_dn_far_double_quad_obstruction_20261019.md),
or [far-fork](arithmetic_dn_far_fork_quad_obstruction_20261019.md)
note. This is an acyclic import: those nonisolated original-unpaid
lemmas do not use a whisker refund, an activation rule or the present
statement. Their proofs in this branch use no ground-field list.

This exhausts the case with no H adjacent corner.

## 6. Exactly one adjacent corner has profile H

By symmetry take c to have H and d to have I/F/B/T. Thus

    c.a=-1/2, c.e=0, c has no other white contact,
    d.b=0.

The opposite angle cannot be a deleted fork/whisker pair: c already
has its white neighbor a and a new marked neighbor could give it
only degree two. If it were a deleted whisker pair, d would also
attach a simply and the other adjacent angle would have H, contrary
to the present case. Therefore it must be originally baseline zero,
in particular c.d=0.

The two opposite native products are x=-b.c and y=-e.d. Neither
is zero, and neither is 1. The first zero would create the affine
D4 star at a; the second is excluded by opposite-wall geometry for
I, or by Section 4 for far profiles. At most one of x,y exceeds 1
by the two-negative-black-block argument.

### 6.1 The profile d is I

The old D_n projection costs and cross cost are

    q_c-1=(x+1)*[(N+2)*x+N-2]/2,
    q_d=(N+2)*y^2/2,
    q_cd=N*y*(1+x)/2.

The two-root Schur determinant is zero. Cancelling its nonzero
factor x+1 gives the exact relation

    z:=4*y^2=2*[(N+2)*x+N-2]/[(N+1)*x-1].       (2)

If x is finite, Section 4 gives x=sqrt(3)/2 and K=Q(sqrt(3)).
The length-two cyclic product makes z integral in K. Direct
quadratic-field conjugation of (2) gives

    Tr(z)=4+4*(7*N-1)/(3*(N+1)^2-4),

strictly between 4 and 5 for N=17,18. This is not the trace of an
algebraic integer, a contradiction.

If x>1, the black index bound forces y<1. The isolated-whisker
finite-label test gives m=5 or 12 for y. Put u=y^2. Inverting (2)
gives

    x=[2*u+N-2]/[2*(N+1)*u-(N+2)].              (3)

The triangle cycle at a,b,c makes 2*x integral. For the two allowed
values of u, (3) has the exact traces

|u|Tr(2*x)|Strict interval at N=17,18|
|---|---|---|
|(3+sqrt(5))/8|4+20*N/(N^2-5)|(5,6)|
|(2+sqrt(3))/4|4*(7*N-1)/(3*(N+1)^2-4)|(0,1)|

Neither is an integer. No square root of u was presumed to belong
to the ground field. This excludes d=I.

### 6.2 The profile d is F,B or T

The far-wall test gives y>1. Hence x cannot exceed 1 and Section 4
forces x=sqrt(3)/2, K=Q(sqrt(3)). Define

    A=(N+1)*x-1,     J=(N+2)*x+N-2.

For d profiles F,B,T respectively put

    (p,h)=(1/2,1/2), (sqrt(2)/2,1), (1,1).

Here q_d=(N+2)*y^2/2+2*p*y+h and the old cross projection is

    q_cd=(1+x)*(N*y/2+p).

The Schur determinant gives, after cancellation,

    2*A*y^2+4*p*(x-1)*y+(h-1)*J-2*p^2*(x+1)=0.

For F, B_F=2*y is in O_K by the simple long cycle, and this is

    A*B_F^2+2*(x-1)*B_F-[(N+3)*x+N-1]=0.        (4)

For T put B_T=2*y; for B put B_B=2*sqrt(2)*y. The corresponding
simple/doubled long cycle puts the indicated B in O_K, and both
equations become

    A*B^2+4*(x-1)*B-4*(x+1)=0.                 (5)

In particular each indicated B belongs to K, although sqrt(2) need
not belong to K in the B case. The discriminants, with x=sqrt(3)/2,
are

    Delta_F=3*N^2+8*N+20+2*(N^2-N-6)*sqrt(3),
    Delta_BT=12*(N+2)+8*(N-2)*sqrt(3).

Their field norms are nonsquares:

|N|Norm(Delta_F)|Consecutive square comparison|Norm(Delta_BT)|Comparison|
|---:|---:|---|---:|---|
|17|197457|444^2<197457<445^2|8784|93^2<8784<94^2|
|18|210496|458^2<210496<459^2|8448|91^2<8448<92^2|

For example Delta_F is 1023+532*sqrt(3) or 1136+600*sqrt(3);
Delta_BT is 228+120*sqrt(3) or 240+128*sqrt(3). If (4) or (5)
had a root in K, its discriminant would be a square in K, whose
norm would be a rational square. The displayed positive integers
are not rational squares. This excludes all three far profiles and
completes the exactly-one-H case.

## 7. Both adjacent corners have profile H

Now c and d both attach a simply and have no other white contacts;
c.e=d.b=0. Finiteness of W+c+d forces c.d=0, since otherwise
{a,c,d} would be a cycle or an indefinite pair. The opposite corner
is therefore also a whole D_n whisker pair: all four marked angles
are deleted whisker units in this potential case.

Each opposite native product x,y is finite or >1, and the finite
case can only be cos(pi/6) by Section 4. Neither product is zero
(affine D4 at a) or 1 (compactness). Both >1 would give two negative
black blocks.

Project off W. All four black roots have square
t=n/[2*(n-1)] and their white cross projection is
h=(n-2)/[2*(n-1)]. Their native four-root matrix consists of the
two blocks [1,-x;-x,1], [1,-y;-y,1], in pairs b,c and e,d.
Its white correction is h times the all-ones matrix. The determinant
lemma, or a direct four-by-four expansion, gives the rank condition

    (N+1)*x*y-x-y-(N-1)=0.                     (6)

There is no division by a zero factor here, since x,y!=1. If both
are finite, x=y=sqrt(3)/2 and (6) is
(7-N)/4-sqrt(3)=0, impossible for N=17,18.

If exactly one is finite, let x=sqrt(3)/2. Equation (6) gives

    y=[x+N-1]/[(N+1)*x-1].

The simple triangle cycle at a,e,d makes 2*y an algebraic integer
in K=Q(sqrt(3)). But

    Tr(2*y)=4*(7*N-1)/(3*(N+1)^2-4) in (0,1),

a contradiction. This exhausts the last case.

## 8. Composition and exact boundary of the result

All possibilities for newly zero adjacent corners were I/F/B/T/H,
and the zero opposite corner was either originally zero or a
distinguished full-D_n pair. The degree tests exclude opposite
fork/whisker pairs; the remaining original-zero, far-whisker and
near-whisker cases have all been considered. Therefore no actual
quad containing the old whisker pair can have all four B* values
zero. Since B* is nonnegative and integer-valued, its total is at
least 1, proving the quad-retention claim.

Together with the separate triangle lemma this yields the following
candidate refund theorem, pending review: subtracting the single
short-whisker-pair baseline unit at every whole D19 or D20 vertex
preserves every small-face payment of any correct angle rule that
dominates the baseline. Its consumers retain the required baseline
units even under the maximum possible fork-unit deletions. Thus it
can compose with a separately proved safe fork activation. It does
not make the unsafe unactivated fork deletions valid.

No published high-dimensional field list was used. The inherited
nonisolated original-unpaid fork lemmas, far finite-product test and
reviewed isolated double-whisker theorem use the arithmetic cyclic
criterion, positive nonphysical forms, finite classification and
rational isotropy; their relevant branches are the ones identified
above. All new traces and discriminant norms are displayed exactly.

This refunds only one unit at a full D20 vertex, taking its bare
charge from 38 to 37. Even a separate one-unit guaranteed fork
activation would leave 36. A uniform target of 35, all other rank-20
forests, native contact coverage and eventual exact/formal verification
remain separate work. The compact arithmetic ceiling is unchanged.
