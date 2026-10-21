# Mixed A/D end consumers in the eight high-label fields

19 October 2026. **New written composition lemma candidate, awaiting
independent review.** Paper and source work only: no scientific job,
enumeration, optimizer, build or Lean process was run. Every earlier
source-bound report is unchanged.

The purpose is narrow. The proposed whole-A_n two-end distance-two
refund and whole-D_n three-unit refund have the SAME disconnected
white graph. Their separate consumer proofs therefore do not compose
automatically. This note classifies and excludes the genuinely mixed
newly unpaid quadrilaterals. It does not assert a global allowance35
or a new dimension ceiling. Both standalone refund sources must also
pass independent review before the composition is promoted.

## 1. Fixed dimension, fields, rule and standalone dependencies

Let n be19 or20 and let P be a compact arithmetic simplest-type
Coxeter polytope of dimension n. Its ground field belongs to

    Q(sqrt(3)), Q(sqrt(6)), Q(sqrt(2),sqrt(3)),
    K9,K11,K15,K16,K20,       Km=Q(cos(2pi/m)).

All have degree at most5. None is Q, Q(sqrt(2)), Q(sqrt(5)),
Q(sqrt(13)), Q(sqrt(17)) or Q(sqrt(21)). Start with the proved
eight-field distance-at-most-two charge, including the extra unit
on an ENTIRE I2(m) component with m>=6. Consider these deletions,
only when the vertex's ENTIRE diagram has the indicated type:

* A_n: the first/third marked pair at each end of the path;
* D_n: the two fork/short-leaf pairs and the short-leaf pair.

The sources for the separate consumer assertions are

* `arithmetic_an_two_end_distance_two_refund_20261019.md`,
  SHA256 `0404c5501536fb9cc09453a214c3ed15b5028a4fe0d839fce6568be0fce817b0`;
* `arithmetic_dn_three_unit_refund_eight_fields_20261019.md`,
  SHA256 `4147a57aac91f048eb81668892725f8f1c87cf86898cb55841a0f51b9378faf6`.

The harmless B-row upper-bound wording correction in the former's
same-endpoint dependency is recorded separately in
`arithmetic_an_same_endpoint_b_row_clarification_20261019.md`.
This note does not alter those sources or substitute for their
independent reviews.

Every triangle still has its unchanged entire high-label I2 corner
of weight2; that corner is neither a whole A_n nor a whole D_n.
Weights stay nonnegative, so larger polygons cause no new obligation.
It remains to check mixed quadrilaterals. A D_n short-leaf-pair
deletion has connected white graph A_(n-2), hence cannot lie in a
face with an A_n end deletion. Thus only D_n fork/leaf deletions
can occur in a mixed face.

## 2. Common-white classification: three selected side types

Put N=n-2 in {17,18}, L=N-1, and s=(N+2)/(4N). The white graph of
either relevant deletion is

    W={e,l1,...,lL}=A1+A_L.

Here e is the isolated white root. All roots in what follows are
unit outward normals. At a selected D_n corner orient the black
roots as a,b, with a the fork and b a short leaf. Thus a contacts
e and l1 simply, b is perpendicular to W, and a.b=-1/2.

The following are precisely the black side types participating in
the two selected roles, relative to this fixed W:

    I: no contact with W;
    L: simple contact only with e;
    F: simple contacts with e and ONE endpoint of the long path.

The F placement is NEAR if the endpoint is l1 and FAR if it is lL.
Both placements must be retained. A selected D fork/leaf corner is
exactly F--I with native bond3. A selected A end-distance-two corner
is exactly F--L with native product0: the full path is
L--e--F--(long white path), and the marked roots are its first and
third nodes. These statements follow from the complete whole-vertex
diagrams and do not discard spectators.

Two F sides cannot be adjacent in an actual finite corner. If their
long contacts use the same endpoint, W plus the two F roots contains
a four-cycle through e and that endpoint. If they use opposite
endpoints, it contains the longer cycle through the entire path.
Adding a native bond between the F roots does not remove a cycle.
No finite Coxeter diagram contains a cycle. Consequently a mixed
quad has at most two F sides, and two F sides are opposite.

If an I side is adjacent to an F side, a nonzero finite native bond
can only be3: any larger bond makes a large branched nonfinite
component. Bond3 is the selected D corner; native product0 is an
original zero. If an L side is adjacent to an F side, every nonzero
native bond makes a cycle with e, hence is impossible at a finite
corner. Their product is therefore0 and that corner is selected A.

It follows that every genuinely mixed newly unpaid quad has one of
the following forms, up to cyclic reversal:

    (i)   F--I--R--L,    one F side;
    (ii)  F--I--F--L,    two opposite F sides.

In (i), the F/I corner is selected D and the L/F corner selected A.
The two R corners must be original zeros. In particular R.I=R.L=0,
and R cannot contact e, since that would put R at marked distance2
from L. Complete finite continuation of the large white path leaves
R isolated or attached simply at an endpoint, by a 4-bond at an
endpoint, or simply next to an endpoint. In (ii), both L corners
are selected A, and the second F/I corner is either original zero
or selected D. These are the three mixed patterns to be checked.

There is no missed L--L zero: two L roots are already at distance2
through e, and their whole vertex is not a selected A or D context.
There is likewise no selected D role involving L and a non-F side.

## 3. The common opposite product and its arithmetic interface

Use cyclic sides a,b,c,d, where a is the old F fork, b the old I
leaf and d the L side. Then

    d.e=-1/2,  d perpendicular to every other old root except b,
    a.d=0,     b.d=-beta.

After eliminating W, the unchanged b has square1, d has square3/4,
and their product remains -beta. They are opposite sides of the
actual compact polygon. The acute two-wall projection argument
excludes a positive-definite two-wall Gram, and the positive-null
case contradicts compactness. Hence

    beta>sqrt(3)/2.                                      (1)

The cycle b-a-e-d-b makes 2beta an algebraic integer in the ground
field. Eliminating the old whole D_n from d gives

    q_d=n*beta^2/2+N*beta/2+n/8.

This independent old-plus-new block has q_d>1 physically. At each
nonidentity embedding its rescaled Gram is positive definite, so

    -N/(2n)-sqrt(n+1)/n < beta'
                           < -N/(2n)+sqrt(n+1)/n.         (2)

Here (2) is a determinant-ratio identity in the cyclic quantities
beta^2 and beta, not an unsupported embedding of arbitrary unit
entries. Moreover K=Q(beta): any nonidentity extension embedding
fixing beta would preserve the physical negative determinant ratio.

Finite native beta has degree at most5. The complete list with
phi(2m)<=10, m>=3, is

    m=3,4,5,6,7,8,9,10,11,12,15.

By (1) only 7,8,9,10,11,12,15 remain. All except9 have a positive
nonphysical cosine conjugate, contradicting (2). For9 the conjugate
cos(7pi/9)<-sqrt(2)/2<-7/10, while (2) lies inside (-7/10,-1/5)
at n=19,20. Thus beta is not finite. Native beta=1 gives a feasible
positive-null pair of actual ambient walls and is forbidden by
compactness. Therefore beta>1.

Put A=2beta+1. It is integral and primitive in K. Equation (2) gives

    -2/5 < A' < 3/5                                    (3)

for every other conjugate. For example the upper comparison follows
from (3n-10)^2>100(n+1), and the lower from
(n+5)^2>25(n+1), both at n=19,20. No field-degree truncation of a
unit-normal square root is involved.

## 4. A short exact integral-interval lemma

Suppose an algebraic integer A generates one of the eight fields,

    3<A<19/5,     -2/5<A'<3/5 otherwise.

Then necessarily

    K=Q(sqrt(3)),       A=2+sqrt(3).                      (4)

Indeed degree at least4 gives a nonzero integer norm of absolute
value below (19/5)(3/5)^3=513/625<1. In degree3, if either small
root is negative the absolute norm is below
(19/5)(2/5)(3/5)=114/125<1. If both small roots are positive, the
norm is1 and the integer trace is4, since it lies strictly between
3 and5. Writing the small roots r,t gives r+t=4-A<1, but
r+t>=2sqrt(rt)=2/sqrt(A)>1, a contradiction.

In degree2 the trace is3 or4. Trace3 gives negative norm -1,-2 or-3,
and hence fields Q(sqrt(13)),Q(sqrt(17)),Q(sqrt(21)), none allowed.
Trace4 gives positive norm1 or2; the fields are Q(sqrt(3)) and
Q(sqrt(2)). The latter is absent from the eight. Degree1 gives Q,
also absent. This proves (4). In particular an additional bound
A<7/2 excludes every possibility.

## 5. Preliminary exclusion of finite a.c in the mixed patterns

Write a.c=-x. In case (i), c=R has no e contact and contacts at most
one long white node. We claim x>1. The x=0 cases are independent
positive-index or affine exclusions inside the old D_n plus c:
an isolated c cannot be perpendicular to a positive n-space; a
simple far endpoint makes D_(n+1); the other simple attachments
give an induced affine D tree, and a long 4-bond gives affine B.
These certificates are listed in Section3 of the source-bound
`arithmetic_dn_fork_whisker_quad_attachment_classification_20261019.md`.

For finite x with a single simple long contact, Section6 of that
source excludes all four positions by its cosine-moment argument
and conjugate positivity, BEFORE using any quad corank equation.
Those arguments therefore apply here unchanged. For no long contact,
v=4x^2 is integral with physical v<4 and all other v'<2/N. Its norm
forces v rational; the old-plus-c determinant ratio is then rational
and negative physically, which would force K=Q. For a long 4-bond,
w=2sqrt(2)x is integral. The source's Section7 conjugate intervals
force w=2, apart from two quadratic values which are not native
finite cosines. Again w=2 makes the negative old-plus-c determinant
ratio rational and forces K=Q. Thus these finite cases are excluded
within the eight fields without importing an isolated-d hypothesis.
Parallel x=1 is forbidden by compactness, proving the claim.

In case (ii) with b.c=0, c is a near or far F. The independent old
D_n+c cost is

    q_c=1-s+2N(x+kappa)^2,
    kappa=1-s (near),       kappa=s (far).

The same source's simple-long/simple-e moment rows leave only
m=4 near or m=5 far, even before their quad-integrality tests.
The physical negative determinant also forces K=Q(x). Thus those
finite cases give precisely Q(sqrt(2)) or Q(sqrt(5)), outside the
eight. Zero x has a positive-null principal cycle. Again x>1.

In case (ii) with b.c=-1/2, the actual triangle a-b-c makes 2x
integral. Its old-plus-c cost factors as

    near: q_c-1=2(x+1)(N(x+1)-1),
    far:  q_c-1=2N(x+1/2)(x+1/2+1/N).

In the near case 2(x+1) has other conjugates in (0,2/N); in the
far case 2x+1 has other conjugates in (-2/N,0). For zero or finite
x, physical bounds below4 or3 make a nonzero integer norm impossible
in every nontrivial degree. The independent block forces K=Q(x),
excluding rational x as well. Compactness excludes x=1. Thus x>1
also in this last pattern.

## 6. One F side: the interval upper bound is too small

In case (i), put c_bar^2=t and a_bar.c_bar=-X after eliminating W.
Then t<=1, X>=x>1, and a_bar^2=s. Because c has no e contact,
c_bar.d_bar=0. The four-black Schur matrix is

    [ s    -1/2   -X    -1/4 ]
    [-1/2   1      0    -beta]
    [-X     0      t     0   ]
    [-1/4  -beta   0     3/4 ].

There are n+2 actual roots in Lorentz dimension n+1, so its determinant
is zero. Eliminating c gives the exact equation

    (X^2/t-s)(beta^2-3/4)=(beta+1)/4.                    (5)

The first factor is strictly greater than 1-s>=49/68. The left side
minus the right side increases with beta>=1. At beta=5/4 it is
strictly greater than

    (49/68)(13/16)-9/16=25/1088>0.

Equation (5) therefore forces beta<5/4, hence A<7/2. The integral
interval lemma contradicts (3). This excludes the entire one-F case,
including every original-zero continuation R.

## 7. Two F sides and one original-zero F/I corner

Suppose now b.c=0. Put X=x+1-s near and X=x+s far. In both cases
X>1. The c/d white projection contributes 1/4, so the projected
matrix is

    [ s    -1/2   -X    -1/4 ]
    [-1/2   1      0    -beta]
    [-X     0      s    -1/4 ]
    [-1/4  -beta  -1/4   3/4 ].

For clarity, compute its corank identity through the old D_n block.
Write

    Q=q_d-1=(N+2)beta^2/2+N*beta/2+(N-6)/8.

The other cost and old-span product are

    q_c-1=2N X^2-s,       C=1/4+N X(beta+1/2).

The rank equation C^2=(q_c-1)Q is equivalent to

    F(X,beta):=2N(beta^2-3/4)X^2
               -(N/2)(beta+1/2)X-1/16-sQ=0.             (6)

For X>=1 and beta>=5/4, F is strictly increasing in both variables.
For the X derivative the minimum lower bound is 19N/8>0. For the
beta derivative use s<1/3 and X>=1 to get the positive lower bound
(9N-2)/3. At the lower corner of this range,

    F(1,5/4)=(47N^2-132N-52)/(128N)>0

for N=17,18. Thus (6) forces beta<5/4. Again A<7/2 contradicts
the interval lemma. Both rooted F placements are covered.

## 8. Two F sides and four selected corners

Finally suppose b.c=-1/2. Both F/I corners are selected D and both
F/L corners selected A. The two F roots now have the same projected
contacts with b,d. Splitting their symmetric and antisymmetric
combinations gives the corank equation

    (X-s)(beta^2-3/4)=(beta+1)/2,                        (7)
    X=x+1-s (near),       X=x+s (far).

Near, X-s>x+1-2s>49/34. At beta=5/4 the left side minus the right
side exceeds 25/544>0 and is increasing. Hence A<7/2 and this branch
is impossible.

Far, X-s=x>1. At beta=7/5 the left side of (7) is strictly greater
than121/100, while the right side is6/5. Thus beta<7/5 and
A<19/5. The interval lemma now leaves the one formal possibility

    K=Q(sqrt(3)),       beta=(1+sqrt(3))/2.

Substitution in (7) forces

    x=(beta+1)/(2(beta^2-3/4))=(3+5sqrt(3))/11.

This is a remaining formal real candidate permitted by the interval
tests, rather than an actual arithmetic model. It fails the separate cyclic-integrality
condition: the actual triangle a-b-c has two simple bonds and hence
2x is integral, whereas

    Tr_Q(sqrt(3))/Q(2x)=12/11

is not an integer. The last mixed contact pattern is therefore
impossible too. No numerical solver or graph census is used.

## 9. Conditional composition and remaining global limits

The common-white classification is complete: a genuinely mixed
newly unpaid quad has one F or two opposite F sides, with the second
F/I corner zero or simple. Sections6-8 excluded all these actual
patterns within the exact whole-vertex/eight-field scope.

Consequently, CONDITIONAL ON the two standalone consumer proofs
passing their full reviews, their A and D deletions compose on actual
two-faces. This new lemma does not assume their desired composition
as a hypothesis; it supplies precisely the mixed cases they omit.
There is no spectator extension. B-type reductions have white graphs
containing a4-bond and hence cannot introduce an additional mixed
A/D zero in these faces; their nonnegative coupons only increase
payments.

The resulting connected A20 and D20 source costs would both be35.
This is not a global rank20 allowance35: the existing balanced B20
and B19+A1 profiles still cost141/4, and every multifactor and
conditional charge interface remains subject to its own audited
coverage. No compact arithmetic ceiling19 is claimed.
