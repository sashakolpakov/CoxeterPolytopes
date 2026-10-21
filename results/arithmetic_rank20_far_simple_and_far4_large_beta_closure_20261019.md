# The two remaining far-end quadrilateral rows at beta at least one

19 October 2026. Root-agent bounds and shared-endpoint norm argument,
independently checked and assembled by batch_guard. This is a written
mathematical proof only. No scientific job, enumeration, Lean, TeX or
build was run. Earlier source-bound reports are unchanged.

## 1. Exact local statement

Let the ENTIRE diagram of an actual finite vertex of a compact
arithmetic hyperbolic Coxeter polytope of dimension20 be

    D_r+T_s, r+s=20, 10<=r<=16,

where T_s is D_s,4<=s<=10, or E6,E7,E8. Its ground field is one of

    Q(sqrt3),Q(sqrt6),K9,K11,K15,K16,K20,K24,
    K_m=Q(2cos(2pi/m)).

Write a for the D_r fork, b,e for its short leaves, and
a-l_1-...-l_L for its long arm. Put N=r-2>=8 and L=N-1.
Consider an ACTUAL quadrilateral with side roots a,b,c,d and complete
white support W=A1(e)+A_L+T_s. Assume

    c.b=c.d=a.d=0,
    b.d=-beta, beta>=1,
    d is orthogonal to every old root except b,l_L.

The two cases are

    F: d.l_L=-1/2,
    B: d.l_L=-sqrt2/2.

Both configurations are impossible over the eight stated ground
fields. In particular, this closes the beta>=1 branches left outside
the two earlier finite-beta reports. The condition beta>=1 is an
explicit hypothesis here; it is not inferred just because b,d are
opposite sides in the quadrilateral. Original-zero other corner
weights motivate the displayed orthogonalities but do not replace
them in the theorem.

The statement does not by itself retain every triangle payment or
every face with several subsequently reduced corner weights. It is
not a simultaneous half-refund theorem or a new dimension bound.

## 2. Arithmetic and actual-face inputs

We use the same arithmetic admissibility interfaces as the earlier
rank20 local proofs. Actual cyclic products in the Cartan-normalized
Gram are algebraic integers in the ground field. A whole finite
rank20 vertex and an additional actual wall give a physically
independent Lorentzian block. After projection against the whole
old vertex, its cost is therefore greater than1 physically and less
than1 at every nonphysical embedding. This last inequality is strict:
an independent root set in the positive-definite conjugate ambient
space has positive-definite Gram.

If that projection cost is a rational expression in an actual cyclic
scalar q, any nonphysical ground-field embedding fixing the physical
Q(q) would repeat the physical cost greater than1. Consequently the
ground field equals Q(q). We use this ONLY with the whole rank20
old star, not with an unspecified embedded component. The same
argument applies when the cost itself is rational: the field must
then be Q.

For actual unit wall normals, write x=-a.c. Acuteness and the
Coxeter property give x=0, x=cos(pi/m),m>=3, or x>=1. Thus every
nonzero value is at least1/2. Every white-supported corner mentioned
below is finite because the quadrilateral and its vertices are
actual. These actual incidences supply both the finite continuation
inventory and the corank relation.

An induced affine diagram on selected actual roots is forbidden in
a compact chamber even if those roots have further outside contacts.
Its positive null combination is nonzero, since its product with an
interior future timelike point is negative. Every other actual wall
has nonpositive product with that combination by acuteness. It is
therefore a feasible ideal ray, contrary to compactness.

## 3. A common decreasing positive-root formula

Project off W. The Gram of a,b,c,d becomes

    [ s   -1/2   -X    -rho ]
    [-1/2  1      0    -beta]
    [-X    0      t    -nu  ]
    [-rho -beta  -nu     u  ],

where s=(N+2)/(4N), t,u>0 and X,rho,nu>=0. Since W+c+d is
a finite actual corner, tu-nu^2>0. Actual corank, with
D=beta^2-u>0, gives

    X=nu*(beta/2+rho)/D
      +sqrt[(t+nu^2/D)*Q(beta)],
    Q(beta)=s+(rho*beta+rho^2+u/4)/D.               (1)

This is the exact positive root derived in
`arithmetic_dn_far_fork_quad_obstruction_20261019.md`, SHA256
`5eb859adebf38253472fea1eced2a72acde9e2eb589c13ed5c68787ddbea840c`.
The constant term in the determinant equation has the sign forcing
the positive root because s>1/4. No other root is discarded by a
numerical choice.

For beta>=1, the right side is nonincreasing in beta. Indeed D
increases; the derivatives of (beta/2+rho)/D and
(rho*beta+rho^2+u/4)/D have negative numerators

    -beta^2/2-2rho*beta-u/2,
    -rho*beta^2-rho*u-2(rho^2+u/4)*beta.

All factors are nonnegative. Thus it is legitimate to substitute
beta=1 to obtain an upper bound in every case below.

## 4. Complete finite continuation inventory for F

For F, W+d=A1(e)+A_N+T_s; the path A_N is
l_1-...-l_L-d. Its rank N is at least8, and c.d=0.
The finite connected diagrams containing this whole path allow
precisely the following path contacts by c:

* simple near endpoint l_1, with e contact absent, simple or4;
* near endpoint l_1 by4, with no e contact;
* simple near-next node l_2, with no e contact;
* simple far-next node l_L, with no e contact;
* simple near endpoint l_1 and a long endpoint of spectator D_s,
  with no e contact.

The last alternative is a bridge; E spectators cannot be used there.
In addition, c may have no path contact, with its complete spectator
contact and e contact retained. With no spectator contact its e
bond can be any finite dihedral bond or zero. If there is a
spectator contact, that contact is the single simple contact from
the full finite spectator table, and its e contact is absent or
simple.

Here is why this inventory is complete. A one-root finite extension
of A_N, N>=8, is a longer A path, a terminal4 B path, or D with
the new root attached next to an endpoint. Attachment farther
inside would have exceptional rank at least9. Two path contacts
make a cycle. Extending the short leaf of the D option by e makes
arms1,2,and at least6, which is not finite; putting a4 inside a
long path is not finite. A path-plus-spectator contact can remain
finite only by extending the spectator's D long arm simply. Any
additional e contact, an E spectator, or a D short-leaf contact
would leave two forks or a forbidden exceptional rank. This is a
classification of the ENTIRE W+c+d, not an embedded subdiagram
whose external contacts were suppressed.

For no-path c with spectator contact, let t_T be its residual
square after spectator projection. The full spectator table gives
t_T<=1/2; if c.e=-1/2, its remaining t=t_T-1/4 is at most1/4.
No numerical or field-specific enumeration is needed for this
finite diagram observation.

## 5. F: all contacts except the shared far endpoint

For F,

    u=(N+1)/(2N), rho=1/(2N),
    D(1)=(N-1)/(2N),
    Q(1)=(N+3)/(2(N-1))<=11/14.                   (2)

### 5.1. No path contact

Write gamma=-c.e. Then nu=0, X=x+gamma/2, and (1),(2) give

    (x+gamma/2)^2/t<=11/14.                       (3)

First suppose c has a spectator contact. If gamma=0, t<=1/2.
The value x=0 fails (1); x=1/2 gives the induced affine D4 star
with center a and leaves b,e,l_1,c. Every larger positive native
x has x^2>=1/2, making the left side of(3) at least1. If gamma
is simple, t<=1/4, so x>0 would make it at least9/4. At x=0,
the roots a,b,e,c,l_1,...,l_5 give affine E8 with arms1,2,5 at a.
This excludes every no-path spectator contact.

Now suppose c has no spectator contact. Then t=1-gamma^2.
At gamma=0, x=0 is impossible and (3) leaves only native labels
3,4,5,6 for x. To justify the cutoff exactly,

    cos^2(pi/7)>cos^2(3pi/20)
                =(1+sin(pi/5))/2>11/14.

The last inequality follows from
sin^2(pi/5)=(5-sqrt5)/8>11/32>16/49, using sqrt5<9/4.
Label3 gives affine D4. For labels4,6 the full old-plus-c
projection cost qc=2N*x^2 is rational and greater than1, forcing
K=Q. For label5 it generates Q(sqrt5), forcing equality with
that field. All these ground fields are excluded from the list.

If gamma>0 and x>0, the smallest native pair is x=gamma=1/2.
Any larger native value in either position makes the ratio exceed
11/14: already replacing x by sqrt2/2 gives a ratio greater than1,
as does replacing gamma by sqrt2/2. Thus only the simple/simple
case could remain, but {a,e,c} is affine A2.

It remains x=0,gamma>0. Equation(3) gives gamma^2<=22/29<11/14,
so again only labels3,4,5,6 remain. Label3 gives the affine E8
just described. Label4 gives affine B4 on a,b,l_1,e,c: the
coefficients are2 at a,e,1 at b,l_1 and sqrt2 at c. For labels5
and6, the full old-plus-c cost is r*gamma^2/2, generating
Q(sqrt5) and Q respectively, and the full-star argument excludes
both. Every no-path subcase is closed.

### 5.2. Simple near endpoint

Here gamma is0,1/2 or sqrt2/2, there is no spectator contact, and

    t=(N+1)/(2N)-gamma^2, nu=1/(2N),
    X=x+(N-1)/(2N)+gamma/2.

At beta=1 the first term J in(1) and the first factor P under
the square root satisfy

    J=(N+1)/(2N(N-1))<=9/112,
    P=N/(2(N-1))-gamma^2<=4/7-gamma^2.

Thus X<=9/112+sqrt22/7<15/16. The latter strict inequality
follows even from sqrt22<6. If x>0 the native gap instead gives
X>=15/16. Hence x=0. For gamma=0, the roots a,l_1,b,e,c,l_2
give affine D5, with coefficients2 at a,l_1 and1 at four leaves.
For gamma=1/2, a,e,c,l_1 is an induced simple4-cycle, affine A3.
For gamma=sqrt2/2 the sharper bound is

    X<=9/112+sqrt11/14<7/16,

using sqrt11<5; whereas X>=(N-1)/(2N)>=7/16. Thus that bond
is excluded as well.

### 5.3. Near endpoint4

There is no e or spectator contact. Direct A-path projection gives

    t=1/N, nu=sqrt2/(2N),
    X=x+sqrt2*(N-1)/(2N).

At beta=1, J<=9sqrt2/112 and P=1/(N-1)<=1/7. Consequently

    X<=9sqrt2/112+sqrt11/(7sqrt2)<7sqrt2/16,

where division by sqrt2 reduces the comparison to
9/112+sqrt11/14<7/16. This contradicts the lower bound for X
even when x=0.

### 5.4. Simple near-next node

Here c.l_2=-1/2 with no e or spectator contact. The projected data
are t=2/N, nu=1/N, X=x+(N-2)/(2N). At beta=1,

    J<=9/56, P=2/(N-1)<=2/7,
    X<=9/56+sqrt11/7<7/8.

If x>0, X>=1/2+(N-2)/(2N)>=7/8, so x=0. But then
a,b,e,l_1,l_2,l_3,c give affine D6, with coefficients2 on
a,l_1,l_2 and1 at its four leaves b,e,l_3,c.

### 5.5. Near endpoint and spectator bridge

The spectator is D_s and the contacts are simple. Here

    t=1/(2N), nu=1/(2N), X=x+(N-1)/(2N).

At beta=1, J<=9/112 and P=1/[2(N-1)]<=1/14. Therefore
X<=9/112+sqrt11/14<7/16, contradicting X>=7/16. The
entire bridge is excluded without any spectator-orthogonality
assumption beyond its actual finite continuation.

## 6. F: the remaining shared far endpoint

The last finite continuation has c.l_L=-1/2 and no e or
spectator contact. Here

    t=u=(N+1)/(2N), nu=(N-1)/(2N), X=x+rho.

At beta=1, D=nu, so the first term in(1) is EXACTLY u and
the first square-root factor is EXACTLY1. Equations(1),(2) imply

    x<=1/2+sqrt[(N+3)/(2(N-1))]<3/2.              (4)

Let A=2x. The actual simple cycle from a along the old arm to c
makes A an integral cyclic scalar. Projection against the FULL
old vertex, not just its long path, gives

    q_c=2N*x^2+2x+1/2=N*A^2/2+A+1/2.             (5)

The value x=0 would make q_c=1/2, contradicting its physical
value greater than1. The full-star argument gives K=Q(A).
At every nonphysical embedding, q_c<1 gives

    (-1-sqrt(N+1))/N<A_i<(-1+sqrt(N+1))/N,
    hence -1/2<A_i<1/4.                           (6)

The weaker displayed endpoints hold strictly even at N=8 because
the original quadratic inequality is strict.

If 0<x<1, the native label and K=Q(x) give exactly the finite
cosine-field cases in Section3 of
`arithmetic_rank20_simple_d_finite_beta_and_far_fork_closure_20261019.md`,
SHA256 `51bde71da4e5af431d07c5f52f81b547634fdf389725a6d8568fc7edff711fc0`.
For clarity their complete list is m=3,4,5,6,7,8,9,10,11,12,15.
Labels3,5 force Q,Q(sqrt5); each other label has a nonphysical
cosine outside(-1/2,1/4), whereas (6) requires its cosine in the
smaller interval(-1/4,1/8). Thus every finite x is excluded.

It remains x>=1. By(4), 2<=A_0<3. Put d=[K:Q], with NO
assumption that a proper subfield would be the ground field.
For d>=3 the nonzero integer norm satisfies

    0<|Norm A|<3*(1/2)^(d-1)<=3/4,

impossible. For d=2, a positive other conjugate gives norm
less than3/4. Therefore the other conjugate is negative, and
its norm is an integer in(-3/2,0), necessarily-1. Its trace
is an integer strictly between3/2 and3, hence2. It follows that
A_0=1+sqrt2 and K=Q(sqrt2), excluded from the eight fields.
Degree1 would force K=Q, likewise excluded. This closes F.

## 7. Complete B continuation and its bounds

For B, W+d=A1(e)+B_N+T_s, with the4 bond at the d end.
Since c.d=0 and N>=8, c either has no path contact, or meets
the near endpoint l_1 simply with no spectator contact and e
contact absent or simple. A fork inside B, a second4, or a
path-plus-spectator bridge containing a4 and a fork cannot be
finite in these ranks. This is the same complete finite
continuation argument used in the finite-beta B report; it does
not depend on the value of beta.

For B,

    u=1/N, rho=sqrt2/(2N),
    Q(1)=(N+2+2sqrt2)/(4(N-1))
          <=(5+sqrt2)/14<1/2.                    (7)

For no-path c, (1) gives

    (x+gamma/2)^2/t<1/2.                          (8)

If a spectator is contacted, t<=1/2 for gamma=0 and t<=1/4
for gamma=1/2. The former case excludes x>=1/2, while x=0
fails(1). In the latter, x>0 is impossible and x=0 gives the
same affine E8 witness with arms1,2,5.

If no spectator is contacted, t=1-gamma^2. At gamma=0,
x=0 is impossible and the only possible positive native value
is x=1/2; it gives affine D4 on a,b,e,l_1,c. If gamma>0,
x>0 would make the left side at least3/4. Thus x=0 and
gamma^2<2/3, leaving precisely labels3,4,5. The first two
give the affine E8 and B4 witnesses already displayed. For
label5, q_c=r*gamma^2/2 forces K=Q(sqrt5), excluded.

Finally suppose c.l_1=-1/2. Write gamma=0 or1/2. Then

    t=(N+1)/(2N)-gamma^2, nu=rho,
    X=x+(N-1)/(2N)+gamma/2.

At beta=1,

    J=(sqrt2*N+2)/(4N(N-1))
       <=(8sqrt2+2)/224<1/8,
    P=N/(2(N-1))-gamma^2<=4/7.

Together with(7), this gives X<1/8+sqrt(2/7)<7/8<15/16.
For x>0 the lower bound is15/16, so x=0. The cases gamma=0
and1/2 give the same induced affine D5 and A3 witnesses as
in Section5.2. This closes B for beta>=1.

## 8. What the composition does and does not prove

The finite-beta F result is in the report with SHA256
`51bde71da4e5af431d07c5f52f81b547634fdf389725a6d8568fc7edff711fc0`.
The finite-beta B result, including the complete exceptional
N=8,label12 corner check, is
`arithmetic_rank20_far4_finite_beta_obstruction_20261019.md`, SHA256
`a93a8d7b15cb49a09172a7bc4f874ef87ab84dc0b05ee743439d3f99459683d8`.
Together with the present proof these close F and B for every
native beta allowed by the stated actual quadrilateral hypotheses.
The beta=0 cases already fail the full old-plus-d strict inequality.

Combining these local rows with the separately reviewed I,T,J and
bridge rows is a finite necessary-contact classification for a
quadrilateral whose three other ORIGINAL corner weights vanish.
It does not establish that a single half-refund retains all other
payments, that every triangle has enough remaining payment, or that
two simultaneous half-refunds cannot make a previously paid corner
zero. Those consumer and composition obligations must be checked
before installing any new global rule. No new arithmetic dimension
ceiling is asserted here.
