# Finite opposites for simple first-d contacts, and the complete far-fork row

19 October 2026. Root-agent finite-opposite lead and earlier exact
Schur formula, independently checked and extended to N>=8 by
batch_guard. Paper only: no scientific job, enumeration, calibration,
Lean, TeX or build was run. Existing source-bound reports are unchanged.

## 1. Scope and statement

The entire old compact finite vertex is D_r+T_s, r+s=20,
10<=r<=16, where T_s is D_s (4<=s<=10) or E6,E7,E8.
Let a be the D_r fork, b,e its short leaves, and
a-l_1-...-l_L its long arm. Put N=r-2>=8,L=N-1. The white
support of an ACTUAL quadrilateral a,b,c,d is

    W=A1(e)+A_L+T_s.

Assume c.b=c.d=0 and that d is orthogonal to every old root
except b and the single stated long-arm node. The two simple
first-d types considered here are

    F: d.l_L=-1/2,
    T: d.l_(L-1)=-1/2.

Write b.d=-beta and a.c=-x. Ambient opposite products are not
assumed dotted. They obey the actual acute Coxeter-wall condition:
each is zero, a finite native cosine cos(pi/m),m>=3, or at least1.

Over the eight target ground fields

    Q(sqrt3),Q(sqrt6),K9,K11,K15,K16,K20,K24,

the finite-beta branch 0<=beta<1 is impossible for BOTH F and T.
For T the entire branch beta>=1 is also impossible, so the T row
is completely excluded. The F beta>=1 and the separate far4 B
rows are not asserted settled here.

## 2. Full old vertex plus d and the actual ground field

Both simple d types have an actual closed path from b through
a and the old long arm to d, with every other bond simple.
Thus 2beta is an integral cyclic product and beta belongs to K.
The FULL rank20 old vertex, including T_s, plus the extra actual
d wall is independent and physically Lorentzian: d is strict at
the simple old vertex. Its projection cost therefore exceeds1
physically and is less than1 at every nonphysical embedding.

The unit D_r inverse gives

    q_F=r*beta^2/2+beta+1/2,
    q_T=r*beta^2/2+2beta+1.                         (1)

For example the simple next-to-far node is standard D_r node2:
its inverse diagonal is4, its product with the leaf inverse is2,
and the leaf diagonal is r/2. These give q_T directly. The far
endpoint is standard node1 and gives q_F. The spectator is
orthogonal to d, so it contributes no omitted projection cost.

If Q(beta) were a proper subfield of K, a nonphysical K embedding fixing
the physical Q(beta) would repeat the full negative old-plus-d
block, or equivalently the value q_d>1. This contradicts strict
nonphysical positivity. Thus

    K=Q(beta).                                    (2)

The contradiction excludes such an extension. This is not an
assumption that the visible cyclic field of an arbitrary embedded
component equals the ground field: the whole old rank20 star is used.

At beta=0 equation(1) has value1/2 or1, contrary to q_d>1.
For every nonphysical embedding, q_F<1 forces

    (-1-sqrt(r+1))/r<beta_i<(-1+sqrt(r+1))/r,

which is contained in(-1/2,1/4) for r>=10. This follows by
squaring sqrt(r+1)<=r/2-1 and <=r/4+1, valid for r>=8.
For T, q_T<1 gives

    -4/r<beta_i<0,

also contained in(-1/2,1/4). Thus BOTH types require the common
strict window

    -1/2<beta_i<1/4                               (3)

at every nonphysical embedding.

## 3. Every finite native beta is excluded in the target fields

Suppose beta=cos(pi/m),m>=3. By(2), the degree is
phi(2m)/2, and it is at most5 in the eight target fields.
The complete possible list is

    m=3,4,5,6,7,8,9,10,11,12,15.

For completeness, if phi(n)<=10 with n=2m even, every prime
factor has p-1<=10, so p is2,3,5,7 or11. The allowed exponents
are at most4,2,1,1,1 respectively. Multiplying their Euler factors
then leaves exactly n=6,8,10,12,14,16,18,20,22,24,30 for n>=6.
This is a finite elementary totient check, not a numerical search.

The native conjugates are cos(k*pi/m) with gcd(k,2m)=1,
modulo the usual k and -k identification. For all cases except
m=3,5 the following single nonphysical conjugate violates(3):

|m|k|forbidden conjugate|
|---|---|---|
|4|3|cos(3pi/4)<-1/2|
|6|5|cos(5pi/6)<-1/2|
|7|5|cos(5pi/7)<-1/2|
|8|3|cos(3pi/8)>1/4|
|9|7|cos(7pi/9)<-1/2|
|10|3|cos(3pi/10)>1/2|
|11|3|cos(3pi/11)>1/2|
|12|5|cos(5pi/12)>1/4|
|15|11|cos(11pi/15)<-1/2|

All displayed k are coprime to2m and differ from +/-1, so
these really are nonphysical conjugates. The negative comparisons
follow by comparing their angles with2pi/3. For the positive
comparisons, cos(3pi/8)>cos(5pi/12), and
cos(5pi/12)=(sqrt6-sqrt2)/4>1/4; the last inequality follows
from sqrt6>1+sqrt2, equivalently3>2sqrt2. No floating-point
comparison is used.

For m=3, equation(2) would force K=Q; for m=5 it would force
K=Q(sqrt5). Neither is one of the eight target fields. In particular
a quartic field merely containing Q(sqrt5) is not an exception:
the full-star argument proves equality in(2). This excludes the
entire finite-beta branch for both F and T before examining c.

## 4. Complete c inventory for T when beta>=1

Now d is type T. The finite opposite corner has

    W+d=A1(e)+D_N+T_s,

where D_N has fork l_(L-1), leaves l_L,d, and long endpoint l_1.
Since N>=8, its only one-root finite extension is a simple bond
to that long endpoint. Short-leaf extensions would exceed E8,
internal contacts create two forks or degree4, and two contacts
create a cycle. The condition c.d=0 is kept.

A c simultaneously joining D_N and T_s would connect two forks
in rank19, which is not a finite diagram. Therefore exactly these
types remain:

* c meets l_1 simply, has no spectator contact, and c.e is0 or-1/2;
* c has no long-arm contact and no spectator contact, and
  c.e=-gamma with gamma zero or any finite native cosine;
* c has no long-arm contact and meets the spectator once; its
  e contact is absent or simple.

The spectator one-root table in
`arithmetic_rank20_fork_spectator_continuation_research_20261019.md`,
SHA256 `95b0102731bb12a43416197260d65ff2998cc0a92db189722bee8ad43bf2d754`,
keeps the entire T_s. Its residual square is at most1/2, or at
most1/4 if the independent simple e contact is present. No white
contact has been silently sealed just because it was not in an
embedded D subdiagram.

## 5. The exact positive-root bound, valid for N>=8

Project off W and write the four-root Gram in order a,b,c,d as

    [ s   -1/2  -X    -rho ]
    [-1/2  1     0    -beta]
    [-X    0     t    -nu  ]
    [-rho -beta -nu    u   ],

where

    s=(N+2)/(4N),       u=2/N,       rho=1/N.

The finite opposite corner gives tu-nu^2>0. Actual corank and
the positive root of the determinant equation give, with D=beta^2-u,

    X=nu*(beta/2+rho)/D
      +sqrt[(t+nu^2/D)*(s+(rho*beta+rho^2+u/4)/D)]. (4)

This formula and its derivation are in
`arithmetic_dn_far_fork_quad_obstruction_20261019.md`, SHA256
`5eb859adebf38253472fea1eced2a72acde9e2eb589c13ed5c68787ddbea840c`.
Only that exact identity is reused; its n=19,20 numerical
inequalities are NOT imported to the present smaller r.

For beta>=1, D>0. Each term in(4) is nonincreasing in beta:
for example the derivatives of (beta/2+rho)/D and
(rho*beta+rho^2+u/4)/D have strictly negative numerators,
respectively -beta^2/2-2rho*beta-u/2 and
-rho*beta^2-rho*u-2(rho^2+u/4)beta. Thus evaluate the upper
bound at beta=1. The last factor becomes

    Q_N=(N+6)/(4(N-2))<=7/12.                      (5)

## 6. A c contact with l_1

Here gamma=0 or1/2 and

    t=(N+1)/(2N)-gamma^2,
    nu=1/N,
    X=x+(N-1)/(2N)+gamma/2.

The beta=1 form of(4) is

    X<=F_N+sqrt[(P_N-gamma^2)*Q_N],
    F_N=(N+2)/(2N(N-2))<=5/48,
    P_N=(N-1)/(2(N-2))<=7/12.

The bound for F_N follows from (N-8)(5N+6)>=0; the other
bounds follow immediately from N>=8. Therefore X<=11/16.
If x>0, the native gap x>=1/2 gives

    X>=1/2+(N-1)/(2N)>=15/16,

a contradiction. Hence x=0. If gamma=0 the selected roots
{a,l_1,b,e,c,l_2} form affine D5, with coefficients2 at a,l_1
and1 at its four leaves. If gamma=1/2, {a,e,c,l_1} form an
induced simple4-cycle, affine A3, with all coefficients1.
Both contradict compactness.

## 7. No long contact, but a spectator contact

Here nu=0 and X=x+gamma/2, so(4),(5) imply

    (x+gamma/2)^2/t<=7/12.                         (6)

If gamma=0, the spectator residual has t<=1/2. The value x=0
is incompatible with the strictly positive right side in(4).
The value x=1/2 gives the affine D4 star centered at a, with
leaves b,e,l_1,c. Every larger native nonzero x is at least
cos(pi/4), so x^2/t>=1, contrary to(6).

If gamma=1/2, then t<=1/4. A nonzero x gives the left side
of(6) at least9/4, impossible. At x=0 the selected roots
{a,b,e,c,l_1,...,l_5} form affine E8, with arms1,2,5 at a.
Thus the spectator-contact subcase is completely excluded.

## 8. No long or spectator contact

Now t=1-gamma^2 and (6) applies with this t. If gamma=0,
x=0 is impossible, x=1/2 gives affine D4, and the only remaining
native value is x=cos(pi/4): the next squared value
cos^2(pi/5)=(3+sqrt5)/8 exceeds7/12. But at x=cos(pi/4)
the full old-plus-c cost is qc=2Nx^2=N, rational and greater
than1. Nonphysical positivity therefore forces K=Q, excluded.

If gamma>=1/2 and x>=1/2, the left side of(6) is at least3/4.
Hence x=0. Solving(6) then gives

    gamma^2<=(N+6)/(2N+4)<=7/10.

Since cos^2(pi/6)=3/4, the only native labels are3,4,5.
Label3 gives affine E8 with arms1,2,5; label4 gives affine B4
on {a,b,l_1,e,c}, with coefficients2 at a,e,1 at b,l_1,
and sqrt2 at c. For label5, the full old-plus-c cost is

    qc=r*gamma^2/2,

which generates Q(sqrt5). Just as in Section2, any proper
extension fixing gamma^2 would repeat the full physically negative
old-plus-c star at a nonphysical embedding. Consequently K must
equal Q(sqrt5), again excluded from the eight-field set.

All affine witnesses here are induced on the stated selected
actual roots. Their positive null sums are nonzero because they
pair negatively with an interior point. Every other outward wall
pairs nonpositively with that sum by acuteness, including spectator
walls not in the witness. Thus a feasible ideal ray would exist,
contrary to compactness; outside contacts do not undo the argument.

## 9. Precisely what is closed

Sections2--3 exclude finite beta for BOTH simple first-d rows F,T.
Sections4--8 exclude beta>=1 for T with the COMPLETE c spectator
inventory. Thus T is closed over all eight target fields, whereas
F with beta>=1 and the distinct far4 B row remain outside this
theorem. The argument includes Q(sqrt6), which needs no special
native-angle table once the full-star field equality is imposed.

This is a theorem about the stated actual quadrilateral and
original-zero orthogonalities. It is not a proof that every other
modified corner retains payment, not a simultaneous half-refund
rule, and not a new rank20 allowance or dimension bound.
