# Complete finite interval for the next L59 edge

20 October 2026. **Paper proof and source-preparation specification.** No
new checker was executed or imported. The explicit finite arithmetic below
selects the long upper mirror; its separately scheduled exact replay is
still pending. This is not a compactness claim.

Use the L59 notation and native E8 basis of the reviewed route note
`compact_h9_golden_discriminant59_edge_route_20261020.md` (SHA256
`725d8775cebec51d8c18b4b0a0d3be0255e8b1472621adb4d3cefa656adb429e`).
The initial sixteen-root certificate has actual generation and separate
replay PASS; it is inherited here, not rerun. Put

    delta=5+12phi=phi*b,   b=7+5phi,
    h=(3+2phi,1+2phi;0),  h^2=-2delta,
    c=(1,1;0),           (h,c)=0.

At the actual E8+A1(c) vertex h, release native node 0. The retained
white diagram is D7 on nodes 1,...,7, together with c. Its edge is

    P(tau)=h+tau*omega0,       omega0^2=4,

and its light endpoint has tau=sqrt(delta/2).

## The eight finite profiles

The D7 arms have lengths 1,1,4. A new short root can extend its long
arm to D8, extend either fork leaf to E8, or join the long arm to c by
two simple bonds to make D9. Otherwise it is orthogonal to D7 and has
with c angle pi/2, pi/3 or pi/5. A long root must be orthogonal to all
short whites by the whole-arrangement prime-59 isolation theorem. No
other finite spherical extension is possible: a finite connected
component containing D7 is D8, E8 or D9, and both relevant rank-nine
white components are already accounted for. Equal short norms allow
only simple or golden bonds; a golden bond cannot belong to a component
containing D7.

Write a new root's binary part as (X/2)h-(k/2)c, where X=x-y is in O.
If it contacts node j of D7 by -1, write its E8 part as

    z=((Z-m_j)/4)*omega0+omega_j,
    m_j=(omega0,omega_j).

If it is D7-orthogonal, write z=(Z/4)omega0. Native E8 inverse identities
give m_7=2, m_1=5, m_2=7 and residual projection costs

    kappa_7=1,      kappa_1=kappa_2=7/4.

In every row its norm equation and contact parameter are

    Z^2-2delta*X^2=R,       tau=delta*X/Z,       X>0, Z>0.

| New component | k | Contact in D7 | R |
|---|---:|---|---:|
| isolated short | 0 | none | 8 |
| A2 on c | 1 | none | 6 |
| I2(5) on c | phi | none | 6-2phi |
| isolated long | 0 | none | 8b |
| D8 | 0 | node 7 | 4 |
| E8 | 0 | node 1 | 1 |
| E8 | 0 | node 2 | 1 |
| D9 bridge | 1 | node 7 | 2 |

The classification is about actual simple walls at a first finite
endpoint. It does not discard arbitrary non-simple mirror angles merely
because they are not fundamental dihedral angles.

## Four congruence exclusions and exact reconstruction

Modulo b, phi=34 and 2delta=0. The residues of R=8,6,6-2phi,2 are
8,6,56,2, all nonsquares in F59: 2 is a nonsquare, 3 is a square and
-1 is a nonsquare. This excludes the isolated-short, A2, I2(5), D9 rows
at all heights. The four literal surviving rows are D8, the two E8
placements, and the isolated long root.

All survivors have k=0. Since h_y=1+2phi is a unit, integral binary
coordinates force X in 2O. For D8, integral E8 covectors force
Z=4T+2 with T in O. For the two E8 rows they force respectively
Z=4T+5 and Z=4T+7. The E8 Cartan matrix is unimodular, so these are
equivalent to, not merely necessary for, integral spatial coordinates.
For the long row, reflection integrality gives Z/(4b) in O. Reducing
its norm equation divided by b then gives

    X^2=45 mod b,     X=a+d*phi => a+34d=+/-24 mod59.

These conditions are applied before square testing. Full root norm and
reflection integrality are checked again after reconstruction.

## A proved small finite upper mirror

There is an explicit integral primitive long root

    w=(38+61phi)h+(133+213phi)omega0,
    w^2=2b,      133+213phi=b(9+14phi).

Indeed Gh=(-delta,delta;0), so Gw/b is integral. Direct norm expansion
gives w^2=14+10phi. Primitivity also follows from the squarefree norm
ideal, or from the coprime norms 41 and 649 of its two scalar
coefficients. It is orthogonal to all eight retained short whites.

Its profile variables are

    Xv=76+122phi,       Zv=532+852phi,
    Zv^2-2delta*Xv^2=8b.

Thus tau_v=delta*Xv/Zv is strictly before the light endpoint. The actual
edge must meet a new simple wall by this finite point. The upper mirror
need not itself be first, and its compatibility with a finite cap list
is never a filter in the first-contact census.

For completeness, a larger E8 upper was found first without any unit
library. Put lambda^2=2delta. Then

    (7+lambda)/(phi^4-lambda)
      =-[(15+23phi)+(1+4phi)lambda]

is integral of relative norm -phi^-5. Squaring and multiplying by phi^5
gives a norm-one unit

    Z*=16715+27044phi,     X*=2392+3870phi.

The identities

    Z*+1=4b(9+14phi)^2,
    X*=(9+14phi)Xv

exhibit the smaller long bisector above. In fact s_w(h)=Z*h+delta*X*omega0,
and s_w(-e0) is the larger fork-at-node-2 E8 root. Consequently that
particular E8 upper cannot be the first endpoint. No independent
chamber symmetry is inferred from this reflected pair.

## Exhaustive tiny coefficient domain and its twenty entries

For any row, contact by tau_v is equivalent to

    8b*X^2 <= R*Xv^2.

The conjugate equation gives

    0 <= R'+2delta'*X'^2,

and hence |X'|<1/2, 1, 13/5 for R=1,4,8b respectively. The physical
bounds are X<25, X<50, X<=Xv<274, respectively. For the first two,
use Xv^2<5000b; explicitly 5000b-Xv^2=14340-8428phi>0. The last
bound follows from phi<81/50. These are conservative bounds; the exact
inequalities above remain the actual acceptance domain.

The exact divisibilities make an even smaller complete census possible.
For D8 put X=2Y and Z=2W. Its equation is

    W^2=1+2delta*Y^2,       0<Y<25,
    (-2delta')Y'^2<=1.

For either E8 placement the identical square equation applies with Y=X;
its bounds are also 0<X<25 and (-2delta')X'^2<=1, in addition to X in
2O. We may therefore exclude both rows using the larger common Y list.
Since -2delta'>4, every such Y=a+d*phi has |Y'|<1/2. Consequently
0<=d<=12. If m=floor((sqrt(5)-1)d/2), all possible a occur in
m-1,...,m+2. The integer formula

    m=(isqrt(5*d*d)-d)//2

is exact for d>=0. Testing these 52 integer pairs by the displayed exact
inequalities gives the following ten values and norms. The norm column
is obtained by expanding E=1+2delta*Y^2=A+Bphi and using
N(E)=A^2+AB-B^2.

| Y | N(1+2delta*Y^2) |
|---|---:|
| 1+phi | 99 |
| 1+2phi | 641 |
| 2+3phi | 2059 |
| 3+5phi | 5771 |
| 4+6phi | 5401 |
| 4+7phi | 5571 |
| 5+8phi | 15489 |
| 6+9phi | 1531 |
| 6+10phi | 20249 |
| 7+11phi | 24131 |

None is an integer square. A square in K has rational norm a rational
square, so all D8 and both E8 profiles are excluded. Using the larger
common list is deliberate: neither the E8 parity nor the different
E8 congruences need be used to obtain this exclusion.

For the long row write the root as A*h+b*t*omega0. The earlier
integrality statements give A,t in O. Its norm equation is

    2b*t^2-phi*A^2=1,
    A^2=(-4+14phi)t^2+1-phi.

The interval bound A<=38+61phi is equivalent, by the physical equation,
to 0<t<=9+14phi. The positive conjugate equation gives
2b'*t'^2<=1. In particular |t'|<9/25, since phi<81/50 gives b'>39/10
and (81/625)(39/5)>1. Thus t=a+d*phi has 0<=d<=16: the conservative
upper follows from t<32, |t'|<9/25 and sqrt(5)>2. The same four a
values around m give 68 exact tests. They leave precisely these ten
entries; the norms again use the displayed integral quadratic formula.

| t | N((-4+14phi)t^2+1-phi) |
|---|---:|
| 1+2phi | 305 |
| 2+3phi | 1181 |
| 3+5phi | 3475 |
| 4+6phi | 1895 |
| 4+7phi | 1189 |
| 5+8phi | 9481 |
| 6+10phi | 11071 |
| 7+11phi | 12659 |
| 8+13phi | 25205 |
| 9+14phi | 1681 |

Only the last norm is square, 1681=41^2. Its radicand is
5165+8357phi=(38+61phi)^2. This has just two square roots in K;
positive contact selects A=38+61phi. It reconstructs exactly w above.
Thus the complete finite-profile interval has one candidate, not merely
one candidate compatible with a postselected finite wall cap.

There is no implicit bound in height, relative-unit exponent or a
floating comparison. The two finite lists are consequences of the
proved interval and conjugate bounds, not an exploratory prefix.

For an integral E=A+Bphi, exact field-square testing needs only integer
square roots. If E=Z^2, N(E)=n^2 is an integer square and N(Z)=+/-n.
For each sign, t^2=2A+B+2N(Z) must be an integer square, while
(t^2-4N(Z))/5 must be the square of the phi coefficient of Z. Checking
both signs and the parity of the resulting two coefficients, followed
by squaring back to E, is an exhaustive square test in O. Integral
closure ensures that a K-square root of an algebraic integer lies in O.

## Selection, the full germ and remaining scope

The interval has only the primitive long root w, so there is no
simultaneous-contact or dominance ambiguity. Its norm and reflection
integrality, all sixteen old wall products and native angles, endpoint
feasibility, and positive definite corner D7+c+w are separate exact
interface checks. In particular postselection compatibility cannot
create the uniqueness assertion. The new future apex can be written

    H=(10+18phi)h+(38+61phi)omega0,
    H^2=-4/phi,        (H,w)=(H,c)=(H,alpha_i)=0 (i=1,...,7).

Since the full chamber edge cannot cross an omitted wall before the
proved integral upper, the complete finite-profile inventory makes w
the first simple wall and gives the full finite germ D7+c+w. The
integral reflection arrangement is locally finite, and the upper is
strictly inside hyperbolic space; hence the first contact is an
ordinary finite corner. A finite-volume or compactness assumption on
the entire chamber is not needed for this bounded edge argument.

The normal pool, spherical-extension completeness and the first-contact
geometry are explicit written premises. A saved program result will
certify the finite arithmetic domain and these exact interface checks.
The combined worker also checks the fixed integral isometry columns in
`compact_h9_golden59_long_facet_period_20261020.md` (SHA256
`699e9922fc0cc3ca6f16d5619e02ab345fbcb39ed25ddb355d32a4120ef14f9d`),
its Gram preservation, frame images and trace. That note supplies the
separate chamber-germ and infinite-coset argument; polynomial equality
checking alone does not prove it. Actual worker execution is PENDING.
No new compact polytope or finite-volume example is claimed here.
