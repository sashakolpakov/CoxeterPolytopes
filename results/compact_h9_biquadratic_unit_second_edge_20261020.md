# The second quartic edge: three mirrors, one simple wall

20 October 2026. Paper proof only; no worker, height search or execution.
No compact dimension-nine example or nonreflectivity claim.

Import the arithmetic, primitive norm pool, full B9 zero germ and first
actual wall from `compact_h9_biquadratic_unit_diagonal_first_edge_20261020.md`,
SHA256 `05b6473eeadb4befa89c82af2f7002c7441ada3d1878c93d4630355a42685efe`.
Write

    s=sqrt2, phi^2=phi+1, u=phi-1, alpha=1+s,
    k=s-1, d=2-s, a=u+s, N9=<-a>+I9.

Embeddings are numbered as in that source: (s,phi),(-s,phi),
(s,1-phi),(-s,1-phi). The physical bilinear form is
`-a*x*x' + sum(y_i*y_i')`. All normals have their inward halfspace
written as a nonpositive scalar product.

## 1. A finite upper and the complete six roles

At the actual B8+A1(r) vertex, release beta2 and retain the whole B7
on coordinates3,...,9 together with `r=(1;phi,0^8)`. The edge is

    P(t)=(1;a/phi,t,0^7), t>0.

The integral primitive normal

    c=(s*alpha*phi; alpha*phi^2, phi*(phi+s), 0^7)       (1)

has norm2. Its first spatial coordinate is a unit, and norm2 divides2,
so it is primitive and reflective. Direct products are

    c.r=-1, c.beta1=-s, c.beta2=-phi*(phi+s),
    c.beta_i=0 for i>=3.

It hits at `T=a/(a+1)`. Indeed `c.P(t)=-a*phi+t*phi*(phi+s)`.
The exact inequality

    d*(a+1)^2-a*phi^2=k^2/phi^2>0                   (2)

puts this hit before light. The plane Gram matrix for r,c is
`[[d,-1],[-1,2]]`, determinant `k^2`, and has native angle pi/8.
These facts give an upper only; firstness is proved below.

Here is the full finite simple-extension inventory. Set
`b=-(v,r)` and let epsilon be the old B7 projection square.

|role|q=v^2|b|epsilon|
|---|---|---|---|
|isolated|1|0|0|
|isolated|2|0|0|
|isolated|d|0|0|
|B8 long-end extension|2|0|1|
|I2(8) at r|2|1|0|
|I2(4) at r|d|k|0|

For clarity, the added golden subfield does not add an attachment to r.
If a native finite angle is pi/m, its squared cosine field is contained
in K. The elementary cyclotomic degree/subfield check leaves only
`m=2,3,4,5,6,8,10`. For q1 or q2, the squareclasses of q*d and
3*q*d exclude labels3,4,5,6 as appropriate: their rational norms in
Q(s) are nonsquares, and hence they are neither squares nor five times
squares there. In a quadratic extension F(sqrt5), an element of F is
a square only in those two cases. For q=d, label3 would require
`d/2`, label5 `d*phi/2`, both nonintegral at (s); label6 requires
sqrt3, outside K. For label10 the relative norm to Q(s) of
`q*d*(phi+2)/4` is five times a square, and is not a square in Q(s).
For label8, the three required products squared are respectively
`1/2,1,d/2` for q1,q2,qd: only q2 gives an integral product, namely1.
For label4, qd gives product k, while the q1 and q2 products fail the
same squareclass test. Thus the only native nonzero contacts with r
are exactly the last two table rows. The nonnative angle3pi/8 is not
an additional simple extension.

The standard finite real Coxeter classification is the other explicit
input to this inventory. Any connected finite extension of the whole
B7 is B8, attaching one norm2 long root with projection square1.
An 8-bond cannot belong to a finite component of rank>=3. The only
new root that could meet r by4 has norm d; it cannot meet an old B7
root except by8 with a long root, again impossible in rank>=3.
This proves all six rows, including absence of a bridge or spectator
contact in either r-contact row.

The other actual zero-height wall beta1 would hit at t=a/phi>T.
Every remaining possible actual closing wall has time x>0, since the
initial full zero-height chamber is B9. Thus it suffices to test all
six rows on the entire closed interval0<t<=T.

## 2. One binary equation and exact bounds

For a tested normal write its first two spatial coordinates as y1,y.
The r product gives `y1=(a*x-b)/phi`. Its other retained projection
has square epsilon. Put

    x=alpha*phi*z, w=z+u*b/s,
    B=a*(2+s)=2phi+s*phi^2, P=2u+s,
    Q=q-epsilon-b^2/d.

Then its norm and the collision bound are exactly

    y^2-B*w^2=Q,             y>=P*w>0.               (3)

For example, before completing the square the norm is
`y^2-B*z^2-2(phi+s)*b*z+u^2*b^2=q-epsilon`.
The cap inequality follows directly by substituting t<=T in
`t=(a*d*x/phi^2+a*b/phi^2)/y`.

The identities

    P^2-B=d*u^4, 1/B=s*phi/2-1,
    d/B=-phi-2+s*phi^2                            (4)

give the complete four-place bounds

    w_0^2<=Q_0*phi^4/d_0,
    w_i^2<=Q_i/(-B_i), i=1,2,3.                  (5)

All three nonphysical B_i are negative. Nothing here assumes a field
square for a candidate norm equation or restricts a chosen height.

For b0 and residual S=q-epsilon=1, summing (5) gives

    Tr(z^2)<=5+3phi+s*(2phi+1)<16.                 (6)

For either contact row, Q=d/2. Put H=s*w=s*z+u*b, an integral
element satisfying `H==u mod(s)` because b is1 or k. For the isolated
qd row put H=z. Its norm equation modulo(s) also forces
`H==u mod(s)`: y is divisible by s, and division by d gives
`H^2==d/B==phi mod(s)` in F4. Both situations therefore satisfy

    H_0^2<=phi^4, H_i^2<=d_i/(-B_i),
    Tr(H^2)<=10+2phi+s*phi^2<18, H==u mod(s).      (7)

For the isolated q2 row, `y^2-B*z^2=2` forces z and y divisible
by s: otherwise the B*z^2 term has odd valuation1 but y^2 has
even valuation. Then x and y1 are divisible by s too, so this
isolated normal is imprimitive. This removes that row at every height.

## 3. Complete small-trace filters

For h=A+Bphi+s(C+Dphi), use the trace form imported from the first
source. Its binary part is

    U(A,B)=4(A+B/2)^2+5B^2.

The only values at most16 are0,4,6,14,16, represented by
`0; +/-1; +/-phi,+/-u; +/-phi^2,+/-(2-phi); +/-2`.
This proves the following finite lists without a coefficient box.

For (6), the complete positive-physical list, up to overall sign, is

    1,phi,u,phi^2,2-phi;
    s,s+1,s-1,s+phi,phi-s,s+u,s-u;
    s*phi,s*u.

The nonphysical bounds in (5) imply

    |z_1|<3/2, |z_2|<5/4, |z_3|<4/5.

They exclude every entry, as recorded explicitly here:

|entries|violated embedding|
|---|---|
|1,u,2-phi|3|
|phi,phi^2|1|
|s,s+1|2|
|s-1,phi-s,s-u|1|
|s+phi,s+u|3|
|s*phi|1|
|s*u|2|

For (7), the even integral trace is at most16. The congruence modulo
s leaves only

    H in {u,phi^2,s+u,s-u}, up to overall sign.

The bounds imply `|H_2|<1` and `|H_3|<3/2`. They exclude u and
s-u by embedding2 and s+u by embedding3. Since H_0>0, the unique
remaining value is

    H=phi^2.                                      (8)

Here are elementary rational seals for every strict enclosure just
used. The bounds `21/13<phi<13/8` and `24/17<s<10/7` follow by
squaring rational bounds for sqrt5 and sqrt2. They give the right
side of (6) less than `893/56<16` and that of (7) less than17.
Furthermore

    -B_1=s*phi^2-2phi >1302/2873>4/9,
    -B_2=2u-s*u^2 >62/91>16/25,
    -B_3=2u+s*u^2 >213/125>25/16.

For (7), `-B_2-d=u*(s-2u)>0`; also
`9*(-B_3)-4*(2+s)=18u-8+s*(9u^2-4)>12/7>0`.
These prove the two asserted H bounds. Thus no decimal or unlisted
square test is part of either finite filter.

## 4. All simultaneous roots and the full face

At (8), `d+B*phi^4=2*[phi*(phi+s)]^2`. Reconstruction in the
three surviving rows yields exactly

|role|normal|
|---|---|
|q2,b1|c|
|qd,bk|r+c|
|qd,b0|alpha*r+s*c|

All meet at t=T. The last two are primitive mirror normals but not
simple chamber walls. Indeed r is actual, and the oriented integral
mirror c has its nonpositive halfspace containing the full chamber.
Their positive combinations vanish only on their common codimension2
intersection. Consequently c is the unique first SIMPLE wall; no
accepted simultaneous mirror was discarded before testing firstness.

This proves the actual whole corner B7+I2(8). Its plane is unimodular
over O, since the Gram determinant k^2 is a unit. One exact perpendicular
time vector is

    T0=alpha*phi^2*(phi+s)*(1;a/phi,a/(phi+s)),
    T0^2=-a, det[T0,r,c]=k.

Thus the original lattice splits integrally as

    N9=(O*T0+I7) orth (O*r+O*c),

with the first summand literally N7. This identifies the FULL induced
mirror arrangement, not only its orthogonal subarrangement. An ambient
mirror crossing the perpendicular hyperbolic space spans, with r,c, a
positive-definite rank3 space. Its integral reflection group is discrete
and physically compact, hence finite. A finite irreducible reflection
group of rank>=3 has no pair of mirrors at angle pi/8. The crossing
mirror must therefore be orthogonal to the plane. Conversely every
intrinsic N7 reflector extends by identity on the unimodular plane.
The actual I2(8) face is exactly the full N7 chamber.

For a precise low-rank calibration, the four N2 normals beta1,beta2,r,c
bound the affine quadrilateral

    0<=y2<=y1<=a/phi, c1*y1+c2*y2<=a*c0.

Its vertices are `(1;0,0)`, `(1;a/phi,0)`, `(1;a/phi,T)` and
`(1;1,1)`. Their squares are respectively `-a`, `-a*d/phi^2`,
`-a*k^2/[phi^4*(a+1)^2]` and `2-a`, all negative. Every nonincident
wall inequality is strict. It is therefore a compact Coxeter
quadrilateral with angles4,2,8,4. The full integral chamber refines
this compact quadrilateral; this proves N2 compact-reflective without
claiming the four-wall polygon is already its full chamber.

The complete edge and recursive face identity consequently give no
nonreflectivity or compactness conclusion for N9. No field-specific
dimension refinement or host closure is asserted.
