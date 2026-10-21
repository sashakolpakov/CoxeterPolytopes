# Uniform outgoing lengths: every center-cutting wall has a volume cost

20 October 2026. **Written conditional theorem; no execution.**
This extends the fixed-central-wall result, without assuming that an
unknown wall is centered on the dodecahedral axis.

Let P be a compact right-angled hyperbolic four-polytope with an actual
regular right-angled dodecahedral facet D. Use outward unit normals and
the Lorentz form with negative time. Put `phi=(1+sqrt(5))/2`. Let n be
the base normal and o its intrinsic future unit center. Assume all twenty
edges leaving D have length `L=arcosh(phi)`. Their opposite actual walls
are therefore

    m_v=sqrt(phi)*v-phi*n,             v a vertex of D.

Define the matched regular-120-cell center

    q0=sqrt(2)*phi*o-phi^(3/2)*n,      q0^2=-1.

**Theorem.** If any actual wall m satisfies `<m,q0>>=0`, its facet has
at least48 vertices, and

    f(P)>=101,   W(P)>=21,   Vol(P)>=23*pi^2/2>34*pi^2/3.

There is at most one such wall. Consequently a polytope in this marked
uniform-length class with volume at most that of the regular120-cell
must contain q0 in its interior. This does **not** prove containment of
the whole120-cell, its universal minimality, or the nonexistence of
larger compact completions.

## 1. Fixed data and normalization

Use the reviewed regular-dodecahedron data

    b_j=h*o+B*r_j,       h^2=phi/2, B^2=1+h^2,
    v=C*o+sqrt(3)*h*w_v, C^2=1+3*h^2, C=phi^2/sqrt(2).

The unit r_j are icosahedral face directions; the unit w_v are the twenty
vertex directions. Their convex hull contains the Euclidean ball of
radius `C/(sqrt(3)*B)`. The five w_v on face j average to that radius
times r_j. An adjacent vertex satisfies
`w=phi*v-sqrt(phi)*b_j` for the released face wall. All these identities
are intrinsic statements about the regular D, not a conjectured list
of the unknown walls of P.

The normals n, twelve b_j and twenty m_v have negative product with q0.
Thus the wall in the theorem is additional. It is not adjacent to D,
since all twelve actual neighbors of D are the b_j. The compact
right-angled normal gap gives `<m,n><-1`; hence uniquely

    m=(z-p*n)/s,  z^2=-1 future, p>1, s=sqrt(p^2-1).

Write `Z=-<z,o>`, `r=sqrt(Z^2-1)`, and
`p0=sqrt(2/phi)=1/h`. The cut condition is exactly

    p>=p0*Z,                       hence p>=p0.           (1)

Throughout, distinct actual walls have product0 or strictly less than-1.

## 2. A cut wall cannot be orthogonal to a vertex roof

Suppose first `<m,m_v>=0`. At v let b_1,b_2,b_3 be its three base face
walls. Write

    z=sqrt(phi)*p*v-sum x_j*b_j,  x_j>=0,
    sum x_j^2=phi*p^2-1,         S=sum x_j.

The cut condition gives `S>=phi*p`. Acuteness with the three adjacent
vertex roofs gives `x_j<=p/phi`. No x_j can vanish, since otherwise
`S<=2p/phi<phi*p`. Thus every x_j is positive and the base-wall gap
gives `x_j>s`. Cauchy also gives `p^2>=3/sqrt(5)`.

For each adjacent roof the gap strengthens its acuteness alternative to

    x_j=p/phi,   or   x_j<(p-s)/phi.

But `(p-s)/phi<s` in this range: the needed inequality is
`(phi^4-1)*p^2>phi^4`, which follows already at `p^2=3/sqrt(5)`
from `3(3phi+1)>(3phi+2)(2phi-1)`. Only the first alternative remains.
The norm equation would then force

    p^2=1/(4phi-6)>phi,

contradicting `p/phi>s`. Hence every m_v is strictly ultraparallel to m:

    cosh d(z,v)>sqrt(phi)*p+s/sqrt(phi)  for all v.      (2)

## 3. Every cut center is very close to o

Average (2) and use (1). Since
`C/p0-sqrt(phi)=1/(2sqrt(phi))`, this gives `s<p/2` and `p^2<4/3`.
In particular

    Z^2<2phi/3<B^2.

The vertex-direction support bound also gives

    min_v cosh d(z,v)<=C*g(Z),
    g(Z)=Z-(h/B)*sqrt(Z^2-1).

The function g decreases on `1<=Z<B` and `g(1)=1`. Combining with (2),

    g(Z)>(sqrt(phi)*p+s/sqrt(phi))/C=:H(p).

It follows that `p<pmax=sqrt(phi)*C-sqrt(3/2)<57/50`. More importantly,

    r<1/20.                                               (3)

Here are small exact bounds for the latter assertion. Use
`809/500<phi<16181/10000`, `sqrt(2)>707/500`, and `C<463/250`.
Then `H(p)>=H(p0)=(sqrt(2)+phi^-2)/C>897/926`. On the other hand,
at `r=1/20` the inequalities `h/B>2/3` and `sqrt(401)<801/40` give
`g(Z)<2323/2400<897/926`. Monotonicity excludes that and every larger r
in the already established range. The pmax bound uses
`sqrt(phi)*C<59/25` and `sqrt(3/2)>61/50`.

For completeness, two cut centers would have
`cosh d(z_1,z_2)<4phi/3-1<2/phi<=p_1*p_2`, using just `p_i^2<4/3`.
Their actual wall product would be positive. Thus there is at most one.

## 4. All neighbors of the cut wall avoid the base face walls

Let M_j be the average of the five prescribed roofs on base face j and
put `gamma=(phi+2)/phi^3`. The exact identity

    b_j=gamma*(M_j-q0/sqrt(phi))+(2h/phi^3)*o

first shows `<m,b_j><0`, hence `<m,b_j><-1`, for all twelve j.

We need the stronger statement that every actual neighbor ell of m also
has `<ell,b_j><0`. Equivalently, it suffices to express b_j modulo m
as a strictly positive combination of the twenty roofs. Starting with

    b_j=gamma*M_j-o/h+gamma*phi*n,

choose mu so

    mu/s=phi*h/(p*sqrt(phi)*C-phi*Z).

The denominator is positive by (1). The residual
`-o/h+gamma*phi*n+mu*m` lies in the linear span of the roofs. It is in
their strictly positive cone provided

    C*((3-phi)*Z-p/h)>sqrt(phi)*B*r.                       (4)

Indeed its total roof coefficient is positive; after division by this
coefficient, (4) places its spatial vertex average strictly inside the
ball `C/(sqrt(3)*B)` contained in `conv{w_v}`. An interior point of that
convex hull has a representation with all twenty coefficients positive.
This also proves strictness when pairing with a nonzero ell orthogonal
to m: the roofs projected to m-perpendicular span that four-space.
Indeed their original span has dimension four, and the positive
denominator above shows that it does not contain m.

The simple bounds `p<57/50`, `r<1/16`, `Z>=1`, `C>9/5`,
`3-phi>69/50`, `1/h<9/8`, and `sqrt(phi)*B<7/4` verify (4):
its left side exceeds `351/2000`, while its right side is below `7/64`.
Consequently

    -<ell,b_j>>=1 for all twelve j.                      (5)

Use the natural intrinsic center

    O_z=(p*z-n)/s.

It is an actual interior point of the facet of m. In fact, for a neighbor
ell, orthogonality to m gives `<ell,O_z>=s*<ell,n><0`: ell is neither n
nor one of the twelve b_j, all of which are disjoint from m. This is
the center used below, **not** the projection of q0.

## 5. A neighbor of height at most33/20 must meet a roof

Let `t=-<ell,O_z>>0` and decompose
`ell=t*O_z+sqrt(1+t^2)*u`, where u is unit and perpendicular to z,n.
Set

    a0=sqrt(phi)*C,   E^2=(11phi+7)/10,   P=p/s.

In the old o-centered spatial coordinates, ell's time is at most
`P*t*Z+r*sqrt(1+t^2)`, while its spatial norm is at least
`sqrt(1+t^2)-P*t*r`. The dodecahedral support bound therefore gives

    max_v <ell,m_v>
    >=-[A(p)+a0*P*(Z-1)+E*P*r]*t+(E-a0*r)*sqrt(1+t^2),
    A(p)=(a0*p-phi)/s.

On the proved range, A(p) decreases and
`A(p)<=phi^(3/2)<103/50`. Also `P<23/10`, `a0<59/25`, and
`157/100<E<8/5`. Use the weaker consequences `r<1/18` and
`Z-1<1/648` of (3). The last display is bounded below by

    -L*t+R*sqrt(1+t^2),
    L=368197/162000,       R=259/180.

This decreases in t. At `t=33/20`, use `sqrt(1+t^2)>77/40` to get
a value greater than `-1058717/1080000>-1`. Therefore, if `t<=33/20`,
the compact actual-wall gap forces `<ell,m_v>=0` for some v.

## 6. The vertex frame excludes that low height

Fix such a v and its three b_j. The frame
`q_v=phi*v-sqrt(phi)*n,b_1,b_2,b_3` is orthogonal in m_v-perpendicular.
By (5), write

    ell=a*q_v-sum y_j*b_j,  y_j>=1,
    sum y_j^2=1+a^2.

There is no off-center error in the height comparison:

    <ell,z>=p*<ell,n>,   <ell,O_z>=s*<ell,n>,
    t=s*sqrt(phi)*a.

Since `s>=phi^-3/2`, the assumed height gives
`0<a<=(33/20)*phi<267/100`. The three adjacent roof products give the
exact alternatives

    y_j=a/sqrt(phi),   or   y_j<a/sqrt(phi)-1/phi.         (6)

Two values of the first type are impossible: the remaining square would
be `1-a^2/phi^3<1`. No first-type value is also impossible. In that case
the three lower bounds force `a>=phi^(3/2)>2`, whereas the norm and (6)
would imply

    0<=3*(a/sqrt(phi)-1/phi)^2-a^2-1
      <(7/8)*a^2-(14/5)*a+1/5<0   on [2,27/10].

The last convex quadratic is negative at both endpoints. Thus exactly
one value, say y_1, is `a/sqrt(phi)`. Put S=sum y_j. Then

    a>=phi^(3/2),    y_2^2+y_3^2=1+a^2/phi^2,
    S/a>=9/5.                                               (7)

Here is a complete small inequality for the last assertion. Write
`b=phi^-1/2`, `d=phi^-1`, `U=b*a-d`, and `T=y_2+y_3`. We have
`2<=T`, `U<2`, and
`y_2^2+y_3^2<=U^2+(T-U)^2`. If `S<=9a/5`, monotonicity in T bounds
the latter expression using `T=(9/5-b)*a`. Subtracting `1+d^2*a^2`
gives

    Q=(56/25+6d-36b/5)*a^2+d*(18/5-6b)*a+2d^2-1.

The rational bounds `309/500<d<6181/10000` and `b>7861/10000` show

    Q<(7217/25000)*a^2-(1725147/2500000)*a
      -11795239/50000000.

This convex polynomial is below `-2/5` at a=2 and below `-1/50` at
a=267/100. It is therefore negative throughout the required interval,
contradicting the norm equation. This proves (7).

Finally set `A=-<ell,o>`, `N=-<ell,n>`, and let W_ell be ell's old
Euclidean spatial component. Then

    A=a*phi*C-h*S>0,      N=a*sqrt(phi),
    A/a<=J=phi*C-(9/5)*h<139/100,
    ||W_ell||^2/a^2<=phi^-3+J^2-phi<9/16.

Positivity of A follows directly by summing (5); alternatively (6)
gives `A/a>=phi*C-3/sqrt(2)=sqrt(2)/phi>0`. The J bound uses
`phi*C<3` and `h>179/200`. Combining (1) and (3),

    s*<m,ell>>=p*(N-h*A)-r*||W_ell||
               >=p0*(N-h*A)-r*||W_ell||>0.

For an explicit strict margin,

    p0*(N-h*A)/a
    >=sqrt(2)*((9/10)*sqrt(phi)-(1/2)*sqrt(5))
    >75649/2000000>3/80,
    r*||W_ell||/a<3/80.

The radical lower bound uses `sqrt(phi)>159/125`,
`sqrt(5)<22361/10000`, and `sqrt(2)>707/500`. This contradicts
`<m,ell>=0`. Every intrinsic wall of the facet of m therefore has
height strictly greater than33/20 at O_z.

## 7. The actual normal fan and the volume bound

At each vertex of this compact right-angled three-facet, the three
intrinsic normals are orthogonal. Their Euclidean Klein directions
give a triangle in its complete normal fan. If all three heights exceed
33/20, the determinant/half-area identity gives

    tan(Omega/2)<100*sqrt(3667)/(1189*sqrt(1489))<2/15,
    Omega<4/15.

The exact last comparison is
`562500*3667=2062687500<2105030569=1189^2*1489`.
The full fan covers area4pi, so `V>15pi>46`. A compact right-angled
three-polyhedron is trivalent, hence V is even and `V>=48`.
The inherited inequalities now give

    f(P)>=2V+5>=101,
    W(P)>=ceil(3(V-20)/4)>=21,
    Vol(P)=pi^2*(f(P)+W(P)+16)/12>=23*pi^2/2.

This is the asserted conditional exclusion below the regular120-cell
volume. The absence of a cut wall implies only that q0 is interior to P.

## Source dependencies and review boundary

Fixed dodecahedral metric, actual-wall conventions, and the normal-fan
identity: `results/compact_ra4_central_ball_completion_volume_cost_20261020.md`,
SHA256 `df5e75b3cfd8466b9b35f36b44f6c6b37f94b25afcdb276e7a443039b1d6da87`.
The special central wall in that theorem is **not** an assumption here.

Compact wall gap, facet capacity, neighboring excess and volume identity:
`manuscripts/compact_ra4_volume_constraints.tex`, SHA256
`afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`.
The universal positive-excess gap is not needed.

All new off-center steps are written above. Root supplied the
no-orthogonal-roof lead and fan threshold; the finite57 lane derived the
uniform cone and height proof. The resource lane independently rederived
the final gap classification, rational quadratic envelope and product
contradiction before this source was assembled. Formal source-bound
reviews remain separate from these mathematical discussions.
