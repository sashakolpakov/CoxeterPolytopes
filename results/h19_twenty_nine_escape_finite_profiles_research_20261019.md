# The twenty-nine-wall escape: fifteen short finite profiles

19 October2026. **Paper derivation only; independent review pending.**
No profile worker, scientific calculation program, search, census, Lean or
build ran. This concerns the NEW29-wall actual edge, not the earlier
eleven-profile edge. It gives24 necessary short profiles, then17 modulo
the discriminant prime and15 after a modulo7 test. The sole possible
long finite profile is impossible at every height. One retained short
profile has an explicit solution, the separately written thirtieth root.
No complete lattice-reflector census or finite-volume claim follows.

## Exact edge and complete white system

The fixed source is the actual29 certificate
`h19_twenty_nine_wall_seed_20261019.json`, SHA256
`d6e7129224a47621c30737d70c87b713b818e8b4154e0b8de2ea5f3ac49cf031`.
Its full source/coordinate bindings and executed volume consequence are
preserved. Use the notation of
`h19_thirtieth_wall_pell_extension_research_20261019.md`, in particular

    D=1+5phi=phi*b, b=4+phi, u=10+17phi, t=12+19phi,
    s=(2u,-phi*u,t*omega1), 7t²-Du²=1,
    V=(1,-phi/2,(Du/(14t))*omega1), V²=-D/(14t²),
    d=(0,1/2,0), E=V+h*d, h=(82-46phi)/77, E²>0.

The fixed whites are s,old1,old3,...,old8. Thus

    W=A1(s)+A1(old1)+A6(chain7--2--3--4--5--6),

with ALL eight normals of norm2. The edge leaves V by freeing old0 and
first meets old17. No long white from a preceding edge is retained here.

The inherited primitive norm alphabet is {2,2b}. A short/long finite
nonright angle is impossible: if g is its product, then g²=b*q_m with
q_m=4cos²(pi/m) in K=Q(sqrt5). The cyclotomic-degree list gives
m=3,4,5,6,10 (the other quadratic cases8,12 give sqrt2,sqrt3).
The rational norms of b*q_m are19,76,19,171,95, none a square. For
short/short roots the only finite nonright products are-1 and-phi.
These are the same explicit field/norm arguments as in the earlier
initial-edge pencil report; they are not inferred from a sampled ledger.

## All short profiles before congruences

For a canonical norm2 root r=(x,y,z), let beta=-(r,s), a=-(r,old1),
and eta_j=-(r,old(j+1)) on the A6 nodes j=7,2,3,4,5,6.
Their finite Schur condition is exactly

    qW=(beta²+a²)/2 + eta^T*A6^-1*eta < 2.          (1)

Each A1 contact is0,1 orphi. There is at most one nonzero A6 contact:
two unit contacts have energy at least2, achieved by the two endpoints,
and larger contacts only increase energy since the inverse is positive.
At a single A6 position i=1,...,6 the energy is i(7-i)/7, with values
6/7,10/7,12/7 in symmetric pairs. A contact phi is impossible because
(6/7)phi²>2. Thus only the following profiles occur:

|A6 contact|A1 contacts (beta,a)|Count|K=14(2-qW)|
|---|---|---:|---|
|none|(0,0)|1|28|
|none|(1,0),(0,1)|2|21|
|none|(phi,0),(0,phi)|2|7(3-phi)|
|none|(1,1)|1|14|
|none|(1,phi),(phi,1)|2|7(2-phi)|
|endpoint i=1 or6|(0,0)|2|16|
|endpoint|(1,0),(0,1)|4|9|
|endpoint|(1,1)|2|2|
|near-end i=2 or5|(0,0)|2|8|
|near-end|(1,0),(0,1)|4|1|
|middle i=3 or4|(0,0)|2|4|

The count is24 INCLUDING the all-zero profile. The corresponding native
E8 nodes are endpoints7,6; near-ends2,5; middle3,4. Two phi A1 contacts
have energy phi²>2. Any phi A1 together with an A6 contact also exceeds2;
two unit A1 contacts together with a near-end or middle contact do so.
This proves completeness without a new historical core enumeration.

## Integral Pell equation and exact reconstruction

Set eta_0=a, eta_1=0 and keep the selected chain eta_j, writing

    z=T*omega1+sum_{j!=1}eta_j*omega_j,
    R=sum_{j!=1}(C^-1)_{1j}eta_j,
    Q0=eta^T*C^-1*eta, ell=2y+phi*x.

The needed inverse row is(7,14,20,16,12,8,4,10). Since
(r,s)=-Du*x+t*(14T+R)=-beta, completion of the exact norm gives

    7t²ell²-D(x+u*beta)²=t²K,
    K=28-7beta²-14Q0+R²=14(2-qW).                 (2)

The second identity is the ordinary Schur complement removing E8 node1;
the remaining block is A1+A6. All these constants lie in O=Z[phi].

Importantly X=(x+u*beta)/t is itself integral, not just field-valued.
N(t)=11 is prime, N(D)=-19, and reduction of(2) modulo(t) forces
x+u*beta=0 modulo(t). Therefore every integral root has

    7ell²-DX²=K,                                  (3)
    x=tX-u*beta,
    y=(ell-phi*x)/2,
    T=(DuX-7t*beta-R)/14.                         (4)

Conversely an integral solution of(3) with both displayed divisibility
conditions in(4) produces an integral norm2 vector with the selected
white contacts. It still needs every current signed native product and
halfspace test; the Pell equation alone does not certify a chamber wall.

The edge products are

    (V,r)=-DX/(14t), (d,r)=ell/2,
    theta=DX/(7t*ell),
    (V+theta*d)²=-DK/(98t²ell²).                   (5)

For the intended halfspace containing V and cutting E, X and ell are
positive and theta<h. The new endpoint is then future timelike because
K>0. All other29 inequalities are strict there apart from the eight
common whites, by their complete V/E ledger. No rejection based solely
on the old strict center is appropriate: if necessary a sufficiently
small positive mixture V+epsilon*c supplies a common strict center.

## Prime(D) and modulo7: twenty-four to fifteen

At(D)=(b), O/(D)=F19 with phi=15. Equation(3) requires
ell²=11K modulo19. For the eleven rows above the corresponding values are

    4,3,7,2,6,5,4,3,12,11,6.

The nonsquares3,2,3,12 eliminate respectively the single-unit A1 row,
the two-unit A1 row, the endpoint/two-unit row, and the near-end/no-A1
row:2+1+2+2=7 labeled profiles. Thus17 remain.

For either no-A6/single-phi row, R=7a. Integrality of T in(4) forces
X=0 modulo7: Du is invertible there, as N(Du)=361 is prime to7.
Write X=7Y. Dividing(3) by7 gives

    ell²-7DY²=3-phi.

Its reduction modulo7 would make3-phi a square. But
N(3-phi)=5 is a nonsquare modulo7 (the nonzero squares are1,2,4).
Multiplicativity of the norm rules this out at every height. These two
additional profiles are eliminated.

The remaining15 are the zero profile; the two no-A6 mixed1/phi profiles;
two endpoint/no-A1 profiles; four endpoint/one-unit profiles; four
near-end/one-unit profiles; and two middle/no-A1 profiles. No all-height
exclusion of those15 is asserted. In fact one is constructed below.

## No long root can close this edge at a finite vertex

A canonical norm2b root must be orthogonal to every short white, by the
finite-angle obstruction above. This leaves just the zero-contact long
profile. In the preceding notation beta=a=eta=0 and X=x/t is integral
by the same prime(t) argument. Reflection integrality is `Gr/b in O^10`.
The free E8 covector is T*e1 and T=DuX/14, so

    T/b=phi*u*X/14 in O.

Since N(phi*u)=19 is coprime to14, X=14Y with Y in O. Also ell=bL,
and the first binary reflection condition is
`phi*(L-x)/2 in O`. Because x=tX is even, L is even: L=2Z.
The norm2b equation is7ell²-DX²=28b. Substitution gives

    bZ²-7phi*Y²=1.                                (6)

Modulo7, taking the norm gives19*N(Z)²=1, hence N(Z)²=3 modulo7,
impossible. Thus there is NO canonical long finite continuation of this
edge, independently of height, unit rescaling or coefficient size.
This does not exclude a long mirror meeting it beyond the light cone,
a wall that deletes V, or a different completion route.

## An exact retained short solution

Take beta=a=0 and the middle A6 contact at native node4. Then R=12,K=4.
The old identity7t²-Du²=1 gives

    X=2u, ell=2t, T=t²-1.

Reconstruction(4) is integral. It is precisely
`Kroot=t*s-t*old0+(omega4-omega1)` from the separate thirtieth-wall
report. That report supplies all29 signed native products, a strict
center, actual facet proof and finite endpoint, not just this necessary
profile equation. It proves a30-wall seed with at least88 finite
vertices on paper, pending the separate fixed checker. No finite volume
or complete reflection census is inferred.

These fifteen profiles belong to the OLD29 escape. After appending that
root, any further escape must again be analyzed using its own actual
white system. No unchanged profile list is claimed for the next edge.
