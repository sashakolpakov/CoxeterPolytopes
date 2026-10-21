# A second Pell mirror at the thirty-wall H19 escape

19 October2026. **New paper construction; exact checker/review pending.**
No scientific program, field/matrix calculation program, root search,
census, Lean process or build ran in this authoring lane. The30-wall
source and existing evidence are unchanged. This proves a proposed31-wall
Coxeter seed with at least89 finite vertices; volume is unresolved.

## Actual current escape

Use the fixed form and roots of
`h19_thirtieth_wall_pell_extension_research_20261019.md`, SHA256
`e62ee4140d132700e8e036d6e4b692943e707b6a121994ae583583debedaf688`.
The coordinator actually generated the30-wall certificate
`h19_thirty_wall_seed_20261019.json`,50,559 bytes, SHA256
`d308740335fabc8cae80224670d777d24d8bcfb2c3484fec4b0ed0bb48758176`.
At this paper handoff, its first separate replay attempt was memory-gated
before launch; that historical deferral is not represented as a PASS.
The generation ledger is09f06707. The fixed source is2f7cb082 and must
remain unchanged.

Its walk stopped after27 finite vertices,23 queued, on an actual feasible
spacelike endpoint from vertex3=p_s, freeing **old1** and first hitting
old17. Put

    phi²=phi+1, b=4+phi, D=1+5phi,
    u=10+17phi, t=12+19phi, 7t²-Du²=1,
    s=(2u,-phi*u,t*omega1),
    v=omega4-omega1=(-1,-2,-2,-1,0,0,0,-1),
    K=t*s-t*old0+v.

With omega_j=C_E8^-1 e_j, that normalized starting vertex and direction
are

    V=(1,-phi/2,(Du/(14t))*omega1), V²=-D/(14t²),
    d1=(0,0;e0/2)=-old1/2, d1²=1/2, (V,d1)=0,
    E1=V+h1*d1, h1=(-46+36phi)/77,
    E1²=15767/11858-499sqrt(5)/847>0.

The common white system is old0,s,old3,...,old8, again2A1+A6 but with
the OTHER binary/simple A1. All30 wall products at E1 are nonpositive.
This note uses this actual new edge, not the previous unchanged profile
list.

## Root, sparse covector and integrality

The coordinator proposed replacing old0 by old1 in the construction:

    L=t*s-t*old1+v=K+t*(old0-old1).                  (1)

Here old1=-e0. It has norm2 and is orthogonal to s and v. Since
v²=2 and(v,s)=-2t,

    L²=2,
    L=(2tu,-phi*tu;(t²-1)*omega1+omega4+t*e0),
    GL=(-Du*t,0;2t*e0+(t²-t-1)*e1+e4).             (2)

The coordinates and covector are integral in O=Z[phi], so the reflection
`x -> x-(x,L)L` is integral. Evenness of the lattice and norm2 imply
primitivity, as in the preceding K proof. Its first coordinate remains
2tu=886+1434phi. Moreover

    (L,K)=2-2t²<-2,

so BOTH Pell mirrors coexist, unlike the earlier mutually positive pair
of bounded-search alternatives.

## Exhaustive thirty-wall native product ledger

Each row follows from the already written exact K products and

    (L,H)=(K,H)+t*((old0,H)-(old1,H)).               (3)

For a native positive-height wall H, the correction in(3) is
`t*((Cz_H)_0-phi*x_H-2y_H)`. The prefix coordinates and existing sparse
covectors give every entry below. Thus no scalar product program is
required for this paper derivation.

|Existing wall|Exact product with L|
|---|---|
|old0,old3,old4,old6,old7,old8,s|0|
|old1|-2t|
|old2|-t²+t+1|
|old5|-1|
|old9|-t(19+28phi)|
|old11|-t(19+29phi)-phi|
|old12|-t(9+13phi)-phi|
|old13|-t(10+18phi)-1|
|old14|-t(10+19phi)-1-phi|
|old15|-t(27+46phi)|
|old16|-t(10+17phi)-1-phi|
|old17|-3t*phi-1-2phi|
|old18|-t(13+20phi)-1-2phi|
|old19|-t(13+22phi)-1-2phi|
|r|-t(95+152phi)|
|k|-t(122+197phi)-2phi|
|w|-t(168+270phi)|
|a|-t(46+74phi)|
|d-root|-t(205+333phi)-1-2phi|
|z|-t(69+117phi)-5-6phi|
|q|-t(3+7phi)-2-2phi|
|h-root|-t(330+536phi)-u|
|j|-t(58+95phi)-5-11phi|
|K|2-2t²|

These are seven right pairs, one angle3 pair,22 ultraparallel pairs.
Every nonzero product except-1 is less than-2 for norm2 old roots;
against r,w,z its magnitude is greater than2sqrt(b)<5. In particular
the old2 and K entries are strictly ultraparallel because t>12.
All previous435 Coxeter pairs are unchanged, giving465 pairs in total.

## Common center and actual facets

For the same strict future timelike center c used in the preceding
checks, (c,old0)=-1/5 and(c,old1)=-1/10. Equations(1) and the checked K
center product give

    (c,L)=-[t*(807+1311phi)+7]/10<0.                 (4)

All31 inequalities are strict at c. Acute orthogonal projection of c
onto each wall is timelike and strict on every other wall, exactly as
in the K proof. Thus all31 inequalities are actual facets; this is not
merely an abstract compatible Coxeter diagram.

## New finite endpoint and every required retention

The new edge products are(V,L)=-Du/(7t) and(d1,L)=t. Hence

    theta=Du/(7t²), P_L=V+theta*d1,
    P_L²=-D/(98t^4)<0.                             (5)

As before theta lies before the light-cone parameter, while E1 is
spacelike, so0<theta<h1. L cuts E1 strictly. The complete V/E1
ledger gives strict inequalities at P_L for every old wall outside
the eight common whites. In particular K is strict there: its product
is constant along this edge, since(K,old1)=0. The active graph at P_L
is E7+2A1, with isolated old0,s and branch old5 with arms
`L`, `old6--old7`, `old4--old3--old8`. It is an actual finite vertex.

For p,p_s,p_q,p_h the old1 inequality is active. Since
`v=old1+2old2+2old3+old4+old8`, their s/v products prove(L,P)<0
as in the K proof; at p_s use its strictly negative old2 product.
The additional P_K=V+theta*(-old0/2) has(L,P_K)=-Du/(7t)<0.

Only p_j needs a separate estimate because old1 is not active there.
Its exact spatial coordinates are

    z_j=A*omega0+B*omega1,
    A=5b/38, B=(-3+4phi)/19.

Consequently

    (p_j,s)=-5(29+31phi)/38,
    (p_j,old1)=-5b/38, (p_j,v)=-(8+21phi)/38,
    (p_j,L)=-[t*(125+150phi)+8+21phi]/38<0.

Thus all six previously named finite vertices survive STRICTLY, retaining
their old nine active walls. P_L is new because its L inequality is0.
The ledger gives L strictly ultraparallel to both r and w; the strict
cap projection witnesses fix the preserved halfspace orientation on
their entire hyperplanes. The two full disjoint41-vertex compact H8
facets therefore survive. All seven named vertices are strictly off
both caps, for the new one by the strict V/E1 cap products. This gives
at least82+7=89 actual finite vertices, with the cap counts inherited
rather than re-enumerated.

A separate small fixed31 checker should independently verify this full
ledger and one bounded edge walk. It may find another escape; absence
of an early escape is not a finite-volume theorem. No unrestricted root
census, compact H9 example, whole-lattice reflectivity or dimension
improvement is claimed here.
