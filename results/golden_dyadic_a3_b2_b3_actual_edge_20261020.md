# A complete edge leaving the forced-roof nineteen-wall escape

20 October 2026. Exact finite paper proof; no worker or prefix search.
The conclusion is one new ACTUAL wall, not compactness of R5.

Import the full M6/H3 splitting and actual walls from the frozen
`results/golden_dyadic_d4_research_20261020_freeze.json`, and the
unit-roof proof and finite19-wall escape from
`results/golden_dyadic_unit_vertex_roofs_20261020.md`,
SHA256 `e0c6b78e7c9c9ece0b77ef75cdfdcdc54d25c5177b1b474a3bcbead9515653e6`.
Use coordinates `(V0;e1,...,e6;a3,s,r)` throughout, with the original
R5 coordinates stated explicitly when useful. Put delta=phi^-3.

## 1. The actual starting vertex, edge and finite upper

The seven common actual white roots are

    A3(beta1,beta2,beta3)+B2(beta5,beta6)+A2(a3,s).

Adding beta4 and the actual roof a2 gives the full positive definite
corner B6+B3. Hence the actual positive acute projection identifies
the finite vertex

    Q0=V0-2phi^-2*omega3, Q0^2=-2delta.

Release beta4, retaining the whole rank8 system

    A3(beta1,beta2,beta3)+B2(beta5,beta6)+B3(a2,a3,s).   (1)

The B3 has norms(2,4,4), with successive bonds4,3. Its Gram is
`[[2,-2,0],[-2,4,-2],[0,-2,4]]`, determinant8. In particular its
inverse entries at the two ends are3/2 and1/2.

Let j=(0;1,1,1,1,0,0) in M6, with j^2=4 and j.Q0=0. The actual
outgoing ray is `Q0+t*j`, t>0. An explicit primitive integral norm2
reflector is

    N=(phi^3;phi,phi,phi,phi,0,0)-2phi*omega3
     =(phi^2;phi,phi,phi,0,0;-phi,0,0,-phi)             (2)

in the two respective coordinate systems. It has N^2=2 and is
orthogonal to ALL eight roots(1). Its product on the ray is

    N.(Q0+t*j)=-2+4phi*t,

so it supplies the finite upper `t0=1/(2phi)`. Indeed
`t0^2=phi^-2/4<delta/2`, the light parameter squared.
Its unit I-coordinates prove primitivity; its D4 sum is-2phi in2O.
Norm2 makes reflection integrality immediate. The normal is only
an upper until the following complete firstness proof.

## 2. All finite profiles at the whole retained rank8 system

Let a possible first new simple root have canonical norm q1,2 or4,
M6 time a, and coordinates x1,...,x6. Actual-wall acuteness at the
old V0 implies a>=0. Define its nonnegative contacts by

    beta_i: single contact d at A3 node k (k0 means none),
    epsilon=1 for the B2 long-tail extension, otherwise0,
    f=-n.a2, h=-n.s.

There is no contact at the B3 middle node a3. The COMPLETE inventory
has the following17 literal rows (paired ends count separately):

|role|q|k,d|epsilon|f,h|R|count|
|---|---:|---|---:|---|---|---:|
|isolated|1|0,0|0|0,0|4|1|
|isolated|2|0,0|0|0,0|8|1|
|isolated|4|0,0|0|0,0|16|1|
|A3 to A4|2|1 or3,1|0|0,0|5|2|
|A3 to D4|2|2,1|0|0,0|4|1|
|A3 to H4|2|1 or3,phi|0|0,0|phi^-4|2|
|A3 to B4|1|1 or3,1|0|0,0|1|2|
|A3 to C4|4|1 or3,2|0|0,0|4|2|
|B2 to B3|2|0,0|1|0,0|4|1|
|B3 to F4|2|0,0|0|1,0|2|1|
|B3 to B4|4|0,0|0|0,2|8|1|
|A3/new/B2 bridge|2|1 or3,1|1|0,0|1|2|

Here are the finite-class completeness details. A3 admits the finite
rank4 end extensions A4,B4,C4,H4, and the middle extension D4.
The B2 norms(2,1) allow only a norm2 long-tail extension with product1.
Other native lengths either give a second4-bond or a nonintegral
pairing. The B3 norms(2,4,4) allow a norm2 short-end extension to F4
or norm4 long-end extension to B4. Contact at its middle creates a
fork with a4-bond, and no finite diagram has that form. A5-bond at
either end gives the nonfinite4,3,5 chain or4,5 subchain.

The only finite join between two retained components is the norm2
end-to-end A3/new/B2 bridge, giving B6. A3/new/B3 would extend F4
beyond rank4 or leave an interior4-bond; B2/new/B3 has two4-bonds
and rank6. Neither is finite. Joining all three is likewise impossible.
Multiple contacts within one component give a cycle. Thus arbitrary
spectators and every finite native length assignment are covered.
The simple-angle alphabet and integral pairings are the inherited
golden norm-pool interface, not arbitrary unweighted graph extensions.

## 3. Exact reconstruction and the complete finite interval

Put `y=x1+x2+x3+x4`. For k>0 the first k coordinates are z+d and
the other4-k are z, where `z=(y-k*d)/4`. For k0 all four equal y/4.
The B2 coordinates are `(x5,x6)=(epsilon,0)`.

Write the P3 component as `p=-h*omega2-v*omega3`. Its a2 product
gives

    v=2phi^-2*a-2delta*f-2phi^-1*h.

Consequently

    n.Q0=-delta*U, U=2a+3f+h,
    y^2-2delta*U^2=R,
    R=4q-k(4-k)d^2-4epsilon-6f^2-2h^2-4fh.          (3)

One way to verify(3) independently is to project onto the edge plane:
the projection is `(U/2)Q0+(y/4)j`, of square
`(y^2-2delta U^2)/4`. The positive white projection has square
`k(4-k)d^2/4+epsilon+(3/2)f^2+(1/2)h^2+fh`, using the displayed
B3 inverse. This also verifies every R in the table.

A first new wall has U>0 and y>0. Contact before or at t0 gives

    y>=2phi^-2*U,
    U^2<=R*phi^6/2,
    U'^2<=R'*delta/2.                              (4)

The last inequality follows since delta'=-phi^3. These are complete
finite contact bounds, not a chosen coefficient height. Ordinary
integrality gives a and all x_i in O. A primitive norm4 reflector
further has a,x_i in2O because the M6 summand is unimodular.

### Profiles orthogonal to the retained B3

For f=h0, U=2a and(4) becomes

    a^2<=R*phi^6/8, a'^2<=R'*delta/8.              (5)

For integer R<=16 these imply `0<a<6`, `|a'|<7/10`.
For a=A+Bphi, therefore B0,1,2: the upper B bound follows from
`(67/10)^2<45`. B0 gives no positive a, B1 permits A0 or1, and
B2 only A1. The complete largest list is `(phi,phi^2,phi^3)`.
Substitution into the EXACT bounds(5) gives

|R|complete a list|
|---|---|
|1,2,4|empty|
|5|phi^2|
|8|phi^2,phi^3|
|16|phi,phi^2,phi^3|

For example a=phi needs R>=8phi, a=phi^2 needs R>=8/phi, and
a=phi^3 needs R>=8, in addition to their other embedding bound.
For R=phi^-4, the trace bound is
`Tr(a^2)<=(phi^2+phi)/8<2`, excluding every nonzero a.

R5 is the A4 row. At a=phi^2 equation(3) gives y=phi^3; its two
end reconstructions have z=phi/2 or(phi-1)/2, both nonintegral.
R8 is isolated q2. At a=phi^2, y^2=8phi^2 is nonsquare because2
is inert. At a=phi^3 it gives y=4phi and z=phi, exactly N(2).
R16 is isolated q4; its three a-values are units and fail a in2O.
Every other B3-orthogonal row is already empty in the list.

### The two retained-B3 contact profiles

For the F4 row, f1,h0,q2,R2, the complete U-strip is

    0<U<=phi^3, |U'|<=sqrt(delta)<1/2,
    U in{phi^2,phi^3}.

This last two-element list follows directly by writing U=A+Bphi:
B0 gives no positive integer, B1 gives A1, and B2 gives A1.
Since U=2a+3 and a>=0, only U=phi^3 remains. It gives
`a=phi^-1,y=2phi`, hence z=phi/2, again nonintegral.

For the B4 row, f0,h2,q4,R8, put `v0=U/2=a+1`. Its strip is
identical to the preceding U-strip, so v0 is phi^2 or phi^3.
The q4 condition a in2O forces `v0=1 mod2`, selecting phi^3.
Thus a=2phi,y=4phi and z=phi, which violates the necessary q4
condition x_i in2O. This row is empty too.

All17 literal profiles have now been exhausted. No smaller field
square, sign choice or simultaneous simple wall is unaccounted for.

## 4. Actual endpoint and remaining scope

The only survivor is the primitive norm2 N(2), and it meets precisely
at the finite upper. Therefore N is the UNIQUE first simple wall.
The endpoint has full finite diagram

    A3+B2+B3+A1(N).

Its positive definite whole Gram and the unique first wall give the
actual finite endpoint, not merely a model with the same local types.

For additional coordinate sanity, N's nonzero products with all19
known normals are

    beta4=u=-phi; w=-4phi-1; h=-phi^7; r=-2phi;
    A(a2)=-38phi-24; b=-7phi-4;
    b0=-2phi^2; A(b0)=-4phi^6;
    W=-4phi^3; Z=-9phi-5.

All other products, including N.a2, vanish. These are native5-bonds
or strict dotted products; N cuts the recorded19-wall escape by4phi.
This proof adds one mathematically forced ACTUAL facet. It does not
assert that twenty known halfspaces close compactly, that the retained
white system is sealed, or that a diagram permutation extends to an
integral chamber period. Those are separate global questions.
