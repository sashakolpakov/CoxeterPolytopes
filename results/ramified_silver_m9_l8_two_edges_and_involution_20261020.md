# Two complete L8 edges and an actual M9 chamber involution

20 October 2026. Paper only. This gives a finite-order chamber symmetry,
NOT an infinite period or a closure of the remaining literal M9 host.

Use s=sqrt2, alpha=1+s, k=alpha^-1=s-1. Import the complete actual
two-edge/facet proof
`compact_h9_ramified_silver_two_actual_edges_20261020.md`, SHA256
`72a99612ee2f2c00c842ba1d73f1a2f3278d9ab2dccc6f135c1342548827fc66`.
It gives the ACTUAL facet of M9 on r=(1;s,0^8), and the exact integral
orthogonal decomposition

    M9=L8 orth <r>, L8=<-alpha>+I8, r^2=2-s.

The intrinsic time vector o maps to alpha*(1;1,0^8) in M9, and its
eight coordinate vectors are the remaining original coordinate units.
The intrinsic primitive norm pool is 1,2,X=2+s, by unimodularity and
the unit-square normalization; its full zero germ at o is B8. Put
beta_i=-e_i+e_(i+1), i1,...,7, and beta8=-e8.

## 1. The first A7 edge: a full D8 vertex

Release beta8 and retain A7=(beta1,...,beta7). The norm2 integral
primitive reflector

    g=(2alpha;(alpha+1)^2,alpha^6)

is a finite upper. With A=5+4s and b=alpha^2 its coordinate sum is
2A, and the edge endpoint is

    P=(A;b^8), P^2=-alpha.

The pairings are g.beta2=-1, all other retained products0. The old
beta8 product is -alpha. The primitive check follows already from
the unit difference between its two spatial levels.

Here is the COMPLETE weighted finite inventory for a first wall at
an A7 face. For a connected profile write the spatial coordinates
as `((y+lambda)^j,y^(8-j))`, and put Z=8y+j*lambda:

|role|norm|j|lambda|residual R in Z^2-8alpha*x^2=R|
|---|---:|---|---|---:|
|isolated|1,2,X|0|0|8,16,8X|
|A8|2|1,7|1|9|
|D8|2|2,6|1|4|
|E8|2|3,5|1|1|
|B8 short|1|1,7|1|1|
|B8 equal-norm|2|1,7|s|2|

This is thirteen literal roles. Finite connected diagrams containing
the whole A7 are exactly A8,B8,D8,E8. The B8 equal-norm normalization
is essential; its product-s with a norm2 endpoint is a4-bond. A normX
new root cannot make a finite higher-rank connection. Multiple contacts
create a cycle, so none is omitted.

Slope comparison with g and conjugate positivity imply

    x>0, x^2<=R*alpha^2, x'^2<=R'*alpha/8.       (1)

For the connected rows the exact complete coefficient lists are:

|R|all possible positive x|
|---|---|
|1|alpha|
|2|alpha, 2+s|
|4|1,alpha,2+s,2alpha|
|9|1,s,alpha,2+s,3+s,2alpha,alpha^2,4+2s,3alpha|

To check completeness write x=m+ns. The bounds on
`n=(x-x')/(2s)` give respectively n0..1, n0..1, n0..2, n0..3;
negative n is impossible. Intersecting the remaining m intervals gives
the table. Safe enclosures are respectively
`x<5/2,|x'|<11/20`; `x<7/2,|x'|<4/5`;
`x<29/6,|x'|<11/10`; `x<29/4,|x'|<33/20`.
The exact inequalities in (1) retain the stated endpoints.

For reproducible HAND arithmetic, if x=m+ns set
C=m^2+2n^2+4mn and N=m^2-2n^2. Then

    Norm(R+8alpha*x^2)=R^2+16R*C-64N^2.

In the displayed x order this gives

|R|rational norms|
|---|---|
|1|49|
|2|164,196|
|4|16,400,656,784|
|9|161,113,1025,1841,257,3089,5921,4049,3969|

A nonsquare rational norm excludes a K-square. At R4,x=alpha, the
square target is 60+40s with norm400: a putative integral square has
constant-coordinate square20 or40, so this case also fails. All other
square hits have the following exact positive Z:

    R1,x=alpha:        Z=A;
    R2,x=2+s:          Z=8+5s;
    R4,x=1:            Z=2+2s;
    R4,x=2alpha:       Z=2A;
    R9,x=3alpha:       Z=3A.

Old-coordinate integrality is EXACTLY `Z=j*lambda mod8O`.
The first three hits fail this in every corresponding role. The last
R9 hit has s coefficient12, also failing. Only R4,x2alpha survives,
and Z2A is congruent to2, not6, modulo8. It gives g and j2.

The isolated rows need no finite list. Their native equation is
`8y^2-alpha*x^2=q`. Modulo2O, alpha*x^2 is only0 or1+s, while
q1 and qX require1 or s, so both fail. At q2 the same reduction
forces the integer coefficient of x=m+ns to be even. Reducing the
s coefficient modulo4 then forces n even. Thus x is in2O, which
contradicts `alpha*x^2=-2 mod4O`. All isolated roles are absent at
every height. Therefore g is the unique actual first simple wall,
and P has the FULL actual D8 germ (g,beta1,...,beta7).

## 2. Its other fork edge: another full D8 vertex

Release beta1, retain W=(g,beta2,...,beta7), an A7 in that order.
Set v=beta1-g, v^2=4, and let T be reflection in v. T is a
K-isometry; it is NOT an integral lattice isometry. It fixes P and
beta2,...,beta7 and exchanges beta1 with g. In particular it carries
the old A7 to W. Put d0=-e7-e8 and

    h=T(d0)=(2alpha^2;A,b^5,(2alpha)^2), h^2=2,
    K=T(o)=(15+10s;13+9s,(7+5s)^7), K^2=-alpha.

These displayed vectors ARE integral, despite T's nonintegrality.
The h coordinate difference is a unit, so h is primitive. It has
product-1 only with beta6 in W, its node6.

For exact interval coordinates take

    R0=(6alpha;2+3s,(4+3s)^7),
    R0^2=8, R0.P=0, R0.beta1=2,
    K=A*P-b*R0.

The forward edge is P-tR0. The h products are
`P.h=-2b`, `R0.h=-2A`, so it hits at t=b/A, before light.

Let a candidate have a=-(P,z)/alpha in O and Y=-(R0,z)>0.
Writing its W projection with any role in the SAME thirteen-role
table gives `Y^2-8alpha*a^2=R`. Its finite upper h gives exactly
the bounds (1), now for a. Thus the SAME complete coefficient and
field-square tables above apply; there is no second hidden census.

For clarity here is the FULL new glue calculation. Apply T to the
candidate and write its old-A7 coordinates as
`(x;(y+lambda)^j,y^(8-j))`, with spatial sum S. Then

    x=A*a-alpha*Y, S=8alpha^2*a-A*Y,
    (Tz,v)=[Y-(8-j)lambda]/4.

Substituting into z=T(Tz) proves that all original coordinates are
integral if and only if

    Y=(8-j)lambda mod4O,
    (4+3s)Y=j*s*lambda mod8O.                    (2)

For instance coordinates2,...,8 have common fractional part
`[j*s*lambda-(4+3s)Y]/8`; the extra first-coordinate difference
is `[Y-(8-j)lambda]/4`. The time congruence is then automatic.
This proves necessity and sufficiency, rather than checking selected
coordinates only.

At the five square hits listed above, (2) gives respectively:

* Y=A: `(4+3s)Y=44+31s` has constant4 modulo8; every role fails.
* Y=8+5s: the product is62+44s; the equal-norm B8 role fails.
* Y=2+2s: the first congruence fails for j2 and6.
* Y=2A: the product is88+62s; exactly j6 survives.
* Y=3A: the product is132+93s; the A8 role fails.

The surviving j6,lambda1 case reconstructs exactly h. No equal-time
alternative is discarded.

For an isolated candidate write `z=aP-(Y/8)R0`. Since a is integral,
the spatial coordinate difference gives Y in4O. Integrality of the
other coordinates then gives Y in4sO, because 4+3s has valuation1
at s. Its norm is consequently `q=4u^2-alpha*a^2`. Modulo2 this
excludes q1 and qX. At q2 it forces a=s*a0; dividing by2 then
gives `1=2u^2-alpha*a0^2`, again impossible modulo2O. Thus all
three isolated roles are absent here too. It follows that h is
uniquely the actual first simple wall, and K has full actual D8
germ `(h,beta7,beta6,beta5,beta4,beta3,beta2,g)`.

## 3. An INTEGRAL involution matching the full actual germs

There is a rational involution S given by

    S(o)=P,
    S(e_i)=-alpha*o-(alpha/2)*sum(e_j)-e_(9-i).

It reverses the old positive A7 simples, sends g to d0, and sends
P to o. It is NOT integral. Defining X0=T*S*T, however, is integral:
let J=(2alpha;1+2s,alpha^7). A direct calculation gives

    J^2=-2alpha^-1,
    X0=J_J * w0(W),
    J_J(x)=-x+2J*(J,x)/J^2.

The first factor is point inversion about the future timelike J;
`2/J^2=-alpha` is integral. The second is the integral longest Weyl
element of the actual A7 on W. Both preserve L8, and their inverses
do too. They commute since J is perpendicular to W. Hence

    X0^2=I, trace(X0)=trace(S)=A-4alpha=1.

The complete actual-germ images are

    X0(g)=beta7, X0(beta1)=h,
    X0(beta2)=beta6, X0(beta3)=beta5,
    X0(beta4)=beta4, X0(beta5)=beta3,
    X0(beta6)=beta2, X0(beta7)=g, X0(P)=K.

The positive future P goes to the positive future K. Thus X0 maps
the FULL inward actual D8 germ at P to that at K; as an integral
arrangement isometry it stabilizes the same full L8 chamber.
This is not inferred from two matching abstract diagrams alone.

Finally extend X0 by identity on the exact summand O*r of M9.
The M9 divisor of r is (s). The complete norm-d crossing table in
the imported source has only nonzero products +/-1 or+/-k, all
units, so r in fact meets every other crossing mirror orthogonally.
The full intrinsic vertex germs above therefore give full ambient
germs D8+A1(r), and the extension stabilizes the ACTUAL M9 chamber.

This is one nontrivial FINITE-order symmetry. It supplies neither an
infinite-order product nor a compactness/nonreflectivity conclusion.
The different divisors of the earlier equal-norm walls r and c are
respected: X0 fixes r, and no r/c transposition was used. The next
genuine obligation is another compatible actual symmetry or a full
compactness argument. No worker, manuscript or Lean expansion is
part of this note.
