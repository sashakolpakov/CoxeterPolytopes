# Two short directions at an all-right finite corner

20 October 2026. Own paper argument; no computation or external
classification is used. This gives a local capacity bound that excludes
one vertex type from the global drift witness in dimensions at least14.
It does not exclude all compact Coxeter polytopes in those dimensions.

## Statement

Let p be a finite vertex of a compact hyperbolic Coxeter n-polytope, with
all incident facets mutually orthogonal. Set

    L_* = asinh(1/sqrt2),    tau_* = tanh(L_*/2)=sqrt3-sqrt2.

At most TWO edges incident to p have length strictly less than L_*.
The count is of edge directions, including the possibility that one
opposite facet closes two directions.

## The actual-wall coordinates

Use future p^2=-1 and outward orthonormal normals n_1,...,n_n. An actual
facet c not through p has the unique expression

    c=a p-sum_j beta_j n_j,   a=-p.c>0,   beta_j=-c.n_j>=0,
    a^2=sum_j beta_j^2-1.

Every nonzero beta_j is cos(pi/m), m>=3, or is strictly greater than1.
This is the Coxeter Gram alternative for actual walls. The value1 is
excluded by compactness: n_j+c would be a nonzero null vector in every
wall halfspace, and hence an ideal point of the polytope. The vector is
future because its product with an interior future point is negative.

Suppose c closes direction i at the adjacent actual finite vertex. The
common n_j, j!=i, together with c have positive definite Gram, so

    delta_i=1-sum_(j!=i) beta_j^2>0,
    cosh^2 L_i=beta_i^2/delta_i,
    a=sqrt(delta_i) sinh L_i.

Consequently L_i<L_* implies

    a^2<1/2,    beta_i^2<(3/2)delta_i.                    (1)

## Complete support census for a short closing wall

All nonzero beta_j^2 are at least1/4. From (1), the support has at most
three elements. Indeed four support elements force the three other
squares to be exactly1/4 (otherwise delta_i<=0); then beta_i^2<3/8
forces beta_i^2=1/4 and a=0, a contradiction.

With three support elements, the released square is below3/4. The only
possible native values are1/4,1/2, and (3+sqrt5)/8. In the first two
cases the other two squares would have to lie respectively in
(3/4,5/6) or (1/2,2/3), intervals containing no sum of two nonzero native
squares. In the last case their sum must be1/2. Thus this profile is
(3,3,5), and

    a^2=h=(sqrt5-1)/8.

For two support elements, if one contact is dotted its square exceeds1,
and shortness forces the other contact to be3 and the dotted square to
be below9/8. Its height satisfies

    1/4<a^2<3/8.                                        (D)

If both contacts are finite, consider the larger one; its edge is the
shorter of the two. The smaller contact must be3 or4. The resulting
profiles are precisely

    (3,m), m>=7, with 0<a^2<1/4;                        (F)
    (4,5), with a^2=h.

For (3,m), the heavy-contact edge is short, and the light-contact edge
can also be short only if a^2<1/12. For (4,5) both directions are short.
An isolated short wall is dotted, with 1<beta_i^2<3/2.

This census only uses the ordered native values
1/4,1/2,(3+sqrt5)/8,3/4,... and the actual finite endpoint condition.
It makes no assumption that every formal profile is realizable.

## Pair compatibility

Let c,d be two DISTINCT actual walls closing short directions, with
heights a,b. Their product is

    c.d=sum_j beta_j gamma_j-ab >-1/2.

All-wall acuteness and the Coxeter gap therefore force c.d=0. Their
supports meet in exactly one coordinate: no common coordinate would
give negative product, while two common coordinates contribute at least
1/2>ab. In their common coordinate,

    beta_j gamma_j=ab.                                 (2)

An isolated dotted short wall cannot coexist with another: its common
contact product exceeds1/2. A (4,5) wall cannot coexist with another,
because its contact product is at least1/(2sqrt2), while ab<sqrt(h/2).

For a (3,3,5) wall, the shared contact must be3 on BOTH walls; any larger
contact already makes its product too large. Every short wall containing
a3 contact has squared height below3/8. But

    h*(3/8)<1/16,

so (2), whose square is1/16, is impossible. Thus this wall also occurs
alone among short closing walls.

Only profiles (F) and (D) remain. Their large contacts have size greater
than sqrt3/2, while every height is below sqrt(3/8). Therefore the shared
coordinate must be the3 contact of BOTH walls. Equation (2) is then
ab=1/4. Two (F) walls have ab<1/4, and two (D) walls have ab>1/4.
Hence there are at most two distinct short walls, one of each type.

A (D) wall closes only its dotted direction at a finite endpoint. If an
(F) wall closes two short directions, its height satisfies a^2<1/12;
compatibility with a (D) wall would instead require
a^2=(1/16)/b^2>1/6. Thus the two-wall case closes at most one short
direction per wall. The one-wall cases close at most two. This proves
the stated capacity bound.

## Consequence for any common-white drift

Use the identity proved in
`results/compact_finite_edge_white_weighted_drift_20261020.md`.
At an A1^n vertex all n common-white diagrams are the SAME A1^(n-1),
so every globally shared positive white weight has a common factor f.
The bound above gives

    A_f(p)>=f(n-2)(sqrt3-sqrt2),
    ||rho_f(p)||=f sqrt n.

For n>=14 the first expression is strictly larger. At n14 this follows
from 353^2>6*144^2 (difference193); the ratio (n-2)^2/n is increasing
for n>2. Thus the vertex supplied by the global drift inequality
A_f<||rho_f|| CANNOT have diagram A1^n in dimensions at least14,
for ANY chosen globally shared positive f.

Other finite vertex diagrams, and compatibility between their white
weights, remain necessary for a dimension exclusion. No new uniform
dimension ceiling is claimed here.
