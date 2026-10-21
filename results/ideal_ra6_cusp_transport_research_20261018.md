# Cusp transport and a triangle bound for entirely ideal six-polytopes

18 October 2026. **Research proof awaiting independent review.** The proposed
new universal bound is

`V >= ceil(3645 F / 187)`, hence `V >= 4425` when `F >= 227`,

for an entirely ideal finite-volume right-angled hyperbolic six-polytope.
Here V counts its ideal vertices and F its facets. The existing internally
reviewed bound was V>=3859. This is a size bound, not an exclusion of
dimension six or an existence claim. No arithmetic or facet-spectrum
restriction is imposed.

## Rectangular cusp coordinates and triangle transport

Write m=n-1. At an ideal vertex choose future null vectors v,z with
v.z=-1 and an orthonormal positive basis e_1,...,e_m orthogonal to both.
After translating the horospherical rectangular box to its center, the
unit outward incident normals have the form

`n_i^+ = e_i+(w_i/2)v`, `n_i^- = -e_i+(w_i/2)v`, with `w_i>0`.

In particular opposite normals sum to w_i v. A cube vertex
epsilon in {+1,-1}^m specifies an ambient edge through v. Call this edge
good if all m two-faces containing it are ideal triangles. Call a cusp
cubic if all its w_i are equal; this is independent of the scale of v.

Put S=sum w_j^2. For the triangle obtained by dropping n_i^epsilon_i
from that edge's supporting walls, its third unit facet normal is

`c_i^epsilon = ((S-2w_i^2)/(4w_i))v +(2/w_i)z`
`               +sum_(j != i) epsilon_j(w_j/w_i)e_j`.                 (1)

Indeed it is orthogonal to the n_j^epsilon_j for j!=i, and has product
-1 with both n_i^+,n_i^-, since all three vertices of its polygon are
ideal. These linear equations and square one determine (1) uniquely.
All these are actual ambient facet normals: a facet cutting a side of
the right-angled two-face is orthogonal to its whole supporting white
forest.

For a good edge its m third walls are mutually orthogonal, and

`n_i^epsilon_i + c_i^epsilon = v_epsilon / w_i`,
`v_epsilon = (S/4)v+2z+sum_j epsilon_j w_j e_j`.                       (2)

The vector v_epsilon is future null, with v.v_epsilon=-2. It is
feasible against every other wall: (2) expresses it as a positive sum
of actual outward normals, each of which has nonpositive product with
any distinct outward normal. The selected walls also give nonpositive
products. It is orthogonal to all the edge-support normals. The null
directions v and v_epsilon are the two endpoints of that edge, and the
2m displayed normals give its full rectangular cusp link at the latter
endpoint. Its new widths, in the normalization (2), are 1/w_i. A good
edge therefore joins cubic cusps either at both ends or at neither end.

## Adjacent good edges constrain the widths

Distinct outward normals in a right-angled chamber have product zero
or at most -1. If their product were in (-1,0), the positive definite
two-root M-matrix would project a strict timelike center into an actual
nonright ridge, a contradiction. The same projection argument shows
that n mutually orthogonal actual facet normals force an actual finite
vertex; this is forbidden in an entirely ideal polytope.

Suppose two good cube vertices differ only in coordinate k. For each
i!=k there are three triangular edges of their i,k square. The third
walls c_i^epsilon and c_i^(epsilon with k reversed) have product

`1 - 2w_k^2/w_i^2`.

Both are orthogonal to c_k^epsilon. They and the n-3 unchanged white
support normals would force a finite vertex if that displayed product
were zero. The normal-product gap consequently makes it at most -1.
Thus w_k>=w_i for every i!=k: any direction joining adjacent good cube
vertices has maximal width.

If the cusp is not cubic, choose a nonmaximal direction. No pair of
good cube vertices lies on an edge in that direction. Its perfect
matching bounds the number g of good cube vertices by 2^(m-1).
If t cube edges correspond to triangular two-faces, each nongood
cube vertex touches a missing edge, so

`2(m*2^(m-1)-t) >= 2^m-g`.

In dimension six this becomes t<=64+g/2. At a noncubic cusp,
g<=16, giving t<=72.

## A moment bound for cubic cusp components

For a cubic cusp scale v so that every w_i=1. Its v is then uniquely
the sum of either pair of opposite unit normals. In dimension six put

`h=z-v/8`. Then `h.h=1/4` and `h.v=-1`.

Across a good edge the new positive basis is
`e_i'=n_i^epsilon_i-v_epsilon/2`. The new complementary future null
vector is `z'=h+v_epsilon/8`. It has product -1 with v_epsilon and
is orthogonal to every e_i'. Thus h=z'-v_epsilon/8: h is constant
throughout each finite good-edge component of cubic cusps.

Choose a future null vector Z perpendicular to this spacelike h,
independently for each component. Set E(v)=-Z.v>0 and
p_i=(Z.e_i)/E(v). The equations Z.h=0 and Z.Z=0 give
sum p_i^2=1/4. Every formal sign vector in (2) has ratio

`r_epsilon=E(v_epsilon)/E(v)=3/2-sum epsilon_i p_i`.

These ratios are positive even if the sign vector does not give an
actual good edge: its formal endpoint is future null and is not
parallel to Z, since its product with h is -1.
With mu=3/2 and d=-sum epsilon_i p_i, uniform sign averages give

`average d=average d^3=0`, `average d^2=1/4`,
`average d^4=3/16-2sum p_i^4 <= 13/80`.

Cauchy gives d^2<=5/4, so mu^2-d^2>=1. By sign symmetry and the exact
reciprocal identity,

`average 1/r = mu * average 1/(mu^2-d^2)`
` <= (3/2)[4/9+(16/81)(1/4+13/80)] = 71/90`.                        (3)

Here 1/(a-x)=1/a+x/a^2+x^2/[a^2(a-x)], with a=9/4 and x=d^2;
the denominator in the remainder is at least one. This is a universal
moment inequality, not a sampled sign or width assertion.

Use the positive function 1/E at each cubic cusp. Actual good neighbors
are a subset of its 32 formal neighbors, so their ratio sum is at most
32*71/90=1136/45. Summing over the finite components, every actual
undirected edge contributes a reciprocal pair with sum at least two.
Therefore the total good degree of all cubic cusps is at most
(1136/45)V_c, where V_c counts the cubic cusps.

## Consequences for face counts

Let T count triangular two-faces. At a six-dimensional ideal cusp the
cube has 80 edges, so sum t=3T. Combining the cubic good-degree bound
with t<=64+g/2 and the noncubic t<=72 gives

`3T <= (3448/45)V_c+72(V-V_c) <= (3448/45)V`.                         (4)

Write a_j for all j-faces and F=a_5. The entirely ideal
Dehn--Sommerville and cusp identities are

`a_1=16V`, `a_3=10V+2a_4-2F`, `a_2=25V+a_4-F`.

They also follow from Euler and summing the existing ideal-five
carrier identity over the facets. Every triangular polygon has three
vertices and every other polygon at least four. Since the total
polygon/cusp incidence is 80V,

`T >=4a_2-80V =20V+4(a_4-F)`.

Together with (4), this gives

`a_4-F <= (187/135)V`.                                                 (5)

Every ideal five-facet has at least56 four-facets by the earlier
universal carrier theorem, so 2a_4>=56F. Thus (5) implies
V>=3645F/187. The completed previous facet bound F>=227 then gives
V>=ceil(827415/187)=4425.

The existing 56 and 227 bounds are explicit internally proved inputs,
not external assumptions invented for this argument. Their proof chains
are documented in
[the earlier ideal-six packet](right_angled_ideal6_carrier_propagation_20261018.md).
The cusp coordinates, width restriction, actual-edge transport, invariant
and moment summation above require independent review before this new
size bound is promoted. The matching small script checks exact model and
coefficient identities; it does not encode the universal geometric proof.
There is no Lean or new dimension-bound claim.
