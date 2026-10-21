# A central collar wall forces a volume cost

20 October 2026. Written proof, awaiting independent review. No scientific
execution, build or formal verification. This is a conditional completion
theorem for one specified local configuration, not a universal minimum-volume
theorem or an obstruction to every compact completion.

## Statement and conventions

Use the Lorentz form with a negative time coordinate. Put

    phi=(1+sqrt(5))/2,       L=arcosh(phi).

Let P be a compact right-angled hyperbolic four-polytope with an actual
regular right-angled dodecahedral facet D. Let n be its outward unit normal
and o its intrinsic future unit center. Suppose that the twenty edges of P
leaving the vertices of D all have length L. At each vertex v of D, the
opposite actual wall then has outward unit normal

    m_v=sqrt(phi) v-phi n.

Suppose P also has the actual wall

    m=(8o-9n)/sqrt(17).

In the exact collar/ball correspondence, this last wall is the ball centered
at o with radius arcosh(9/8). The theorem concerns completions retaining
these specified walls; no sufficiency of this finite prefix is assumed.

**Theorem.** The facet of m has at least fifty vertices. Consequently

    f(P)>=105,       W(P)>=23,       Vol(P)>=12 pi^2,

which is strictly greater than the regular 120-cell volume 34 pi^2/3.
The conclusion allows larger-volume compact completions and does not
identify all their other walls.

The only inherited four-dimensional inequalities used below are

    f(P)>=2V(F)+5,
    W(P)>=ceil(3(V(F)-20)/4),
    Vol(P)=pi^2(f(P)+W(P)+16)/12.

In particular the separate positive-excess gap W=0 or W>=17 is not needed.

## 1. Fixed dodecahedral data

Write the twelve intrinsic face normals and twenty vertices of D as

    b_j=h o+B r_j,          h^2=phi/2,       B^2=1+h^2,
    v=C o+sqrt(3)h w_v,     C^2=1+3h^2.

Here the r_j and w_v are the unit icosahedral and dodecahedral directions
in the Euclidean three-space perpendicular to o and n. These constants
follow by placing the regular polyhedron at its center: three mutually
orthogonal b_j pass through each v, and o=Cv-h(b_1+b_2+b_3). Adjacent
vertices satisfy w=phi v-sqrt(phi)b_j for the face wall released along
their edge. This also proves the displayed formula for the outgoing wall
m_v and its finite endpoint phi v-sqrt(phi)n.

Introduce

    O=(9o-8n)/sqrt(17),          X=8sqrt(phi) C=sqrt(160phi+96),
    A=(9sqrt(phi) C-8phi)/sqrt(17),
    E^2=(11phi+7)/10,
    k=h sqrt(17)/(X-9),         lambda=4sqrt(2/phi)/(X-9).

O is a future unit point on m's supporting hyperplane. The projections
of m_v and b_j to that hyperplane are

    m_v^U=A O+sqrt(phi)*sqrt(3)h w_v,
    b_j^U=(9h/sqrt(17))O+B r_j.

The central wall is strictly ultraparallel to all thirty-three prescribed
walls n,b_j,m_v. For n and b_j this follows from

    <m,n>=-9/sqrt(17),        <m,b_j>=-8h/sqrt(17).

For the remaining walls,

    <m,m_v>=-(X-9phi)/sqrt(17)<-1.

Indeed X^2-(9phi+sqrt(17))^2=79phi-2-18phi sqrt(17)>4phi-2>0,
using sqrt(17)<25/6. Thus none of these thirty-three walls is adjacent
to m. Every neighbor of m is an additional wall of a putative completion.

For use in the estimates, the following intentionally loose exact bounds
will suffice:

    19/10<A<21/10,       3/2<E<8/5,
    k>3/8,              9/20<lambda<23/50.                 (1)

Here is rational arithmetic verifying them. Use 161/100<phi<81/50 and
94/5<X<377/20. They give the A and E bounds directly, and

    k^2>5474/38809>9/64.

Since lambda^2=32/[phi(X-9)^2], its lower and upper bounds are equivalent
respectively to

    1458X>27137+160phi,
    9522X<173633+4640phi.

The first follows from X>94/5 and phi<81/50; the second follows already
from X<19 and phi>8/5. All quantities squared here are positive.

## 2. All neighbors of the central wall avoid the base face planes

Let ell be any actual outward unit wall normal adjacent to m. Thus ell
lies in U=m-perpendicular, and acuteness gives <ell,m_v><=0 for all v.
The twenty projected m_v^U span U and sum to 20A O. Hence

    t=-<ell,O>>0:                                            (2)

if the sum of these nonpositive products were zero, every product would
be zero, contradicting their spanning property. This proves that O lies
strictly inside every intrinsic wall of the actual compact facet F_m.
Thus O is an intrinsic interior point of F_m. Nonadjacent parent walls
restrict to redundant strict inequalities on the whole hyperbolic plane
of a compact right-angled facet, by the compact acute-wall gap.

More is true. Average the five m_v^U whose vertices lie on the j-th face
of D. Since their average spatial direction is C r_j/(sqrt(3)B),

    b_j^U=gamma*average_face(m_v^U)+delta O,
    gamma=B^2/(sqrt(phi)h C)>0,
    delta=[-9+8B^2 sqrt(phi)/C]/(h sqrt(17))>0.

For the last sign, squaring the positive terms gives

    64 B^4 phi-81 C^2=(77/2)phi-1>0.

As O is a positive multiple of the sum of all twenty m_v^U, this is a
strictly positive combination of a spanning collection. Therefore

    <ell,b_j><0 for all twelve j.

These are distinct actual walls of the compact right-angled P. Its
normal gap now gives the important uniform inequalities

    -<ell,b_j>>=1.                                         (3)

This step uses all the prescribed vertex walls, not an assumed list of
the unknown boundary of F_m.

## 3. Every intrinsic wall height at O exceeds 7/4

Write ell=tO+sqrt(1+t^2)u, where u is a unit Euclidean direction.
The Euclidean convex hull of the twenty w_v has face planes
r_j dot x=C/(sqrt(3)B); it contains the ball of this radius. Thus some v
satisfies u dot w_v>=C/(sqrt(3)B), giving

    max_v <ell,m_v> >= -At+E sqrt(1+t^2).                  (4)

The expression on the right decreases with t, because A>E. If t<=7/4,
(1) and sqrt(65)>8 imply that it is greater than

    [E sqrt(65)-7A]/4 > -27/40 > -1.

Every distinct actual pair of compact right-angled walls has product
zero or strictly less than -1. Hence ell must be orthogonal to at least
one m_v.

Fix this v and let b_1,b_2,b_3 be the three face walls of D through it.
Set q_v=phi v-sqrt(phi)n. The orthogonal frame q_v,b_1,b_2,b_3 lies in
m_v-perpendicular, with q_v^2=-1 and b_j^2=1. By (3), write

    ell=a q_v-y_1 b_1-y_2 b_2-y_3 b_3,
    y_j>=1,       S=y_1+y_2+y_3>=3.

Expanding m in this frame gives

    m=-mu m_v+beta q_v-(8h/sqrt(17))(b_1+b_2+b_3),
    beta=(8phi C-9sqrt(phi))/sqrt(17)>0.

The equation <ell,m>=0 and then the product with O give, exactly,

    a=cS,       c=4sqrt(2)/(X-9),       t=kS,
    1=y_1^2+y_2^2+y_3^2-c^2 S^2.                         (5)

The three neighbors w of v in D satisfy
m_w=sqrt(phi)q_v-phi b_j. Acuteness against these actual walls implies

    y_j<=lambda S,       c^2=phi lambda^2.                (6)

Still assuming t<=7/4, (1) gives S<14/3. Put M=lambda S. Then
2M>S-1, since lambda>9/20 and S<14/3. A convex maximum over
1<=y_j<=M and sum y_j=S yields

    y_1^2+y_2^2+y_3^2<=M^2+(S-M-1)^2+1.

For clarity, its possible maximizing extreme points are permutations of
(M,S-M-1,1), or (S-2,1,1). The alternative with two coordinates M is
excluded by 2M>S-1. If the second type occurs, its value is no greater
than the displayed bound, since x^2+(S-1-x)^2+1 increases for
x>=(S-1)/2 and M>=S-2>=(S-1)/2.

Substituting in (5)-(6), and using phi>8/5 and (1), would force

    1 <= [1-2lambda+(2-phi)lambda^2]S^2-2(1-lambda)S+2
       < (37/200)S^2-(27/25)S+2.

The last convex quadratic is at most the larger of its endpoint values
on [3,14/3], namely 17/40 and 89/90. Both are below one. This
contradiction proves

    -<O,ell>>7/4 for EVERY intrinsic wall of F_m.          (7)

## 4. The normal fan forces at least fifty vertices

At a vertex of F_m, its three actual intrinsic wall normals are mutually
orthogonal. Write them as ell_i=t_i O+sqrt(1+t_i^2)u_i. Their Euclidean
unit directions satisfy

    u_i dot u_j=t_i t_j/sqrt((1+t_i^2)(1+t_j^2))>0.

Use the Klein chart centered at O. The Euclidean outward facet normals
are precisely these u_i. Its complete normal fan partitions the unit
two-sphere into one spherical triangle per vertex of F_m. For the
triangle of the displayed three directions, let its area be Omega.
The elementary spherical triangle half-area identity gives

    tan(Omega/2)=abs det(u_1,u_2,u_3)/(1+sum_{i<j}u_i dot u_j).

There is no branch ambiguity: all three pair products are positive,
so this is the small spherical triangle, of area less than pi.
Set x_i=t_i^2. The determinant lemma for a diagonal matrix plus rank one
gives

    det(u_1,u_2,u_3)^2=(1+x_1+x_2+x_3)/prod_i(1+x_i).

This decreases in each x_i. By (7), x_i>49/16; also the denominator in
the half-area identity is at least 1+3(49/65)=212/65. Consequently

    tan(Omega/2) <= 4sqrt(163)/(53sqrt(65)) <32/265,
    Omega<64/265.

The last rational comparison is simply 163/65<64/25. Summing over the
complete normal fan, not merely over the known local vertices, yields

    4pi < V(F_m)*64/265,
    V(F_m)>265pi/16>795/16.

Thus V(F_m)>=50. The inherited capacity and neighboring-excess theorems
now give f(P)>=105 and W(P)>=ceil(90/4)=23. The volume identity gives
Vol(P)>=12pi^2 as asserted.

## Dependencies and scope

The fixed regular-dodecahedron metric is used explicitly, together with
actual-wall acuteness and the compact zero-or-below-minus-one normal gap.
The fan is the ordinary Euclidean normal fan of the bounded Klein image
of the actual facet; its vertices are not assumed to be among any finite
prefix of known vertices. The spherical half-area identity is elementary
spherical trigonometry. No four-dimensional recognition or minimum-volume
theorem, complete collar classification, or unproved continuation is used.

Pinned inherited inputs, unchanged:

* `manuscripts/compact_ra4_volume_constraints.tex`, SHA256
  `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`:
  compact normal gap, facet capacity, neighboring excess, volume identity.
* `results/compact_ra4_exact_dodecahedral_collar_cover_20261020.md`, SHA256
  `c1a3f94c880d37a6180a7b83840eee05c27883eeb8161d77f63ddaa26b942669`:
  the interpretation as a complete open-ball cover and its actual-facet
  convention. Its existence claim is not used to assume a completion here.

The choice of a regular D is a hypothesis of this conditional theorem.
The separate universal farthest-dodecahedron result supplies such facets,
but does not supply the twenty equal outgoing lengths or the central wall.
Accordingly this closes a concrete nonregular collar branch below the
120-cell volume; it does not yet close the unrestricted volume problem.
