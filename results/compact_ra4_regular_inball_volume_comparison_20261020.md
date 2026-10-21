# The regular inball implies the 120-cell volume bound

20 October 2026. Written conditional theorem, awaiting independent review.
No executions, searches, builds or formal verification. The comparison
of spherical simplices below is proved directly, rather than assumed as
a spherical isodiametric theorem.

Put phi=(1+sqrt(5))/2 and

    T0=phi^(3/2),       R4=asinh(T0).

**Theorem.** If a compact right-angled hyperbolic four-polytope P contains
a ball of radius R4, then

    Vol(P)>=34pi^2/3.

Equality holds exactly for the regular right-angled 120-cell, which does
contain that ball. There is no hypothesis about a marked facet, its
outgoing edge lengths, a common center for several dodecahedral facets,
or the shapes of the other facets.

The unproved global step is that an arbitrary P must contain such a
ball. The farthest-wall inequality max_n(-<q,n>)>=T0 does not give the
required minimum over all walls. No universal minimum-volume theorem
is deduced here.

## 1. A direct spherical-simplex volume comparison

Let u_1,...,u_4 be linearly independent unit vectors in Euclidean R4,
and let G be their Gram matrix. Suppose

    G_ij>=c for i!=j,       0<c<1.

Their positive cone cuts out a spherical tetrahedron in S3. Write its
three-volume as Omega(G). Let G_c=(1-c)I+c*11^t be the Gram matrix of
the regular spherical tetrahedron with pair products c. Then

    Omega(G)<=Omega(G_c).                               (1)

Here are the two elementary estimates and an exact integral proving it.
The trace of G is four, and its greatest eigenvalue lambda satisfies

    lambda>=1+3c

by the Rayleigh quotient of the vector (1,1,1,1). Arithmetic-geometric
mean on its other three positive eigenvalues gives

    det G<=lambda*((4-lambda)/3)^3
          <=(1+3c)(1-c)^3=det G_c.                      (2)

The last scalar function decreases for lambda>=1, since its derivative
is 4(1-lambda)(4-lambda)^2/27. On the simplex

    Delta={a_i>=0:sum_i a_i=1},

the other estimate is pointwise:

    a^t G a>=c+(1-c)sum_i a_i^2=a^t G_c a.              (3)

Use coordinates (a_1,a_2,a_3) on Delta, with
a_4=1-a_1-a_2-a_3. The exact cone integral is

    Omega(G)=sqrt(det G)*integral_Delta(a^t G a)^(-2) da. (4)

To verify its normalization, put U=(u_1,...,u_4). On the cone U R_+^4,
the Gaussian integral of exp(-|y|^2/2) is 2Omega(G), since
integral_0^infinity exp(-r^2/2)r^3 dr=2. Change variables y=Ux and
then x=s*a with a in Delta. The latter Jacobian is s^3, and

    integral_0^infinity exp(-s^2(a^t G a)/2)s^3 ds
       =2(a^t G a)^(-2).

Dividing by two proves (4). Thus there is no unspecified shape-dependent
factor. Combining (2), (3) and (4) proves (1). Equality is possible only
for G=G_c: if an off-diagonal entry is larger than c, (3) is strict on
the interior of Delta, while (2) never increases its prefactor.

## 2. Apply the comparison to every actual vertex

Choose the center q of a contained ball of radius R4. For every outward
unit actual wall normal n, its positive height satisfies

    t_n=-<q,n>>=T0.

Conversely these inequalities are exactly the containment condition for
the ball, because P is the intersection of its inward halfspaces.

A compact right-angled four-polytope is simple. At an actual vertex its
four unit incident normals n_i are mutually Lorentz orthogonal. In the
Euclidean four-space q-perpendicular decompose

    n_i=t_i q+sqrt(1+t_i^2)u_i.

Then u_i are unit vectors, and for i!=j

    u_i dot u_j=t_i t_j/sqrt((1+t_i^2)(1+t_j^2))
                >=T0^2/(1+T0^2)=phi/2=:c.

Their Gram matrix is positive definite: before normalizing the spatial
vectors it is I+t t^t. These u_i are precisely the Euclidean facet-normal
directions of the bounded Klein image of P centered at q. Hence their
positive cones, one at each actual vertex, form its complete normal fan.
They cover S3 with disjoint interiors and total volume 2pi^2.

By (1), each normal tetrahedron has volume at most the regular value for
c=phi/2=cos(pi/5).

## 3. The regular comparison value and the volume conclusion

The explicit classical comparison model is the regular compact
right-angled 120-cell T with its 120 regular dodecahedral facets. Its
600 vertices follow also from 120*20=4*600. At its center all wall
heights equal T0. Indeed projection onto a facet sends the center to
the center of its regular dodecahedron, whose intrinsic wall-height
square is h^2=phi/2. Therefore

    h^2=T0^2/(1+T0^2),       T0^2=phi^3.

Its normal fan consists of 600 congruent regular spherical tetrahedra
of pair product phi/2: the ordinary regular 600-cell triangulation of
S3. This classical regular model and its facet count are the explicit
external metric/combinatorial input for the comparison; their use does
not presume that an unknown P is assembled from such cells. In particular

    Omega(G_c)=2pi^2/600.

The complete actual fan of P consequently gives

    2pi^2=sum_vertices Omega(G_v)
           <=v(P)*(2pi^2/600),
    v(P)>=600.

The inherited compactRA4 identities are

    W=v-5f>=0,
    Vol(P)=pi^2(v-4f+16)/12.

Since v-4f=(v+4W)/5>=v/5, they imply

    Vol(P)>=pi^2(600/5+16)/12=34pi^2/3.

Equality forces v=600 and W=0. The reviewed zero-excess recognition
then identifies the entire P with T. Conversely T attains equality.
The theorem makes no assertion of uniqueness for a prescribed ball
center in an arbitrary P; that is unnecessary for volume recognition.

## Scope and pinned source interface

`manuscripts/compact_ra4_volume_constraints.tex`, SHA256
`afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`,
supplies the actual-wall and simplicity conventions, the compact volume
and excess identities, the regular-120-cell comparison data, and the
zero-excess geometric recognition. Its three-dimensional metric input
identifies a compact right-angled all-pentagonal facet with the regular
dodecahedron. The spherical comparison (1)-(4) is proved above, not
imported from that source.

The hypothesis is an inscribed ball of radius R4. It is stronger than
the existing maximum-wall-radius lower bound. Showing that every compact
right-angled P has this inradius, or dealing with smaller inradius by a
different volume argument, remains open. The result also does not prove
a second-volume theorem or a global garland decomposition.
