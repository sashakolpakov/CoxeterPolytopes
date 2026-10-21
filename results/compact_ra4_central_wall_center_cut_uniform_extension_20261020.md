# Uniform volume cost for central walls reaching the 120-cell center

20 October 2026. Written mathematical addendum, awaiting independent review.
No scientific jobs, searches, builds or formal verification. The fixed-wall
edition `compact_ra4_central_ball_completion_volume_cost_20261020.md` remains
unchanged at SHA256
`df5e75b3cfd8466b9b35f36b44f6c6b37f94b25afcdb276e7a443039b1d6da87`.

## 1. Parameter range and theorem

Retain its notation: P is compact right-angled in H4, D is an actual regular
right-angled dodecahedral facet with center o and unit outward normal n,
and all twenty outgoing edges have length L=arcosh(phi). Thus the actual
opposite vertex walls are m_v=sqrt(phi)v-phi n. Put

    h^2=phi/2,       B^2=1+h^2,       C=phi^2/sqrt(2),
    a=sqrt(phi)C=phi^(5/2)/sqrt(2).

Let an additional actual central wall be

    m=(o-pn)/r,        p=cosh(rho)>1,        r=sinh(rho).

In the collar formulation this is the ball of radius rho centered at o.
Compatibility with the base face planes forces rho below D's inradius.
Compatibility with the twenty vertex balls then forces strict disjointness:
the orthogonal alternative p=C/sqrt(phi) exceeds that inradius, while
the disjoint alternative is

    p<p_max=a-sqrt(3/2)=cosh(R_D-arcosh(sqrt(phi))).       (1)

Tangency is excluded by compactness. The regular 120-cell center matched
to D and this vertex star is

    q_0=sqrt(2)phi o-phi^(3/2)n.

Consequently m reaches or cuts q_0 exactly when

    p>=p_0=sqrt(2/phi).                                  (2)

**Theorem.** Throughout p_0<=p<p_max, the actual facet F_m has at least
fifty vertices. Every such compact completion has

    f(P)>=105,      W(P)>=23,      Vol(P)>=12pi^2>34pi^2/3.

Thus every central wall that reaches the matched 120-cell center is
excluded from every completion at or below the 120-cell volume. This
does not address the remaining interval 1<p<p_0 or off-center balls.

## 2. The projected cone works throughout the admissible range

The future unit point on m's wall is O=(po-n)/r. The projections of the
twenty vertex walls and twelve base face walls to m-perpendicular are

    m_v^U=A O+sqrt(phi)*sqrt(3)h w_v,      A=(pa-phi)/r,
    b_j^U=(ph/r)O+B r_j.

Their spatial data are the same regular dodecahedral and icosahedral
directions as in the fixed-wall proof. In particular, for every intrinsic
neighbor normal ell of m, the spanning sum sum_v m_v^U=20A O proves
t=-<ell,O>>0 and makes O an actual intrinsic interior point of F_m.

The face-average cone identity from that proof becomes

    b_j^U=gamma*average_face(m_v^U)+delta O,
    gamma=B^2/(sqrt(phi)h C),
    delta=[B^2 sqrt(phi)/C-p]/(hr)>0.                     (3)

The sign holds even at p_max, since

    B^2 sqrt(phi)/C-p_max=sqrt(3/2)-sqrt(2/phi)>0.

Thus every unknown neighbor ell still satisfies

    -<ell,b_j>>=1 for all twelve j.                      (4)

As before, this follows from the strict positive combination of a spanning
set, all-wall acuteness and the compact normal gap. It does not assume a
list of vertices or other walls of the completed central facet.

## 3. A uniform height bound

Set

    lambda=1/[sqrt(2phi)(a-p)],
    k=h r/(a-p)=phi r lambda,
    E^2=(11phi+7)/10.

The function A decreases with p: A'=(p phi-a)/r^3<0 by (1).
At the ends of the interval in the theorem,

    A(p_0)=phi^(3/2),      A(p_max)=sqrt(3/2)phi,
    lambda(p_0)=1/sqrt(5), lambda(p_max)=1/sqrt(3phi),
    k(p_0)=1/sqrt(5phi).

Both lambda and k increase with p; for the latter,

    k'=h(ap-1)/[r(a-p)^2]>0.

In particular 19/10<A<21/10, 3/2<E<8/5, and

    1/sqrt(5)<=lambda<1/2,       k>=1/sqrt(5phi).          (5)

Suppose some intrinsic neighbor has t=-<ell,O><=sqrt(3). The dodecahedral
support estimate from the fixed-wall proof gives

    max_v <ell,m_v> >= -At+E sqrt(1+t^2)>-27/40>-1.

For completeness, A>E makes the right side decreasing in t; use
t<7/4 and evaluate at 7/4 with the stated rational A,E bounds. Therefore
ell is orthogonal to at least one m_v, by the compact normal gap.

Use the exact local frame at that v. With y_j=-<ell,b_j>>=1 and
S=y_1+y_2+y_3>=3, the same expansion now gives

    t=kS,       y_j<=lambda S,
    1=sum_j y_j^2-phi lambda^2 S^2.                       (6)

These follow directly from <ell,m>=0 and the three adjacent vertex walls
m_w=sqrt(phi)q_v-phi b_j, with q_v=phi v-sqrt(phi)n. Assumption t<=sqrt(3)
and (5) imply

    3<=S<=S_0=sqrt(15phi)<5.

Since lambda>=1/sqrt(5)>2/5, we have 2lambda S>S-1. The elementary
three-variable convex maximum from the fixed-wall proof therefore yields

    1<=Q_lambda(S)
      =[1-2lambda+(2-phi)lambda^2]S^2-2(1-lambda)S+2.

For S>=3 and lambda<1/2 this expression strictly decreases with lambda:

    partial_lambda Q=2S[1-S+(2-phi)lambda S]<0.

It follows that 1<=Q_lambda(S)<=Q_(1/sqrt(5))(S), where

    Q_(1/sqrt(5))(S)=(9/5-phi)S^2-(12-4phi)S/5+2.

This is convex in S. Its values at both ends of [3,S_0] are below one:

    Q(3)=11-33phi/5<1,
    Q(S_0)=12phi-13-(4/5)(3-phi)sqrt(15phi)<1.

The second inequality is exact: both sides of

    (4/5)(3-phi)sqrt(15phi)>12phi-14

are positive, and the difference of their squares is 4(60phi-97)>0.
The last sign follows from sqrt(5)>67/30. This contradicts (6), proving

    -<O,ell>>sqrt(3) for EVERY intrinsic neighbor of m.    (7)

## 4. The complete facet fan and the volume conclusion

The normal-fan argument uses the actual compact facet F_m. At any one of
its vertices, let ell_i=t_i O+sqrt(1+t_i^2)u_i be its three intrinsic
unit normals. They are mutually Lorentz orthogonal, so

    u_i dot u_j=t_i t_j/sqrt((1+t_i^2)(1+t_j^2)).

Their spherical normal triangle has area Omega. With x_i=t_i^2>3,

    tan(Omega/2)=abs det(u_1,u_2,u_3)/(1+sum_{i<j}u_i dot u_j),
    det(u_1,u_2,u_3)^2=(1+sum_i x_i)/prod_i(1+x_i).

The determinant decreases in each x_i, while the denominator is at least
13/4. Hence tan(Omega/2)<sqrt(10)/26 and Omega<sqrt(10)/13. These
triangles partition the whole unit sphere, giving

    V(F_m)>52pi/sqrt(10)>156/sqrt(10)>49.

The final comparison follows by squaring: 156^2>10*49^2. Thus V(F_m)>=50.
The facet-capacity and neighboring-excess inequalities and volume identity
then give f>=105, W>=23 and Vol(P)>=12pi^2 as claimed.

## 5. A genuine limitation for smaller central balls

It is not valid to assert (7) uniformly as p decreases to one from these
local wall constraints alone. For sufficiently small r, adjoin the twelve
normals

    ell_j=h O+B r_j.

Their Gram matrix is exactly that of the regular dodecahedron, and all
are orthogonal to m. Their products with each n,b_i,m_v tend to minus
infinity as r tends to zero. Thus they, the central wall and the original
thirty-three walls form a finite acute right-angled prefix for sufficiently
small r. It has a common strict point just inside D. All old walls restrict
to redundant strict inequalities on m's plane, so the actual central facet
of this prefix is exactly a regular dodecahedron with wall heights h<sqrt(3).

No compactness or finite volume of this prefix is asserted. This is a
counterexample only to a uniform local height inference. The smaller-radius
completion question needs an additional global argument, and the theorem
above does not silently assume it.

Dependencies are precisely the unchanged fixed-wall proof and its two
pinned inherited compactRA4 inputs: the actual-facet collar convention,
compact normal gap, capacity, neighboring excess and volume formula.
No universal minimum, decomposition or whole-120-cell recognition is used.
