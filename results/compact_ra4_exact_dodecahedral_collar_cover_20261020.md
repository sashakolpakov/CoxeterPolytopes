# Exact collar and orthogonal-ball-cover reduction for the compact volume problem

20 October 2026. Written mathematical reduction, awaiting independent
review. No executions, numerical searches, builds or Lean runs.

Let P be compact and right-angled in H4, and fix an actual facet F with
outward unit normal n. The argument applies to any compact intrinsic F;
the new farthest-facet lemma permits choosing F to be the fixed regular
dodecahedron D. It gives an exact fixed-domain formulation of the
remaining universal volume problem, not its solution.

## 1. The whole polytope is a variable normal collar

Let b_j be the neighboring actual wall normals of F. They are
orthogonal to n and are exactly the intrinsic face normals of F.
Every point of the inward halfspace of n has the unique expression

    y(t,x)=cosh(t)x-sinh(t)n,     x in H3=n-perpendicular, t>=0.

For such a neighboring wall,

    <y(t,x),b_j>=cosh(t)<x,b_j>.

Thus every point of P projects into the ACTUAL F, not merely its
supporting plane. Conversely the b_j impose precisely x in F.

Every remaining wall m is disjoint from F and has <m,n><-1 by the
compact acute-wall gap. Write uniquely

    m=sinh(L_m)v_m-cosh(L_m)n,       L_m>0,

with v_m a future unit timelike vector in n-perpendicular. To justify
the future choice, its timelike projection has negative product with
any relative interior point of F, because m is strict there. Acuteness
against each b_j gives <v_m,b_j><=0. Hence v_m belongs to F.

The inequality of m becomes

    tanh(t) <= g_m(x),
    g_m(x)=tanh(L_m) cosh dist_F(x,v_m).

Put g(x)=min_m g_m(x). Compactness implies g(x)<1 at every x in F:
otherwise the whole normal ray y(t,x), t>=0, would satisfy all wall
inequalities. Since F is compact and there are finitely many walls,
g is continuous and uniformly below one. It is strictly positive.
Consequently the entire P is exactly

    { y(t,x): x in F, 0<=t<=h(x) },     h(x)=artanh g(x).

Its volume is therefore exactly

    Vol(P) = integral_F [sinh h(x)+(sinh^3 h(x))/3] dVol_F(x),

because the normal metric is dt^2+cosh^2(t) ds_F^2. This is a full
polytope description, not a check at finitely many vertices.

## 2. Equivalent orthogonal balls inside the fixed facet

Define rho_m>0 by

    cosh rho_m = coth L_m,
    sinh rho_m = 1/sinh L_m.

Then g_m(x)<1 precisely when x belongs to the open intrinsic hyperbolic
ball B(v_m,rho_m). The compactness condition just proved says that
these finitely many open balls cover ALL of F, including its boundary.

They satisfy the following exact compatibility conditions.

* For every face plane of F, the center v_m lies on that plane, or
  its inward distance from it is strictly greater than rho_m. Indeed
  <m,b_j>=-sinh(L_m)sinh dist(v_m,wall(b_j)) is either zero or <-1.
  Thus a ball meets a boundary plane orthogonally through its center,
  or stays strictly away from that plane.
* For distinct balls with centers v_i,v_j, their associated walls
  give either

      cosh d(v_i,v_j)=cosh rho_i cosh rho_j,

  or

      d(v_i,v_j)>rho_i+rho_j.

  The first case is exactly orthogonality of the two sphere boundaries;
  the second is strict disjointness. These formulas follow directly
  from the scalar products of the displayed wall normals. Tangency is
  excluded by compactness, not silently treated as orthogonality.

For example, at an intersection point x of the two spheres, their
intrinsic normal vectors are v_i-cosh(rho_i)x and
v_j-cosh(rho_j)x. Their product is
cosh(rho_i)cosh(rho_j)-cosh d(v_i,v_j), verifying the first assertion.

## 3. Converse: a precise construction/minimization model

Conversely take finitely many centers v_i in a compact right-angled
H3 polyhedron F and radii rho_i>0 satisfying both compatibility
conditions above, and suppose their OPEN balls cover F. Embed its
normal hyperplane in H4 and adjoin the unit normals

    m_i=(v_i-cosh(rho_i)n)/sinh(rho_i)

to n and the intrinsic normals b_j of F. Their pair products are
zero or strictly below -1, by exactly the preceding identities.
The intersection of the inward halfspaces is the exact normal collar
with

    g(x)=min_i cosh d(x,v_i)/cosh rho_i <1.

It is compact, convex and full-dimensional. All its actual ridge
angles are right. Each candidate m_i really is a facet: at its own
center v_i its profile is 1/cosh rho_i<1, while any other profile is
at least cosh rho_i>1, from the pair compatibility condition. Hence
m_i uniquely cuts the normal ray there and remains the unique roof
wall over a nearby open subset of the interior of F. The base and
its neighboring b_j walls are also actual, by the unchanged base F
and sufficiently small positive normal heights. This prevents a
construction claim based on redundant or purely formal mirrors.

Thus compact right-angled H4 polytopes with a marked dodecahedral
facet are equivalent to these finite compatible strict ball covers of
the SAME fixed D, through an exact geometric correspondence.

## 4. The remaining cost inequality and a verified limitation

For F=D, a proof that every such compatible cover has the displayed
volume integral at least 34pi^2/3 would prove the universal 120-cell
minimum. Equality would still require a geometric recognition argument;
one must not insert the desired cover or polytope in advance.

The current pointwise wall-radius theorem only controls the largest
wall distance from each interior point. It does not imply a lower
bound on every collar height. The outgoing edge lengths s_v at the
twenty vertices v of D satisfy

    coth(s_v)coth(s_w)<=phi

when v,w are adjacent. These conditions alone do not identify all
centers/radii or the roof g. In particular a shorter outgoing edge
can be accompanied by longer neighboring edges. A proof based solely
on replacing every s_v by the regular pentagon length is invalid.

Likewise concavity of g in homogeneous coordinates supplies only a
weak lower envelope from vertex values; it does not recover the full
120-cell roof in the interior. The cover conditions across the whole
D, including all pair and boundary compatibilities, are the missing
global mechanism. No new universal volume bound or second-volume
claim is made in this note.

Inputs, preserved byte-identically:

* `manuscripts/compact_ra4_volume_constraints.tex`,
  SHA256 afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f:
  compact acute-wall gap and actual facet inheritance.
* `results/compact_ra4_farthest_dodecahedron_wall_radius_20261020.md`,
  SHA256 d082e4a918d4504a59ac200f32bbb141f1f2cb52d701006eaea7ba5b0e5bcc22:
  the universal choice of a marked D facet and the scope of the
  wall-radius bound.
