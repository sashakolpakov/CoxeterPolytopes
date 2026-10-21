# Independent hostile audit: finite-area K3 Coxeter triangles

**Verdict (2026-08-11): PASS under the exact full-trace and cyclic-ground-field
hypotheses.**  The theorem includes ideal vertices and proves a finite field
set, not merely degree at most six.  Its cyclic integrality argument remains
valid for nonprimitive choices of facet normals.

## 1. Rank and all-place incidence bounds

The order module has rank three, while its primitive full trace is the
entire Picard lattice of a complex projective K3.  Therefore

\[
                         3d=\rho(X)\le20,
                         \qquad d\le6.
\]

At the distinguished place, two adjacent walls of a finite-area Coxeter
triangle meet at a finite Coxeter angle or at an ideal vertex.  Their pair
cycle (p_{ij}) consequently lies in ([0,4]).  At every definite place,
strict Cauchy--Schwarz gives (0\le\sigma(p_{ij})<4).  Thus the weak bound

\[
                         |\sigma(p_{ij})|\le4
\]

is correct at every embedding; equality at an ideal distinguished vertex
is deliberately retained.

For the directed triple cycle (c), symmetry gives

\[
                         c^2=p_{01}p_{12}p_{20}.
\]

All pair coordinates are nonnegative at every real place, so
(|\sigma(c)|\le8).  No compactness-only strict inequality is used.

## 2. Integrality without primitive normals

Put (T_i=1-s_i\in\operatorname{End}_{\mathcal O}(M)).  Over the field,
(T_i) is rank one:

\[
                  T_i(x)=\alpha_i(x)e_i,
                  \qquad \alpha_i(e_i)=2.
\]

Although a rescaling of (e_i) can make an individual directed Cartan
entry (a_{ij}=\alpha_i(e_j)) nonintegral, the closed products are traces
of integral endomorphisms:

\[
 \operatorname{Tr}(T_iT_j)=a_{ij}a_{ji}=p_{ij},\qquad
 \operatorname{Tr}(T_0T_1T_2)=a_{01}a_{12}a_{20}=c.
\]

The trace of an endomorphism of a finite projective order-module lies in
(mathcal O).  Hence (p_{ij},c\in\mathcal O) without assuming that the
normal is primitive or that either directed entry is integral.  Expanding
(s_i=1-T_i) recovers exactly

\[
 \operatorname{Tr}(s_is_j)=-1+p_{ij},
\]

and

\[
 \operatorname{Tr}(s_0s_1s_2)
 =-3+p_{01}+p_{12}+p_{20}-c.
\]

The executable audit independently checks these identities and invariance
of every cycle under arbitrary rescaling of the three normals.

## 3. Closed walks and the actual ground field

Every closed directed walk factors into simple directed cycles: cut at two
occurrences of a repeated vertex and recurse.  On three vertices, the
simple cycles have only these forms:

- length one, giving (a_{ii}=2);
- length two, giving one of the (p_{ij}); or
- length three, giving (c) or its reverse, equal to (c) by symmetry.

Therefore all cyclic Cartan products lie in
(mathbf Z[p_{01},p_{12},p_{20},c]).  This also handles walks such as
(0,1,0,2,0), which factor as (p_{01}p_{02}) but defeat the weaker
"global alternation or triangle" shortcut.  The frozen carrier hypothesis
states literally that the complete cyclic packet generates (k).

## 4. Primitive grid and discriminant

A minimal generating subpacket has at most two elements: three strict field
adjunctions would force degree at least eight.  For joint generators
(\theta_1,\theta_2), every unordered pair of the at most six embeddings
forbids at most one coefficient (m) in
(\theta_1+m\theta_2).  There are at most ({6\choose2}=15) such pairs,
so one of (m=0,\ldots,15) is primitive.

Every packet coordinate has conjugate magnitude at most eight.  The chosen
primitive algebraic integer therefore has conjugate magnitude at most
(16\cdot8=128), and pairwise differences are at most (256).  The
squared Vandermonde gives

\[
 |\operatorname{Disc}(k)|
 \le |\operatorname{Disc}(\mathbf Z[\theta])|
 \le256^{d(d-1)}\le256^{30}.
\]

Hermite--Minkowski supplies a literal finite set of fields.

The theorem neither constructs the carrier nor proves essential
surjectivity.  It covers compact and ideal-vertex finite-area triangles;
truncated or higher-facet cases remain separate.

The independent executable audit is
[audit_k3_n2_finite_area_triangle_field_bound_independent_hostile.py](../scripts/audit_k3_n2_finite_area_triangle_field_bound_independent_hostile.py),
with separate tests in
[test_k3_n2_finite_area_triangle_field_bound_independent_hostile.py](../tests/test_k3_n2_finite_area_triangle_field_bound_independent_hostile.py).
