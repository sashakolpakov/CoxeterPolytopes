# Farthest dodecahedral facets and a sharp wall-radius inequality

20 October 2026. Written proof awaiting independent review. No jobs,
numerical experiments, classification census or formal verification.

Let P be any compact right-angled hyperbolic four-polytope. For an
interior point q and each outward unit facet normal n, put

    t_n(q) = -<q,n> = sinh dist(q,wall(n)),
    T(q) = max_n t_n(q),       phi=(1+sqrt(5))/2.

The distance is to the supporting hyperplane. This convention matters;
it does not assert that the nearest point of q on every supporting
hyperplane lies in its facet.

**Proposed theorem.** Every facet attaining T(q) is a regular
right-angled dodecahedron. For every interior q,

    T(q)^2 >= phi^3.

Equality occurs if and only if P is the regular right-angled 120-cell
and q is its center. In particular the farthest-wall height function
of an arbitrary P is determined pointwise by its dodecahedral facets.

This is a universal metric restriction and equality recognition, not
a volume comparison. In particular no inscribed 120-cell or universal
minimum-volume theorem is inferred.

## 1. Actual projection and the strict spherical code

Use the Lorentz form (-,+,+,+,+). Distinct actual wall normals have
nonpositive products. Let n attain T=T(q)>0 and define

    q'=(q+Tn)/sqrt(1+T^2).

It is a future unit timelike point on n's wall. For every other wall m,

    <q',m>=(<q,m>+T<n,m>)/sqrt(1+T^2)<0.

Consequently q' lies in the relative interior of the actual facet F.
The intrinsic outward unit normals of F are precisely those parent
normals m_i adjacent to F, all orthogonal to n. Their heights satisfy

    0<t'_i=-<q',m_i>=t_{m_i}(q)/sqrt(1+T^2)
             <=T/sqrt(1+T^2)<1.

In the intrinsic H3 write m_i=(t'_i,sqrt(1+t_i'^2)u_i), with u_i unit
in R3. Acuteness gives, for i distinct from j,

    u_i.u_j <= t'_i t'_j/sqrt((1+t_i'^2)(1+t_j'^2)) <1/2.

The directions are distinct, since equal directions would make the
parent scalar product strictly positive.

## 2. An elementary quartic gives at most thirteen directions

This deliberately avoids importing the three-dimensional kissing
number twelve. Set

    f(x)=(x+1)(x+2/5)^2(x-1/2).

It is nonpositive on [-1,1/2]. Its Legendre expansion in dimension
three is

    f = (7/50) P0 + (23/50) P1 + (107/175) P2
                        + (13/25) P3 + (8/35) P4,
    f(1)=49/25.

All coefficient signs are positive. Here

    P0=1, P1=x, P2=(3x^2-1)/2,
    P3=(5x^3-3x)/2, P4=(35x^4-30x^2+3)/8.

For completeness their required positivity is elementary. If u is unit
in R3, define the trace-free symmetric tensors

    A2(u)=u^2-I/3,
    A3(u)=u^3-sym_3(u I)/5,
    A4(u)=u^4-sym_6(u^2 I)/7+sym_3(I I)/35.

The symmetrizations are unnormalized sums over the indicated distinct
positions: three terms u_a delta_bc for rank three, six terms
u_a u_b delta_cd, and three terms delta_ab delta_cd for rank four.
Tracing the last tensor gives

    u^2 - (7u^2+I)/7 + (5I)/35 = 0.

Because each A_j(v) is trace-free, contracting it with A_j(u) is the
same as contracting with u^j. Their scalar kernels are respectively

    x^2-1/3,       x^3-3x/5,       x^4-6x^2/7+3/35.

These are positive multiples of P2, P3 and P4. The P1 kernel is the
ordinary scalar product. Thus each sum over ordered pairs of one
of P1,...,P4 is nonnegative, being a squared tensor-sum norm times a
positive constant. For N unit directions with pair products <1/2,

    (7/50)N^2 <= sum_ij f(u_i.u_j) <= N f(1)=(49/25)N,

so N<=14. If N=14, equality forces every off-diagonal f(u_i.u_j)=0.
Strictness excludes the root 1/2; the remaining roots are -1 and -2/5.
All pair products would be <=-2/5. But then

    ||sum_i u_i||^2 <= N-(2/5)N(N-1)<0,

an impossibility. Hence N<=13.

## 3. Thirteen faces cannot occur here

A compact right-angled three-polyhedron has trivalent finite vertices,
and every face is a polygon with at least five sides. Euler and
incidence imply at least twelve faces. The independently proved
polygon capacity lemma says that a k-gonal face requires at least
2k+2 total faces. Its proof counts the face, its k neighbors, the k
distinct perpendicular-edge opposites, then excludes equality 2k+1
by an ear in the complementary dual disc.

Thus a right-angled three-polyhedron with at most thirteen faces can
have no face with six or more sides. Every face is pentagonal, and
Euler gives exactly twelve faces. Its metric is the regular
right-angled dodecahedron D, by the explicitly stated three-dimensional
twelve-face recognition input in the current compact paper.

Apply this to the actual F from Section 1. This proves the farthest
facet assertion for every maximizing n, not only for one selection.
In particular, at every interior q a dodecahedral facet exists.

## 4. The intrinsic dodecahedral constant

In intrinsic H3 place the center o of D at (1,0). Its twelve wall
normals have form

    m_i=(h,sqrt(1+h^2) r_i),

where the r_i are the unit vertices of the regular icosahedron. Adjacent
ones have scalar product 1/sqrt(5), and the corresponding walls are
orthogonal. Hence

    h^2/(1+h^2)=1/sqrt(5),
    h^2=1/(sqrt(5)-1)=phi/2.

The directions occur in antipodal pairs. For any future unit point
u=(a,w) in D's wall, a>=1. The two heights for an opposite pair sum
to 2ah. Therefore

    max_i(-<u,m_i>) >= ah.

Equality a=1 is equivalent to u=o. This is only a bound on the maximum,
not a claim that each individual height is at least h.

## 5. Sharp radius and equality propagation

Take the maximizing facet F at q and its actual foot q'. It is D.
Let a be the time coordinate of q' relative to D's center. Combining
the last inequality with Section 1 gives

    a sqrt(phi/2) <= max_i t'_i <= T/sqrt(1+T^2).

Since a>=1,

    T^2 >= (phi/2)/(1-phi/2) = phi^3.

If equality T^2=phi^3 holds, both inequalities are equalities and a=1.
The actual foot is the center of F. All twelve intrinsic heights are
then sqrt(phi/2), so every neighboring parent wall has height exactly

    sqrt(phi/2)*sqrt(1+T^2)=T.

Each neighboring facet is therefore also a maximizing facet at the
same original q, hence is D. The same equality argument applies anew
to each of them. Connectedness of parent facet adjacency propagates
this to every facet. The all-D continuation proof identifies P with T,
the regular right-angled 120-cell. Once its matched facet and its
centered foot are fixed, q is the unique inward normal point at the
displayed distance; by symmetry it is the center of T.

Conversely at the center of T, symmetry makes all parent heights equal
to T0. Projection to a facet is its center, and

    sqrt(phi/2)=T0/sqrt(1+T0^2)

gives T0^2=phi^3. Thus the equality is attained exactly as claimed.

## 6. Useful boundary of the argument

A non-dodecahedral facet can never be a farthest facet from any
interior point. It does not follow that the facet is redundant or
that all facets are dodecahedra: a facet need not be farthest from any
interior point. Likewise the maximum-height inequality does not imply
a lower bound on the minimum height, and does not produce an inscribed
ball of that radius. Either invalid inference would falsely settle
the unrestricted volume problem.

For a facet known to be D and vertices v,w joined by one of its edges,
the outside endpoint wall normals have forms

    r_v=sinh(s_v)v-cosh(s_v)n,
    r_w=sinh(s_w)w-cosh(s_w)n.

Since -<v,w>=cosh L=phi, their acute product only implies

    coth(s_v)coth(s_w)<=phi.

It does not imply s_v>=L separately: a shorter edge may be compensated
by longer edges at its neighbors. This is the exact obstruction to
claiming a whole 120-cell from a single dodecahedral facet.

Source interfaces (unchanged):

* manuscripts/compact_ra4_volume_constraints.tex,
  SHA256 afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f:
  compact normal gap, polygon capacity, twelve-face metric input,
  all-D recognition.
* results/right_angled_dimension8_descending_height_proof_20261020.md,
  SHA256 9470fa1fa3a973a12c412fbc1257537655536dba4ffc4e9ac23c320615babdbe:
  actual maximum-height projection and direction construction.

Novelty is not asserted: a targeted primary-literature comparison of
this particular farthest-facet/radius statement remains to be done.
The argument and precise scope above are independent written proofs,
not an inference from the number of reviews of earlier results.
