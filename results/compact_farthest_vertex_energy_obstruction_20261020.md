# Farthest-vertex energy and profile-adaptive weights

20 October 2026. Own paper proof. No computation, classification, new
uniform dimension ceiling, TeX build or Lean verification is asserted.

## 1. One vertex satisfies every positive weighted constraint

Let P be a compact hyperbolic Coxeter polytope. Choose ANY future unit
timelike vector z and choose an actual vertex p maximizing E=-z.p over
the finite vertex set. At p let n_i be its outward unit normals, B their
positive definite Gram matrix, omega_i their dual vectors, and

    d_i=(B^(-1))_ii,   u_i=omega_i/sqrt(d_i),
    tau_i=tanh(L_i/2)>0.

The adjacent actual vertex along edge i is

    p_i=cosh L_i*p-sinh L_i*u_i.

Maximality of E gives

    -z.p_i=E cosh L_i+sinh L_i*z.u_i<=E,
    tau_i<=x.u_i,   x=-(z-E*p)/E.                       (1)

Here x is spatial at p and ||x||^2=1-E^(-2)<1. Thus the SAME vertex p
satisfies, for EVERY nonzero vector of nonnegative local weights w,

    sum_i w_i tau_i <= x.(sum_i w_i u_i)
                       < ||sum_i w_i u_i||.           (2)

The strict inequality follows from ||x||<1 and the fact that the u_i
are a basis. No weight compatibility between adjacent vertices is
needed in this argument. Weights may depend on the vertex type and the
actual profile/length data, provided the resulting certificate is
proved at every possible vertex type. This does not assert that every
vertex can be the farthest vertex for some z.

The common-white weighted undirected identity in
`results/compact_finite_edge_white_weighted_drift_20261020.md` remains
valid, but its shared-edge weights are a sufficient global mechanism,
not an extra restriction on the pointwise farthest-vertex test.

## 2. Exact pointwise energy

Let K be the Gram matrix of the unit dual directions:

    K_ij=(B^(-1))_ij/sqrt(d_i d_j).

Define the covering energy

    H_K(tau)=min { ||y||^2 : y.u_i>=tau_i for every i }.

The minimum exists and is unique: the feasible set is nonempty (the
u_i form a basis), closed and convex, and the squared norm is coercive
and strictly convex. Equation (1) implies the necessary condition

    H_K(tau)<1                                           (3)

at the chosen farthest vertex. For any w>=0, (2) gives the elementary
certificate

    (w.tau)^2/(w^T K w)<1.

If K^(-1)tau>=0, take w=K^(-1)tau. The exact energy is then

    H_K(tau)=tau^T K^(-1)tau.

Indeed the vector with equalities y.u_i=tau_i gives this norm, and the
displayed nonnegative weight gives the matching lower bound by
Cauchy--Schwarz. When some coefficients are negative one must retain
the inequality constraints and choose the active subset; the same
unconstrained expression is not automatically the minimum.

If the finite vertex diagram is a forest of orthogonal components,
their spatial spans are orthogonal. Feasibility separates by component,
and the exact energies ADD. Consequently any proved componentwise
lower estimates, including compatible exceptional-profile counts, may
be summed before comparison with1.

## 3. A rank-two block formula

For I2(m), put c=cos(pi/m). The unit dual Gram is

    K=[[1,c],[c,1]].

Its covering energy at half-tangents t_1,t_2>0 is exactly

    t_2^2,                               if t_1<=c*t_2;
    t_1^2,                               if t_2<=c*t_1;
    (t_1^2+t_2^2-2c*t_1*t_2)/(1-c^2),   otherwise.     (4)

In the first branch y=t_2*u_2 is feasible and already has the required
norm. The second branch is symmetric. In the remaining branch
K^(-1)t has nonnegative coefficients, so Section2 applies. In
particular, for all cases,

    H_K(t_1,t_2)>=(t_1+t_2)^2/[2(1+c)].

For A2 the interior formula is (4/3)(t_1^2+t_2^2-t_1*t_2), with branch
ratios1/2 and2. For an A1 component the energy is simply t_i^2.

## 4. A strict consequence already available

The actual-wall argument in
`results/compact_right_corner_short_edge_capacity_20261020.md` proves
that an A1^n corner has at most two directions with

    L_i<L_*=asinh(1/sqrt2),
    tau_*=tanh(L_*/2)=sqrt3-sqrt2.

At a farthest A1^n vertex, the dual directions are orthonormal, so (3)
becomes sum_i tau_i^2<1. On the other hand the capacity theorem gives

    sum_i tau_i^2 >= (n-2)(5-2sqrt6).

For n>=12 this exceeds1: already
10(5-2sqrt6)>1 follows from49^2>6*20^2 (difference1).
Therefore a farthest vertex in a compact Coxeter polytope of dimension
at least12 CANNOT have diagram A1^n. This improves the shared-linear-
weight witness threshold14 for this one family, with no new support
classification.

This is a pointwise vertex obstruction, not a compact dimension
theorem. Other finite vertex diagrams must still be excluded at the
SAME farthest-vertex test. The explicit arbitrary-rank short-edge star
in `results/connected_arithmetic_lanner_triangle_rank_barrier_20261020.md`
still prevents replacing full compatible profile analysis with a
blanket prohibition of one short A_n edge.
