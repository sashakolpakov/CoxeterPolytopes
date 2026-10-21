# Farthest-point pressure from actual spherical faces

20 October 2026. Paper derivation only, in discussion with root. No
enumeration, numerical optimization, proof job or dimension claim.
This specializes the actual mixed-projection lemma, SHA256
`3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`,
in `compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`.

## 1. A general inequality retaining actual finite faces

All normals below are UNIT outward normals of actual facets. Let p be
an actual vertex farthest from a fixed strict interior point
`z=E(p-x)`, where `x in p-perp` and `||x||<1`. Thus `x.n_j>0` for
each old incident normal n_j. For an actual outside wall write

    c=a p-v_c,    a>0,    ||v_c||^2=1+a^2.

Choose any old subdiagram J for which the ACTUAL normals J+c have a
positive definite Gram. With B_J the old unit Gram, put

    beta_j=-c.n_j>=0,
    w_J=sum_j (B_J^-1 beta)_j n_j,
    delta=1-beta^T B_J^-1 beta>0.

The inverse has nonnegative entries. Hence w_J has nonnegative old
normal coefficients, and `x.w_J>0` whenever beta is nonzero. The
orthogonalized wall `c+w_J` has norm delta and the projection of p to
the common finite face is

    q=[p-v_J]/sqrt(1-R_J),
    v_J=a(v_c-w_J)/(a^2+delta),
    R_J=a^2/(a^2+delta).

This is an ACTUAL feasible point: the positive acute projection
argument applies to the selected actual facets. No intrinsic Coxeter
property of arbitrary metric faces is assumed. Farthest-point
comparison therefore gives the exact additional constraint

    x.(v_c-w_J) >= H(a,delta),
    H(a,delta)=a/[1+sqrt(delta/(a^2+delta))].           (1)

In particular `x.v_c>H(a,delta)` when beta is nonzero. The vector
inequality in (1) is stronger than discarding the positive correction.
It applies to nonisolated caps too; there v_c is the complete spatial
projection, not only the released contact's dual direction.

## 2. Correct isolated-contact floor

For an isolated contact i put `D=||omega_i||^2`,
`u_i=omega_i/sqrt(D)`, beta=cos(theta), and

    v_c=beta sqrt(D) u_i,
    a^2=beta^2 D-1,
    r=tanh L_i=sqrt(1-1/(beta^2 D)).

Taking J={i}, with delta=sin^2(theta), (1) implies

    x.u_i > r/[1+tan(theta)/sqrt(D-1)].               (2)

Its right side S has square

    S^2=(D-1)/D *
        [sqrt(D-1)-tan(theta)]/[sqrt(D-1)+tan(theta)].

For D>=6 and theta<=pi/5 this exceeds the correct long-edge floor

    tau_long^2=(sqrt(D)-1)/(sqrt(D)+1).

Indeed the comparison reduces to `tan(theta)<F(sqrt(D))`, where

    F(q)=(2q+1)sqrt(q^2-1)/(2q^2+2q+1).

F is increasing for q>1: the derivative sign after squaring reduces
to `2q(4q^3+18q^2+18q+5)>0`. At q=sqrt(6),
`F(q)^2>8/15>5-2sqrt(5)=tan^2(pi/5)`; the first inequality is
`331>116sqrt(6)`, and the second follows from `sqrt(5)>67/30`.

This is a bound on x.u_i, not an assertion that the actual short edge
itself has length at least arcosh(sqrt(D)). In particular the floor
is a square-root inequality for tau, not the displayed expression
incorrectly interpreted as tau itself.

## 3. Much smaller H4 residuals, including nonisolated caps

At a whole A_n vertex use its unit chain normals. Suppose an actual
outside wall has contact beta=cos(pi/5) with n_i, and is orthogonal to
`n_(i+1),n_(i+2)`. The three old nodes on that side together with c
form the ACTUAL spherical H4 diagram. Consequently

    delta_5=1-(3/2)beta^2=(7-3sqrt(5))/16,
    w_+=beta[(3/2)n_i+n_(i+1)+(1/2)n_(i+2)].

Equation (1) holds with this very small delta_5. Any contacts of c
outside these three old nodes remain in v_c and a; they are not set
to zero. Thus this applies, for example, to a released label5 with
another contact far down a white arm.

If c is also orthogonal to `n_(i-1),n_(i-2)`, there are TWO actual
H4 projections. Put

    w_- =beta[(3/2)n_i+n_(i-1)+(1/2)n_(i-2)],
    y_bar=v_c-(w_++w_-)/2.

Averaging (1), and using the exact unit-chain products, gives

    x.y_bar >= H(a,delta_5),
    ||y_bar||^2=a^2+1-(15/8)beta^2,
    H(a,delta_5)^2 < a^2+1-(15/8)beta^2.             (3)

The last inequality follows from `||x||<1`. It is valid for isolated
and nonisolated walls satisfying the stated local orthogonality;
it does not assume away any far white contacts. Its exact ingredients
are `v_c.w_+=v_c.w_-=(3/2)beta^2`,
`||w_+||^2=||w_-||^2=(3/2)beta^2`, and
`w_+.w_-=(3/4)beta^2`.

The useful remaining obligation is to combine the vector demands
(1), or their two-sided form (3), across a COMPLETE collection of
actual opposites. A single scalar edge floor still loses this data.
For example the previously recorded formal A14 relaxation already
meets every correct central long-edge floor: its central tau_i^2
are at least `149769/280000>1/2`, whereas the floor squares are
less than 1/2 since D_i<=15/2<9. That formal data is not a native
Coxeter configuration. Neither it nor the present projection
inequalities settle whole-A_n exclusion or change a dimension bound.
