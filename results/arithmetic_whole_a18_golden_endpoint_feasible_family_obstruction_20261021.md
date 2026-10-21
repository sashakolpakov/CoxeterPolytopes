# A18 golden weak endpoint: the full two-wall family does not close

21 October2026. Exact manual supporting obstruction, pending full review.
Root proposed the feasible family; the source-audit author checked its
geometry and derived the convexity reduction and explicit counterprofile.
No jobs, enumeration, Lean, TeX, manuscript edit or dimension claim.

The context is the pending quadratic inventory
`arithmetic_whole_a18_quadratic_weak_cap_inventory_and_remaining_interface_20261021.md`
(SHA256 `4710025219a3d8a2d5a485de27d3a9f6962fa61289ef9e6120ee880ee4eec618`).
The counterprofile below is proved independently of that inventory's
unreviewed classification. Actual-point comparison and the N19 notation
are those of `compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`
(`3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`)
and `arithmetic_whole_a18_point52_actual_field_degree_reduction_20261021.md`
(`b57aa75adbde5c1a1b6f78354c58bba47d03b2047a923438b133be9aee5f2707`).

## The actual feasible family and its exact information

At a whole-A18 vertex use UNIT old normals, h_i=x.omega_i,
D_i=2i(19-i)/19, ||x||<1 and x.n_i>0. Suppose an ACTUAL isolated
dotted opposite is c=A p-gamma omega_i, gamma>1. Put u=gamma^(-2).
For every lambda in[1/gamma,gamma], xi=c+lambda*n_i is future timelike:

    xi^2=1-2gamma*lambda+lambda^2<0,
    xi.c=1-lambda*gamma<=0, xi.n_i=-gamma+lambda<=0.

Its time coefficient is A>0. Against every other ACTUAL wall each
summand has nonpositive product. It is therefore an actual feasible
point. Farthest comparison, division by gamma and v=lambda/gamma give

    h_i-v(x.n_i)>=T(v),
    T(v)=sqrt(D_i-u)-sqrt(2v-v^2-u), v in[u,1].       (1)

No desired height bound was supplied as a premise. The bare relation
is >=, including the two endpoints. With g(v)=sqrt(1-u-(v-1)^2),
g''(v)=-(1-u)/g(v)^3<0, so T is convex. Its chord through u,1
lies ABOVE T. The left side of(1) is affine in v. Consequently the
ENTIRE family(1) is equivalent to JUST its two endpoint inequalities:

    h_i-u(x.n_i)>=sqrt(D_i-u)-sqrt(u(1-u)),
    h_i-(x.n_i)>=sqrt(D_i-u)-sqrt(1-u).              (2)

Thus an interior parameter cannot strengthen both endpoint feet.

## An explicit counterprofile to the proposed relaxed exclusion

Take i5 and beta=phi^4=(7+3sqrt5)/2, so
u=4/beta=14-6sqrt5 satisfies u^2-28u+16=0. Exact polynomial tests give
.5835<u<.5836. All decimals here are rational numbers.

Define w0 by its old-dual heights: interpolate linearly from0 to
H6=2.346, then set

| i | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H_i | 2.346 | 2.454 | 2.524 | 2.558 | 2.558 | 2.524 | 2.454 | 2.346 | 2.195 |

and interpolate linearly from14 to19. Every listed height STRICTLY
exceeds sqrt(D_i)-.52 by squaring. The reference is concave and

    ||w0||^2=(1/2)[2.346^2/6+2(.108^2+.070^2+.034^2)
                         +.151^2+2.195^2/5]=.969566.

Set Y=omega5-u*n5-(5/6)omega6 and mu=.187. Direct dual products give

    Y.omega_j=0 (6<=j<=18), w0.Y=0,
    ||Y||^2=S=5/3-2u+u^2 in(.84,.841).

For x0=w0+mu*Y this yields

    ||x0||^2<.969566+.187^2(.841)=.998974929<.999,
    h5=1.955+.187(5/3-u)<2.157553,
    x0.n5=.187(1-u)>0.

All heights at6..18 retain their w0 values. The only changed old-normal
products are at4,5,6, where Y has products u/2,1-u,u/2-5/6.
The w0 product at6 is .1415; after the shift it is still>.0402.
All old-normal products of x0 are therefore nonnegative.

Write r=sqrt(D5-u). The rational u interval gives
r<2.6048, sqrt(u(1-u))>.4929 and sqrt(1-u)>.645. Hence

    T(u)<2.1119,
    x0.(omega5-u*n5)=1.955+.187S>2.11208;
    T(1)<1.9598,
    x0.(omega5-n5)=1.955+2(.187)/3>2.0796.

Both endpoints, hence the ENTIRE family, hold strictly.

To make every chamber sign strict, let d have heights d.omega_j=j(19-j),
so d=sum_(j=1)^18 omega_j, d.n_j=1 and ||d||^2=1140. Put
x=x0+10^(-6)d. Every old-normal product is now positive. The family
pressures increase by10^(-6)(70-v)>0, while all other central .52
floors remain strict. Since every x0 height is<3,

    ||x||^2<.999+2*10^(-6)(18*3)+1140*10^(-12)
           <.999109<1.

Also h5<2.158<sqrt(D5)-.52, so5 is still weak. The point lies on
the interior side of the proposed cap too: h5<sqrt(D5-u), equivalent
to gamma*h5<A. Reflection gives the identical obstruction at14.

This is a SAME-VECTOR counterprofile to a deduction from the complete
two-wall family, all other central .52 floors and strict old chamber
signs. It is NOT a compact ambient polytope, an actual complete wall
arrangement, or a counterexample to the intended whole-A18 exclusion.
It proves that paid-graft/KKT optimization using ONLY those hypotheses
cannot close this endpoint. A new actual-wall/helper constraint, a
stronger conditional floor elsewhere, or a global arithmetic/lattice
obstruction is necessary. No further research or computation is claimed.
