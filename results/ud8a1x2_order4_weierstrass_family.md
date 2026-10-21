# An explicit D8+A1^2 order-four carrier for the H10/15 section

## Positive theorem

Let \(s=t^2\). For parameters \(r,\alpha,\beta,\gamma\), put

\[
\begin{aligned}
p_0&=-3r^2+3-\alpha-\beta,\\
q_1&=9r-2\alpha r-\beta r-r^3-8+3\alpha-2\gamma,\\
q_0&=2r^3-2+\alpha-\gamma-q_1,
\end{aligned}
\]

and

\[
\begin{aligned}
P(s)&=-3s^3+\alpha s^2+\beta s+p_0,\\
Q(s)&=2s^4-\alpha s^3+\gamma s^2+q_1s+q_0.
\end{aligned}
\]

Let \(X_{r,\alpha,\beta,\gamma}\) be the minimal resolution of

\[
 y^2=x^3+P(t^2)x+tQ(t^2),
\]

with

\[
 \sigma(x,y,t)=(-x,iy,-t).
\]

On a nonempty Zariski-open parameter set, the singular fibers are

\[
 I_4^*+2I_2+10I_1.
\]

For a very general member,

\[
 NS(X)=U\oplus D_8(-1)\oplus A_1(-1)^2.
\]

The primitive non-symplectic order-four automorphism fixes every \(D_8\)
component and exchanges the two \(A_1\) components. Thus its Picard action
is the exact action of Brandhorst--Hofmann row 0.4.4.2. After the already
certified Kähler/Weyl chamber transport, this gives an explicit
algebro-geometric realization of R12-D8A1x2-C2a-H10-15, a finite-volume
noncompact Coxeter polytope in \(H^{10}\) with 15 facets.

The abstract-polytope lineage is still only bounded-corpus novelty; this
theorem claims a new explicit algebraic realization in the project, not
worldwide priority for the polytope or the ambient order-four K3 family.

## The simultaneous fiber constraints

Parity makes the equation invariant:

\[
 \sigma^2(x,y,t)=(x,-y,t),\qquad
 \sigma^*(dt\wedge dx/y)=-i\,dt\wedge dx/y.
\]

The coefficients solve two independent inverse fiber-design problems.

At \(t=\infty\), use \(u=1/t\), \(X=u^4x\), and \(Y=u^6y\). Then

\[
 Y^2=X^3+u^2\bar P(u^2)X+u^3\bar Q(u^2),
\]

where

\[
\begin{aligned}
\bar P(z)&=-3+\alpha z+\beta z^2+p_0z^3,\\
\bar Q(z)&=2-\alpha z+\gamma z^2+q_1z^3+q_0z^4.
\end{aligned}
\]

The constant and linear terms of
\(4\bar P(z)^3+27\bar Q(z)^2\) vanish identically, while its
\(z^2\)-coefficient is

\[
 -9(\alpha^2-12\beta-12\gamma).
\]

Away from this coefficient's zero locus, the orders are

\[
 (\operatorname{ord}_u a_4,\operatorname{ord}_u a_6,
   \operatorname{ord}_u\Delta)=(2,3,10),
\]

so the fiber is exactly \(I_4^*\), with root lattice \(D_8(-1)\).

At \(s=1\), the definitions give

\[
 P(1)=-3r^2,\qquad Q(1)=2r^3,\qquad
 Q'(1)=-rP'(1)-r^3.
\]

Consequently

\[
 D(s)=4P(s)^3+27sQ(s)^2=(s-1)^2R(s),
\]

with \(\deg R=5\). For \(r\ne0\) and off one further hypersurface, the two
fibers over \(t=\pm1\) are exactly \(I_2\). Their singular cubics are

\[
\begin{array}{c|c}
t=1 &(x-r)^2(x+2r),\\
t=-1 &(x+r)^2(x-2r),
\end{array}
\]

so \(\sigma\) visibly exchanges the two resolved \(A_1\) components. The
remaining five \(s\)-roots lift to ten \(I_1\) fibers, and

\[
 10+2\cdot2+10=24.
\]

The desired open set is nonempty. At
\((r,\alpha,\beta,\gamma)=(1,0,1,0)\),

\[
 P=-3s^3+s-1,\qquad Q=2s^4-s+1,
\]

and

\[
 R=108s^5-36s^3+31s-4.
\]

Its discriminant is \(6477927509458944\ne0\), \(R(0)=-4\ne0\), and

\[
 \operatorname{Res}(R,P)=-1440747\ne0.
\]

So the residual roots are simple, nonzero, and do not collide with a zero
of \(a_4\).

## Why the D8 action is the identity

This is the local gate that cannot be replaced by a rank count: \(D_8\) has
an outer involution which could exchange its two spinor nodes.

In the coordinates above,

\[
 \sigma(X,Y,u)=(-X,iY,-u).
\]

Here is the calculation in the actual Weierstrass germ.  Translate the
double Tate root by \(V=X-u\). Direct expansion gives

\[
\begin{aligned}
Y^2={}&V^2(V+3u)+\alpha u^4V+\beta u^6V+p_0u^8V\\
     &+(\beta+\gamma)u^7+(p_0+q_1)u^9+q_0u^{11}.
\end{aligned}
\]

This identity is checked symbolically from the displayed global
Weierstrass coefficients. The action is exactly

\[
 (Y,V,u)\longmapsto(iY,-V,-u),
\]

and the right-hand side is anti-invariant, as required by
\((iY)^2=-Y^2\).

Now resolve the branch plane in its actual charts. On the first blowup put
\(q_0=V/u\). Then \(q_0\) is fixed and the strict branch on \(u=0\) is

\[
 q_0^2(q_0+3)=0.
\]

The transverse branch is separated at the individually fixed point
\(q_0=-3\), while the two tangent branches retain the fixed center
\((u,q_0)=(0,0)\). On the next chart put \(q_1=V/u^2\). Here
\(q_1\mapsto-q_1\), and the special strict branch is \(3q_1^2=0\), so its
common center remains fixed. On the third chart put \(q_2=V/u^3\). Now
\(q_2\) is fixed and the two tangent branches separate at the roots of

\[
 3q_2^2+\alpha q_2+\beta+\gamma=0.
\]

Their discriminant is

\[
 \alpha^2-12\beta-12\gamma\ne0,
\]

exactly the open coefficient forcing the \(I_4^*\) fiber. Thus both
separation points are fixed individually. Every blowup center, every
exceptional divisor, and every separated branch point in the actual
double-cover resolution is invariant. Hence the exceptional \(D_8\)
diagram is fixed vertex by vertex and its spinor nodes are not exchanged.

## Exact generic NS and Mordell--Weil saturation

The fiber and zero section together with the reducible-fiber components give

\[
 L=U\oplus D_8(-1)\oplus A_1(-1)^2\subset NS(X),
\]

of rank 12 and determinant \(-16\).

The parameters are recovered exactly from the normalized coefficients:

\[
 r=-\frac{3Q(1)}{2P(1)},\qquad
 \alpha=[s^2]P,\qquad
 \beta=[s]P,\qquad
 \gamma=[s^2]Q.
\]

The first identity is valid on the defining open \(r\ne0\), since
\(P(1)=-3r^2\) and \(Q(1)=2r^3\). The Jacobian of the recovered
\((r,\alpha,\beta,\gamma)\) is exactly \(I_4\).

The unique \(I_4^*\) fiber fixes infinity. A base automorphism commuting
with \(t\mapsto-t\) is then \(t\mapsto ct\), and preserving
\(\{1,-1\}\) forces \(c^2=1\). A short-Weierstrass homothety preserving the
leading Tate coefficients \((-3,2)\) satisfies
\(\lambda^4=\lambda^6=1\), hence \(\lambda^2=1\). Thus the normalized
coordinate stabilizer is finite.

The marked fibration and section determine integral classes
\(F,O\in NS(X)\) with \(F^2=0\), \(O^2=-2\), and \(F\cdot O=1\).
They are locally constant on a connected fiber of the forgetful map. Fixed
\(F,O\) recover the pencil \(|F|\) and section \(O\), up to the finite
coordinate stabilizer just computed. Since \(H^0(X,T_X)=0\), K3
automorphisms are discrete as well. Hence the marked-to-unmarked map has
zero-dimensional fibers and this is genuinely a four-dimensional K3
family.

Since \(T\) is a \(\mathbf Q(i)\)-module, \(\rho\) is even. A strict
generic jump from \(\rho\ge12\) would have \(\rho\ge14\), whose primitive
order-four period ball has dimension at most three. This contradicts the
four-dimensional family, so the very general Picard rank is 12.

Rank alone does not settle saturation because \(L\) has determinant 16.
Use the Shioda height formula. A nonzero component class in the \(D_8\)
discriminant has local correction at most 2, and each \(I_2/A_1\)
correction is at most \(1/2\). Thus a nonzero torsion section would satisfy

\[
 0=h(P)=4+2(P\cdot O)-\sum_v\operatorname{contr}_v(P)
       \ge 4-(2+1/2+1/2)=1,
\]

a contradiction. Mordell--Weil torsion is trivial. Shioda--Tate gives rank
zero, so the trivial lattice is saturated and

\[
 NS(X)=L
\]

for a very general member.

The rational witness proves the fiber open is nonempty but is not itself
asserted to have Picard rank 12 without a separate arithmetic reduction
computation.

## Exact BH action and Coxeter handoff

Let \(R_+,R_-\) be the two \(I_2\) roots. Since \(D_8\) is fixed and these
roots are swapped,

\[
\begin{aligned}
NS(X)^\sigma&\cong U+D_8(-1)+\langle-4\rangle,\\
(NS(X)^\sigma)^\perp&\cong\langle-4\rangle.
\end{aligned}
\]

The fixed lattice has rank 11, signature \((1,10)\), determinant 16, and
nontrivial Smith factors \(2,2,4\). Its discriminant quadratic histogram is

\[
 \{0:4,\ 1:4,\ 3/4:2,\ 7/4:6\}.
\]

The coinvariant is rootless. The fixed-plus-coinvariant lattice has the
literal index-two glue
\((R_++R_-+R_+-R_-)/2=R_+\). These data match the pinned
Brandhorst--Hofmann 0.4.4.2 certificate, which proves integral conjugacy of
the full Néron--Severi action pairs and transports its invariant Kähler
class to Roulleau's complete nef chamber.

The full characteristic polynomial is

\[
 (\lambda-1)^{11}(\lambda+1)(\lambda^2+1)^5.
\]

The BH fixed-locus fingerprint is

\[
\begin{array}{c|c}
\sigma &8\text{ isolated points}+2\text{ rational curves},\\
\sigma^2&4\text{ rational curves}+1\text{ genus-3 curve}.
\end{array}
\]

The square's second line also follows from the 2-elementary invariants
\((r,a)=(12,4)\).

The independent literal Vinberg audit gives:

| invariant | value |
|---|---:|
| hyperbolic dimension | 10 |
| facets | 15 |
| finite vertices | 42 |
| cusps | 10 |
| elliptic edges | 294 |
| failed edge/parabolic audits | 0 |

Thus this is an explicit realization of a **finite-volume noncompact**
hyperbolic Coxeter polytope.

The canonical digest of the complete finite-vertex, ideal-component, edge,
and parabolic audit is
00977c1aa97e62df0ac6bbd22e635cd3b37dbf38e22b3cf52770ecd7a42c15fb.

## Lineage and reproduction

The parity construction lies inside
[Artebani--Sarti Example 6.3](https://arxiv.org/abs/1102.4436). The present
normal form is the inverse-designed specialization which adds the
\(I_4^*/D_8\) fiber by a second-order leading-discriminant cancellation and
adds \(A_1^2\) by an exchanged finite double root.

Artifacts:

- [certificate](../scripts/certify_ud8a1x2_order4_weierstrass_family.py)
- [machine archive](../data/enumerations/ud8a1x2_order4_weierstrass_family.json)
- [tests](../tests/test_ud8a1x2_order4_weierstrass_family.py)

Rebuild and test with:

    PYTHONPATH=.:src python3 scripts/certify_ud8a1x2_order4_weierstrass_family.py --check-archive
    PYTHONPATH=.:src pytest -q tests/test_ud8a1x2_order4_weierstrass_family.py

For an external replay, add
--brandhorst-hofmann-zip /path/to/K3Groups.zip to the certificate command.
