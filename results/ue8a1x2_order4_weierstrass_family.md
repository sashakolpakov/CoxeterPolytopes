# An explicit order-four Weierstrass carrier for the `H10/12` section

## Positive theorem

Put (s=t^2), and define

\[
\begin{aligned}
A(s)&=-3+\alpha(s-1)+\beta(s-1)^2,\\
C(s)&=2-(\alpha+1)(s-1)+\gamma(s-1)^2+\delta(s-1)^3.
\end{aligned}
\]

Let (X_{\alpha,\beta,\gamma,\delta}) be the minimal resolution of

\[
 y^2=x^3+A(t^2)x+tC(t^2),
\]

and let

\[
 \sigma(x,y,t)=(-x,iy,-t).
\]

On a nonempty Zariski-open set of parameter space, the singular fibers are

\[
 II^*+2I_2+10I_1.
\]

For a very general member,

\[
 NS(X)=U\oplus E_8(-1)\oplus A_1(-1)^2.
\]

The automorphism (sigma) is primitive non-symplectic of order four.  It is
the identity on (U+E_8(-1)) and exchanges the two (A_1) fiber components.
This is the exact Picard action of Brandhorst--Hofmann row `0.4.4.1`.
After the already certified Kähler/Weyl chamber transport, the family is
therefore an explicit algebro-geometric realization of
`R12-E8A1x2-C2-H10-12`: a finite-volume noncompact Coxeter polytope in
(H^{10}) with 12 facets.

This is a new equation-level Coxeter realization in this project.  The
abstract polytope is known, and no worldwide priority for the underlying
order-four K3 family is claimed.

## Why the fibers are forced by the equation

The equation has the required parity:

\[
 A((-t)^2)=A(t^2),\qquad (-t)C((-t)^2)=-tC(t^2).
\]

Thus the image of its right-hand side under
((x,t)\mapsto(-x,-t)) is its negative, exactly matching
((iy)^2=-y^2).  Moreover

\[
 \sigma^2(x,y,t)=(x,-y,t),
\]

so (sigma) has order four.  On the holomorphic two-form,

\[
 \sigma^*(dt\wedge dx/y)=-i\,dt\wedge dx/y.
\]

It is therefore primitive non-symplectic.

Write (u=s-1).  The reduced discriminant is

\[
 D(s)=4A(s)^3+27sC(s)^2=u^2Q(u),
\]

where (Q) has degree five when (delta\ne0), and

\[
 Q(0)=-9\bigl(\alpha^2+6\alpha-12\beta-12\gamma+9\bigr).
\]

Hence, away from the displayed hypersurface, (t=1) and (t=-1) are
exactly type-(I_2) fibers.  This is also visible directly:

\[
\begin{array}{c|c}
t=1 & x^3-3x+2=(x-1)^2(x+2),\\
t=-1 & x^3-3x-2=(x+1)^2(x-2).
\end{array}
\]

The automorphism exchanges the two singular points and therefore exchanges
the two resolved (A_1) components.

After homogenizing over ([T:U]\in\mathbf P^1), the Weierstrass
coefficients have orders at least four and exactly five at (U=0), while
the discriminant has order ten.  Thus the fiber at infinity is (II^*),
with root lattice (E_8(-1)).  The (E_8) graph has no nontrivial graph
automorphism, so all of its component classes are fixed by (sigma).

The remaining quintic (Q(s-1)) gives five values of (s); on the generic
open set each has two distinct nonzero square roots in (t), giving ten
(I_1) fibers.  The Euler check is

\[
 10+2\cdot2+10\cdot1=24.
\]

The open set is provably nonempty.  At
((\alpha,\beta,\gamma,\delta)=(0,0,0,1)),

\[
 Q(u)/27=u^5+u^4-2u^3+2u^2+5u-3.
\]

This quintic has discriminant (283904\ne0), and its value at
(s=0), i.e. (u=-1), is (-4\ne0).  Since (A=-3) in this witness,
no residual discriminant root is additive.  Squarefreeness and the
noncollision conditions therefore define a nonempty Zariski-open locus with
the stated fiber configuration.

## Exact generic Néron--Severi lattice

The fiber, zero section, and nonidentity components of the reducible fibers
give

\[
 L=U\oplus E_8(-1)\oplus A_1(-1)^2\subseteq NS(X),
\]

of rank 12 and determinant (-4).

The normal form has four genuine moduli, as follows from an explicit inverse
to its coefficient map.  With derivatives taken in \(s\) at \(s=1\),

\[
 \alpha=A'(1),\qquad
 \beta=\frac{A''(1)}2,\qquad
 \gamma=\frac{C''(1)}2,\qquad
 \delta=\frac{C'''(1)}6.
\]

The Jacobian of these four recovered quantities with respect to
\((\alpha,\beta,\gamma,\delta)\) is exactly \(I_4\).

There is no hidden continuous coordinate quotient.  The unique \(II^*\)
fiber fixes infinity.  A base automorphism commuting with \(t\mapsto-t\)
and fixing infinity is \(t\mapsto ct\); preserving the normalized pair
\(\{1,-1\}\) forces \(c^2=1\).  A short-Weierstrass homothety preserving
\(A(1)=-3,C(1)=2\) satisfies
\(\lambda^4=\lambda^6=1\), hence \(\lambda^2=1\).  The normalized
coordinate stabilizer is therefore finite.

Finally, forgetting the marked fibration and section has
zero-dimensional fibers, not merely countable fibers.  They determine
integral classes \(F,O\in NS(X)\) with

\[
 F^2=0,\qquad O^2=-2,\qquad F\cdot O=1.
\]

Integral classes are locally constant on a connected fiber of the
forgetful map.  For fixed \(F,O\), the pencil \(|F|\) and the curve \(O\)
recover the elliptic fibration and its section up to the already-removed
base coordinates.  Moreover \(H^0(X,T_X)=0\), so the automorphism group of
a K3 surface is discrete.  Thus no positive-dimensional forgetful fiber
remains, and the underlying K3/automorphism family is four-dimensional.

Because (sigma) is primitive of order four, its transcendental lattice is
a (mathbf Q(i))-module.  Hence both (operatorname{rk}T) and
(ho=22-operatorname{rk}T) are even.  We already have (ho\ge12).
If the very general Picard rank were larger, it would be at least 14, so
(operatorname{rk}T\le8).  The primitive-eigenvalue period ball would then
have dimension at most (8/2-1=3), contradicting the four-dimensional
family.  Therefore the very general Picard rank is exactly 12.

There is no saturation ambiguity.  A proper same-rank even overlattice of
(L) would have index two and determinant (-1), hence would be even
unimodular of signature ((1,11)).  Such a lattice cannot exist because
(1-11) is not divisible by eight.  Consequently

\[
 NS(X)=L
\]

for a very general member.  Shioda--Tate then gives Mordell--Weil rank zero,
and the same saturation argument excludes torsion.

This proves a geometric-generic theorem.  It does not assert that the
particular rational witness above has Picard rank 12; certifying one named
arithmetic member by good reduction would be an optional refinement.

## Exact Picard action and glue

Let (R_+) and (R_-) denote the two (A_1) roots.  The action is

\[
 R_+\longleftrightarrow R_-.
\]

Thus

\[
\begin{aligned}
NS(X)^\sigma&=U+E_8(-1)+\langle R_++R_-\rangle
             \cong U+E_8(-1)+\langle-4\rangle,\\
(NS(X)^\sigma)^\perp&=\langle R_+-R_-\rangle
             \cong\langle-4\rangle.
\end{aligned}
\]

The coinvariant is rootless.  The fixed lattice has rank 11, signature
((1,10)), determinant 4, Smith factors (1^{10},4), and discriminant
quadratic histogram

\[
 \{0:1,\ 1:1,\ 7/4:2\}.
\]

The primitive glue is literal.  If (f=R_++R_-) and (c=R_+-R_-), then

\[
 R_+=(f+c)/2,
\]

so the orthogonal fixed-plus-coinvariant sublattice has index two in
(A_1^2).  These are exactly the fixed lattice, rootless
(langle-4\rangle) coinvariant, discriminant action, and primitive glue
of the pinned `0.4.4.1` certificate.  That certificate proves integral
conjugacy of the full Néron--Severi action pairs, not merely equality of
rank and determinant.

The generic cohomological characteristic polynomial is

\[
 (\lambda-1)^{11}(\lambda+1)(\lambda^2+1)^5.
\]

The (H^2)-traces of (sigma,sigma^2,sigma^3) are (10,2,10), so the
fixed-locus Euler characteristics are (12,4,12).  The complete BH row has
the finer fingerprint

\[
\begin{array}{c|c}
\sigma & 8\text{ isolated points}+2\text{ rational curves},\\
\sigma^2 & 5\text{ rational curves}+1\text{ genus-4 curve}.
\end{array}
\]

The second line also follows from the standard 2-elementary formula for
(S(\sigma^2)=NS(X)), whose invariants are ((r,a)=(12,2)).

## Lineage: an inverse-designed specialization

[Artebani--Sarti's order-four classification](https://arxiv.org/abs/1102.4436)
contains the ambient construction.  Their Example 6.3 starts with

\[
 y^2=x^3+a(t)x+b(t),\qquad a(-t)=a(t),\quad b(-t)=-b(t),
\]

and the same automorphism ((x,y,t)\mapsto(-x,iy,-t)).  Their Example 7.2
explains geometrically that two nodes exchanged by the base involution add
an (A_1^2) summand to the invariant lattice of the square.

The present normal form is obtained by designing the remaining root lattice
from the desired Coxeter action:

- force orders ((\ge4,5,10)) at infinity to insert a (II^*/E_8) fiber;
- force a double root of the discriminant in (s=t^2) to insert the
  exchanged (I_2/I_2), hence (A_1^2), pair;
- retain exactly four moduli, the period-domain dimension required for
  generic Picard rank 12.

This is precisely the fiber-root inverse-design criterion in working form:
the algebraic equation is built from the required root summands and diagram
action, and the Picard-rank/moduli calculation proves that no hidden divisor
classes spoil the target lattice.

## Coxeter consequence

The equation supplies the algebraic carrier; the pinned chamber certificate
supplies the Kähler/Weyl transport to Roulleau's complete nef chamber.  The
independent literal Vinberg audit gives:

| invariant | value |
|---|---:|
| hyperbolic dimension | 10 |
| facets | 12 |
| finite vertices | 17 |
| cusps | 2 |
| elliptic edges | 99 |
| failed edge/parabolic audits | 0 |

Therefore the resulting Coxeter polytope is **finite-volume noncompact**.
The presence of two cusps is part of the positive theorem, not evidence for
compactness.

The complete Vinberg output—not only these counts—is pinned by canonical
audit SHA-256
d20267d1d38090f54d650d20d88d71448aa0393686cf64705e9780f91c37bd5a.

## Reproduction

The offline certificate is
[`scripts/certify_ue8a1x2_order4_weierstrass_family.py`](../scripts/certify_ue8a1x2_order4_weierstrass_family.py),
with machine artifact
[`data/enumerations/ue8a1x2_order4_weierstrass_family.json`](../data/enumerations/ue8a1x2_order4_weierstrass_family.json).

Run:

```bash
PYTHONPATH=.:src python3 scripts/certify_ue8a1x2_order4_weierstrass_family.py \
  --check-archive
PYTHONPATH=.:src pytest -q \
  tests/test_ue8a1x2_order4_weierstrass_family.py
```

For a full external replay of Brandhorst--Hofmann row `0.4.4.1`:

```bash
PYTHONPATH=.:src python3 scripts/certify_ue8a1x2_order4_weierstrass_family.py \
  --check-archive --brandhorst-hofmann-zip /path/to/K3Groups.zip
```

The pinned database is version 2.1, DOI `10.5281/zenodo.10410278`.
