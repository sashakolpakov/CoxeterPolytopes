# Finite volume forces a mixed core--tail attachment

**Status (2026-08-10): exact real-hyperbolic lemma.**  This note proves only
a necessary attachment condition.  It does not bound the tail rank or the
number of mixed facets.

Let

\[
                         V=V_0\mathbin\perp E
\]

be a real quadratic space, where \(V_0\) has signature \((1,p)\) and
\(E\) is negative definite of dimension \(m>0\).  Use the hyperboloid
normalization

\[
 \mathbf H(V)=\{x:(x,x)=1,\ x\text{ lies in the chosen positive cone}\}.
\]

## 1. Pure roots give infinite volume

### Lemma 1

Let \(P\subset\mathbf H(V)\) be a full-dimensional polyhedron cut out by
root inequalities.  If every defining root lies entirely in \(V_0\) or
entirely in \(E\), then \(P\) has infinite hyperbolic volume.

### Proof

Away from the measure-zero locus with zero \(E\)-coordinate, every point of
\(\mathbf H(V)\) has a unique expression

\[
                 x=\cosh(t)u+\sinh(t)v,
 \qquad t>0,\quad u\in\mathbf H(V_0),\quad
 v\in S(E):=\{v:(v,v)=-1\}.
\]

The hyperbolic metric and volume form are

\[
 \begin{aligned}
 ds^2&=dt^2+\cosh^2(t)\,ds^2_{\mathbf H(V_0)}
                 +\sinh^2(t)\,ds^2_{S(E)},\\
 d\operatorname {vol}
   &=\cosh^p(t)\sinh^{m-1}(t)\,dt\,
       d\operatorname {vol}_{\mathbf H(V_0)}
       d\operatorname {vol}_{S(E)}.
 \end{aligned}
\]

An inequality from a root in \(V_0\) depends only on \(u\), and an
inequality from a root in \(E\) depends only on \(v\).  Thus, up to the
measure-zero set \(t=0\),

\[
                         P=(0,\infty)\times P_0\times Q,
\]

where \(P_0\subset\mathbf H(V_0)\) and \(Q\subset S(E)\) have positive
measure because \(P\) is full dimensional.  Therefore

\[
 \operatorname {vol}(P)
 =\operatorname {vol}(P_0)\operatorname {vol}(Q)
   \int_0^\infty\cosh^p(t)\sinh^{m-1}(t)\,dt=\infty.
\]

The integral also diverges in the endpoint case \(p=0,m=1\), when its
integrand is \(1\).  Hence a finite-volume full Coxeter chamber for a
nonzero tail must have at least one genuinely mixed core--tail facet. \(\square\)

Equivalently, if all generating reflections are pure, the reflection group
is the product of its core group with a finite negative-definite reflection
group, and its hyperbolic quotient has the same divergent radial factor.

## 2. The exact projection-span consequence

Let \(\Delta_T\) be the pure-tail defining roots and put

\[
                   W_T=\langle s_\alpha:\alpha\in\Delta_T\rangle.
\]

This is a finite subgroup of \(O(E)\) in the integral crystallographic
setting.  Let \(\Delta_M\) be the mixed defining roots, and define

\[
 U=\operatorname {span}_{\mathbf R}
       \{w(\operatorname {pr}_E r):w\in W_T,\ r\in\Delta_M\}
       \subseteq E.                                                     \tag{2.1}
\]

### Lemma 2

If \(P\) has finite volume, then

\[
                              \boxed{U=E}.                              \tag{2.2}
\]

### Proof

The space \(U\) is \(W_T\)-invariant.  If \(E'=U^\perp\ne0\), every mixed
root lies in \(V_0\oplus U\) and hence its reflection fixes \(E'\)
pointwise.  Every pure-core root also lies in \(V_0\oplus U\).

It remains to locate the pure-tail roots.  Since each of their reflections
preserves \(U\), a rank-one reflection root \(\alpha\in E\) must satisfy
either \(\alpha\in U\) or \(\alpha\in U^\perp=E'\): if some \(u\in U\)
has \((u,\alpha)\ne0\), the formula for \(s_\alpha(u)\) forces
\(\alpha\in U\); otherwise \(\alpha\in E'\).  Consequently every defining
root is pure for the new orthogonal decomposition

\[
                         V=(V_0\oplus U)\mathbin\perp E'.
\]

Lemma 1 would then give infinite volume, a contradiction.  Thus
\(E'=0\), proving (2.2). \(\square\)

Condition (2.2) is the strongest projection statement supplied by this
argument.  It concerns the **\(W_T\)-saturation** of the mixed projections,
not their raw linear span and not their number.  In particular, for the
free exceptional tail

\[
 e_1-e_2,\ldots,e_{m-1}-e_m,e_m,
\]

the pure-tail reflections generate \(W(B_m)\), the full signed-permutation
group.  The real \(W(B_m)\)-module \(E\) is irreducible: the orbit of any
nonzero vector spans all of \(E\).  Hence in this case (2.2) forces only
that at least one mixed facet have nonzero tail projection.  One must not
infer \(m\) independent mixed facets, or even that the raw projections of
the displayed mixed facets span \(E\).

More generally, (2.2) says that the mixed attachments must meet every
\(W_T\)-invariant tail factor.  Converting that necessary coupling into a
tail-rank bound requires additional geometric input; finite volume alone
does not provide it.
