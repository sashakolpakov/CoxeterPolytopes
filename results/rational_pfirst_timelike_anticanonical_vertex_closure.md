# Timelike anticanonical vertices have exact flex-cluster carriers

**Status (2026-08-11): proved classification-free exact carrier theorem for
ordinary physical chambers on \(I_{1,n}\), \(3\le n\le8\).  Twice-\((-4)\)
walls and chambers with \(-K\notin P\) remain OPEN.**

This is the nearest non-cusp continuation of the rank-ten anticanonical
cusp theorem.  Here

\[
L=I_{1,n},\qquad K=-3H+E_1+\cdots+E_n,\qquad
F=-K,
\]

with \(3\le n\le8\).  Thus

\[
                         F^2=K^2=9-n>0.                \tag{0.1}
\]

The anticanonical class is timelike, not an ideal ray.

The result is intrinsic and P-first: relative reflection maximality and
ordinary physical canonical colors force a unique chamber up to the
stabilizer of \(K\).  A direct truncation of the flex cubic base cluster
then constructs a rational surface, independently proves \(F\) nef, and
returns the exact full nef cone.

No reflection or weak-del-Pezzo classification is used.

## 1. Surface-free theorem

### Theorem 1.1 (timelike anticanonical-vertex closure)

Fix \(3\le n\le8\), \(L=I_{1,n}\), and the characteristic vector \(K\)
above.  Let \(P\) be a finite-volume integral Coxeter chamber in the chosen
positive component.  Assume:

1. every oriented facet is ordinarily physical:
   \[
   C^2=-1,\ K.C=-1,\qquad\text{or}\qquad
   C^2=-2,\ K.C=0;                                     \tag{1.1}
   \]
2. \(W(P)\) has no proper reflection overgroup inside \(O^+(L)\); and
3. \(K\notin P\).

Then there is a smooth projective rational surface \(S/\mathbf C\) and a
marking

\[
                       \operatorname{NS}(S)\cong L
\]

which carries \(K_S\) to \(K\), realizes every facet by a distinct smooth
rational curve, and satisfies

\[
                         P=\mathbf P(\operatorname{Nef}S). \tag{1.2}
\]

Up to an integral reflection isometry fixing \(K\), the facets are

\[
 E_n,\qquad E_i-E_{i+1}\ (1\le i<n),\qquad
 H-E_1-E_2-E_3.                                        \tag{1.3}
\]

There is one \((-1)\) facet and \(n\) \((-2)\) facets.  The vector
\(F=-K\) is their finite timelike vertex.

The proof occupies Sections 2--5.

## 2. The standard maximal simplex

Let \(P_n\) be the chamber cut out by (1.3).  Its canonical pairings are

\[
(-1,0,\ldots,0),
\]

and the roots form a determinant-\(-1\) basis of \(L\).  The square-\(-2\)
Gram determinant has absolute value \(9-n\), so it is negative definite.

The Coxeter diagram consists of:

- an order-\(4\) edge from \(E_n\) to \(E_{n-1}-E_n\);
- for \(n=3\) only, a second order-\(4\) edge from \(E_3\) to the strict
  flex-line root;
- the coordinate-difference chain;
- the root \(H-E_1-E_2-E_3\), attached to
  \(E_3-E_4\) when \(n\ge4\), and orthogonal to the difference chain when
  \(n=3\).

Exact dualization gives only finite or ideal vertices.  In particular:

\[
\begin{aligned}
v_0&=F=(3,-1,\ldots,-1),&v_0^2&=9-n,\\
v_1&=H-E_1,&v_1^2&=0.
\end{aligned}                                           \tag{2.1}
\]

All other dual rays have positive square.  Hence \(P_n\) is a finite-volume
Coxeter simplex.  The anticanonical vector satisfies

\[
F.E_n=1,\qquad
F.(E_i-E_{i+1})=F.(H-E_1-E_2-E_3)=0,                  \tag{2.2}
\]

so it is a finite vertex, not a cusp.

### Proposition 2.1 (all host reflections)

The side reflections of \(P_n\) generate every integral hyperbolic
reflection of \(I_{1,n}\).  Therefore \(P_n\) is relatively
reflection-maximal inside \(O^+(L)\).

### Proof

For a primitive vector \(u\), unimodularity makes its divisibility one.
If its reflection is integral, then \(u^2\mid2\); hence

\[
                              u^2=-1\quad\text{or}\quad-2. \tag{2.3}
\]

The reflections in \(E_n\) and the coordinate differences generate all
signed permutations of the \(E_i\).  Up to these operations and overall
sign, write

\[
u=(a,-m_1,\ldots,-m_n),\qquad
a\ge0,\quad m_1\ge\cdots\ge m_n\ge0.                   \tag{2.4}
\]

Let \(q=-u^2\in\{1,2\}\) and
\(s=m_1+m_2+m_3\).  If \(s\le a\), then

\[
\sum_{i=1}^n m_i^2
\le m_1^2+m_2^2+(n-2)m_3^2
\le s^2
\le a^2,                                               \tag{2.5}
\]

contradicting \(\sum m_i^2=a^2+q\).  To see the middle inequality, put

\[
m_1=k+A,\qquad m_2=k+B,\qquad m_3=k.
\]

The difference between its right and left sides is

\[
          (9-n)k^2+4k(A+B)+2AB\ge0.                    \tag{2.6}
\]

Thus \(s>a\).  Each \(m_i\le a\) when \(a>0\), so \(s\le3a\).  Equality
forces the simple root \(H-E_1-E_2-E_3\).  Otherwise its reflection changes
the time coefficient to \(2a-s\), whose absolute value is strictly smaller
than \(a\).  Reorder, change overall sign if needed, and repeat.

At \(a=0\), the norm equation leaves one unit exceptional coordinate for
square \(-1\), or two for square \(-2\).  Signed permutations carry these
to \(E_n\) and a coordinate difference.  The descent therefore generates
every integral reflection.  \(\square\)

The doubled weak-log host has the same underlying group of integral
isometries, so this is also the adapted-host maximality needed for exact
nef return.

## 3. Physical coloring rigidifies the canonical marking

The ordinary colors give the homogeneous wall inequalities

\[
F.C=1\quad\text{on a \((-1)\) facet},\qquad
F.C=0\quad\text{on a \((-2)\) facet}.                  \tag{3.1}
\]

This does **not** by itself choose the positive sheet, so we do not yet
assert \(F\in P\).

Relative maximality implies that \(W(P)\) contains every integral
reflection, since adjoining any one of them stays discrete in \(O^+(L)\).
By Proposition 2.1, an oriented reflection isometry \(w\) carries \(P\) to
\(P_n\).  Put

\[
                              G=wF.
\]

Because \(w\) carries the oriented facet halfspaces to the standard
oriented halfspaces, (3.1) says that \(G\) belongs to the full homogeneous
inequality cone of \(P_n\).  It is integral, characteristic, and
\(G^2=9-n\).

The remaining issue is whether the translating reflection changes \(K\).
The following elementary lemma rules that out.

### Lemma 3.1 (unique reduced characteristic vector)

The only characteristic vector \(G\) in the homogeneous standard
inequality cone satisfying

\[
                              G^2=9-n                  \tag{3.2}
\]

is

\[
                              G=F=(3,-1,\ldots,-1).
\]

### Proof

Write a characteristic vector in the standard chamber as

\[
G=(a,-m_1,\ldots,-m_n),\qquad
a\ge m_1+m_2+m_3,\quad
m_1\ge\cdots\ge m_n\ge0.
\]

Every coordinate is odd, so put

\[
a=2\ell+3,\qquad m_i=2k_i+1.
\]

The chamber inequality gives

\[
                    \ell\ge k_1+k_2+k_3=:s,            \tag{3.3}
\]

while (3.2) becomes

\[
             \ell(\ell+3)=\sum_{i=1}^n k_i(k_i+1).     \tag{3.4}
\]

With \(k_1=k+A,\ k_2=k+B,\ k_3=k\), the upper bound obtained from the
first three \(k_i\)'s is

\[
\begin{aligned}
s(s+3)
&-\left[k_1(k_1+1)+k_2(k_2+1)
       +(n-2)k_3(k_3+1)\right]\\
&=(9-n)k(k+1)+4k(A+B)+2AB+2A+2B.                       \tag{3.5}
\end{aligned}
\]

For \(n\le8\), this is strictly positive unless
\(k=A=B=0\).  The remaining \(k_i\) are at most \(k_3\), and
\(\ell(\ell+3)\ge s(s+3)\).  Equality in (3.4) therefore forces

\[
\ell=k_1=\cdots=k_n=0,
\]

which gives \(G=F\).  \(\square\)

Apply Lemma 3.1 to \(G=wF\).  It gives \(G=F\), so \(w\) fixes \(F\) and
\(K=-F\).  Only now, because the standard \(F\) is the finite vertex of
\(P_n\), do we conclude that the original \(F\) lies in the chosen positive
component and is the finite vertex of \(P\).  This proves the surface-free
rigidity statement and the facet list (1.3) without assuming the desired
sheet orientation.

The hypothesis \(K\notin P\) is consistent and nonvacuous: the unique
\((-1)\) facet has \(K.E_n=-1\).

## 4. Explicit carrier from the flex cluster

Consider the plane cubic pencil

\[
                    Y^2Z=X^3+XZ^2+\lambda Z^3.         \tag{4.1}
\]

At its flex \(O=[0:1:0]\), with local coordinates \(x=X/Y,\ z=Z/Y\), the
base ideal is

\[
                         (z^3,\ z-x^3-xz^2).            \tag{4.2}
\]

The complete base cluster has length nine.  Its first three points follow
the flex line.  After those blow-ups, (4.2) becomes, up to local units and
a regular coordinate change,

\[
                              (x^6,w).
\]

The remaining six points form a successive exceptional chain.

For \(3\le n\le8\), blow up only the first \(n\) points and call the surface
\(S_n\).  The transformed pencil has class

\[
                          3H-\sum_{i=1}^nE_i=-K_{S_n}. \tag{4.3}
\]

It has residual base ideal

\[
                              (x^{9-n},w),              \tag{4.4}
\]

but no fixed curve.  A divisor class carrying a pencil without a fixed
curve is nef: for any irreducible curve, choose a pencil member which does
not contain it and use nonnegative local intersection.  Thus

\[
                             -K_{S_n}\ \text{is nef}.   \tag{4.5}
\]

This is the independent geometric nef certificate; no proposed-wall
inequality is used.

The chronological blow-up geometry makes the following distinct smooth
rational curves effective:

\[
E_n,\qquad E_i-E_{i+1}\ (1\le i<n),\qquad
H-E_1-E_2-E_3.                                         \tag{4.6}
\]

The first is the last exceptional curve.  The differences are the strict
transforms of the preceding exceptional curves, and the final class is the
strict flex line.

## 5. Exact full-nef return

Let \(C\) be any irreducible negative curve on \(S_n\).  Since
\(F=-K_{S_n}\) is nef, adjunction gives

\[
                         C^2=2p_a(C)-2+F.C.             \tag{5.1}
\]

The right side shows that \(p_a(C)=0\) and leaves only:

\[
\begin{array}{c|cc}
C^2&F.C&K.C\\ \hline
-2&0&0\\
-1&1&-1.
\end{array}                                             \tag{5.2}
\]

Every possible hidden negative curve therefore has an integral reflection
in \(L\).  Proposition 2.1 and relative maximality put that reflection in
the side group of \(P_n\).  Its mirror cannot cross the fundamental
chamber interior.  An ample class has positive intersection with every
effective curve, fixing the correct side.  Curves of nonnegative square
are handled by the closed future-cone Lorentz inequality.

It follows that every class in the cone over \(P_n\) is nef.  Conversely a
nef class lies in the positive closure and is nonnegative on the effective
facet curves (4.6).  Therefore

\[
                      P_n=\mathbf P(\operatorname{Nef}S_n). \tag{5.3}
\]

Transporting the marking by the \(K\)-fixing reflection from Section 3
realizes the original chamber \(P\).  This completes Theorem 1.1.

## 6. What this closes, and what remains open

For ranks \(4,\ldots,9\), the ordinary timelike-anticanonical lane is now
closed:

\[
\boxed{
\begin{gathered}
\text{relative reflection maximality}
+\text{physical \((-1)/(-2)\) colors}\\
\Longrightarrow
\text{canonical rigidity + explicit rational carrier + exact nef return}.
\end{gathered}}
\]

These are genuine non-cusp examples:

\[
K\notin P,\qquad -K\in P,\qquad (-K)^2>0.
\]

The theorem does not cover:

- twice-primitive physical \((-4)\) facets, for which
  \(F.B=-2\) and hence \(-K\notin P\);
- maximal chambers with canonical escape but \(-K\) outside the chamber;
- rank ten inside this timelike argument; it is handled by the separate
  characteristic-cusp theorem;
- higher ranks inside this descent argument; the combined ordinary theorem
  now excludes them directly by the spacelike anticanonical wall-cone
  escape; or
- nonunimodular adapted hosts.

Thus the genuinely open property-4 lane left after the combined theorem is
the twice-\((-4)\), exponent-two adapted-host boundary (including its
nonunimodular hosts and effectivity gate).

The exact rank-by-rank replay is
[scripts/certify_rational_pfirst_timelike_anticanonical_vertex_closure.py](../scripts/certify_rational_pfirst_timelike_anticanonical_vertex_closure.py).
