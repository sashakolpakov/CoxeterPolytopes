# Independent hostile audit: compact \(\mathbf H^3\) general-fork field finiteness

**Verdict (2026-08-12): the repaired theorem passes.**  The first owner
draft did not pass: its invariant formula had the wrong sign, its norm
estimate omitted the selected-place correction, it sent the label-three
fork through a tail chart to which it did not belong, and it inferred a
nonzero determinant from linear independence alone.  The terminal owner
version repairs those points.  This audit independently reconstructs the
invariant arithmetic, the uniform label cutoff, the Fekete step, and the
global compact-three-dimensional case split.

The theorem proved is field finiteness for compact arithmetic Coxeter
chambers in \(\mathbf H^3\), in the Vinberg quadratic-form/cycle-field
scope.  It does not prove chamber finiteness, a dimension bound, or a facet
bound.  It does not invoke Nikulin's global ground-field theorem, Nikulin's
fixed-field/fixed-rank reflective-lattice theorem, ABSW, or a classification
of compact Coxeter polyhedra.

## 1. The arithmetic coordinate survives root rescaling

For a root line represented by \(r_i\), put

\[
 c_{ij}=\frac{2q(r_i,r_j)}{q(r_i,r_i)},\qquad
 t_{ij}=c_{ij}c_{ji},\qquad
 z=c_{12}c_{23}c_{31}.
\]

Under \(r_i\mapsto\lambda_i r_i\), the directed entry changes by
\(c_{ij}\mapsto\lambda_j\lambda_i^{-1}c_{ij}\).  Thus \(t_{ij}\) and
\(z\) are invariant closed Cartan products.  If \(s_i\) is the integral
reflection and \(T_i=1-s_i\), then a closed product is

\[
 \operatorname {Tr}(T_{i_1}\cdots T_{i_r}).
\]

It is therefore integral even when the stable \(\mathcal O_k\)-module is
projective, no individual Cartan entry is integral, and a square-\((-2)\)
normalization is unavailable over \(k\).

For the determinant calculation only, symmetrize at a real place:

\[
 G=\begin{pmatrix}-2&u&a\\u&-2&b\\a&b&-2\end{pmatrix},
 \qquad
 \frac{\det G}{2}=u^2+a^2+b^2+uab-4.                    \tag{1.1}
\]

There is a sign trap here.  With the displayed definition of \(c_{ij}\)
and diagonal entries \(-2\), one has \(c_{ij}=-G_{ij}\).  Consequently
\(z=-uab\), and the invariant formula is

\[
 \boxed{\beta=t_{12}+t_{13}+t_{23}-z-4\in\mathcal O_k}. \tag{1.2}
\]

The unrepaired formula with \(+z\) is false in this convention.  Formula
(1.2) proves that \(\beta\), rather than an auxiliary square root \(u\), is
the legitimate integral function on the diagonal-torus quotient of the
marked quadric arrangement.

## 2. Sign separation and the whole field

At the selected place, walls 1 and 2 are ultraparallel, so their span
\(U\) is nondegenerate of signature \((1,1)\).  The ambient space has
signature \((1,3)\), hence \(U^\perp\) is negative definite.  Wall 3 is
independent of \(U\): if its normal \(e\) lay in the span of the other two
normals, then the two sections \(H_e\cap H_{r_1}\) and
\(H_e\cap H_{r_2}\) would be the same hyperplane in \(H_e\), contrary to
their being distinct sides of the facet.  Its
\(U^\perp\)-component is therefore nonzero.  Thus the three-wall packet has
signature \((1,2)\), not merely rank three, and

\[
                         0<\beta_0<256.                  \tag{2.1}
\]

The upper bound uses \(2<u<14\) and \(a,b<2\).  At every nonselected
embedding the three independent conjugate roots lie in a definite space,
so their Gram matrix is negative definite and

\[
                         \beta^\tau<0.                   \tag{2.2}
\]

If \(F=\mathbf Q(\beta)\) were a proper subfield of \(k\), the selected
embedding of \(F\) would have a second real extension to the totally real
separable extension \(k/F\).  That extension fixes the positive value
\(\beta_0\), contradicting (2.2).  Hence

\[
                         \boxed{k=\mathbf Q(\beta)}.      \tag{2.3}
\]

This is a single-coordinate full-field statement; no unproved integrality
of \(u,a,b\) is being used.

## 3. The exact conjugate interval

Let

\[
 s_m=4\sin^2(\pi/m),\quad s_n=4\sin^2(\pi/n),\quad
 E=\mathbf Q(s_m,s_n)\subset k.
\]

The inclusion follows from \(s_m=4-t_{13}\) and
\(s_n=4-t_{23}\).  As a quadratic in \(u\), (1.1) has discriminant

\[
 (4-a^2)(4-b^2)=s_ms_n
\]

and minimum \(-s_ms_n/4\).  Therefore every nonselected extension of a
base embedding \(\rho:E\hookrightarrow\mathbf R\) satisfies

\[
 -\frac{s_m^\rho s_n^\rho}{4}\ \le\ \beta^\tau<0.       \tag{3.1}
\]

The lower inequality is not strict: equality at the quadratic minimum is
compatible with negative definiteness.  The closed interval (3.1) is the
one required by Fekete and changes none of the estimates.

For fixed \((m,n)\), the product of interval capacities is

\[
 \prod_{\rho:E\hookrightarrow\mathbf R}
 \frac{s_m^\rho s_n^\rho}{16}
 =\frac{|N_E(s_ms_n)|}{16^{[E:\mathbf Q]}}<1,            \tag{3.2}
\]

because every conjugate of either \(s_d\) is strictly between 0 and 4.
The relative Fekete theorem permits one fixed exceptional embedding, here
\(0\le\beta_0\le256\).  It bounds \([E(\beta):E]\).  All conjugates are
then bounded, so the coefficients of the relative minimal polynomial lie
in a bounded Minkowski box in \(\mathcal O_E\).  Hence there are only
finitely many \(\beta\) and fields for fixed labels.

## 4. Correct uniform cutoff for both labels

Write \(D=[E:\mathbf Q]\), \(e=[k:E]\), and
\(L_\rho=s_m^\rho s_n^\rho/4\).  Exactly one of the \(e\) extensions over
the selected base embedding \(\rho_0\) is exceptional.  Consequently the
correct norm estimate is

\[
\begin{aligned}
 1\le |N_{k/\mathbf Q}(\beta)|
 &<256L_{\rho_0}^{e-1}\prod_{\rho\ne\rho_0}L_\rho^e\\
 &=\frac{1024}{s_{m,0}s_{n,0}}
   \left(\frac{|N_E(s_ms_n)|}{4^D}\right)^e.             \tag{4.1}
\end{aligned}
\]

The factor in front was absent from the first draft.  Since
\(s_{d,0}>16/d^2\), and \(M=\max(m,n)\), (4.1) gives

\[
 |N(\beta)|<4M^4R^e,qquad
 R=\frac{|N_E(s_ms_n)|}{4^D}.                             \tag{4.2}
\]

Put \(D_d=[\mathbf Q(s_d):\mathbf Q]=\varphi(d)/2\) and

\[
 \gamma(d)=|N_{\mathbf Q(s_d)/\mathbf Q}(s_d)|
 =\begin{cases}p,&d=p^r,\\1,&\text{otherwise},\end{cases}
 \qquad h(d)=\gamma(d)^{2/\varphi(d)}.
\]

The tower formula gives exactly

\[
                         R=\left(\frac{h(m)h(n)}4\right)^D. \tag{4.3}
\]

One has \(h(d)\le3\) for all \(d\ge3\), with equality only at \(d=3\).
For a prime power this is the elementary inequality
\(p^{2/\varphi(p^r)}\le3\); otherwise \(h(d)=1\).

The bound \(\varphi(d)^2\ge d/2\) implies, for \(M\ge10^6\),

\[
 \log h(M)\le\frac{2\log M}{\sqrt{M/2}}<\log(13/12).     \tag{4.4}
\]

The left side is decreasing after \(e^2\).  At the endpoint,
\(\log10<7/3\), \(\sqrt{500000}>700\), and
\(\log(13/12)>1/13\) give the stronger comparison
\(\log h(M)<1/25<1/13<\log(13/12)\).  Thus

\[
 R<(13/16)^D,qquad D\ge\varphi(M)/2>\sqrt{M/2}/2.        \tag{4.5}
\]

At \(M=10^6\), the exact integer inequality

\[
                 16^{353}>4(10^6)^4 13^{353}             \tag{4.6}
\]

shows \(4M^4(13/16)^D<1\).  The continuous upper envelope is decreasing
thereafter: \(\log(16/13)>3/16\) puts its derivative threshold below
\(512/3<1000=\sqrt{10^6}\).  Equations (4.2)--(4.6) contradict the
nonzero integral norm.  Therefore

\[
                         \boxed{m,n<10^6}.                \tag{4.7}
\]

This direct argument includes \(m=3\) or \(n=3\).  A label-three fork is
not a tail, and no such reduction is needed.

## 5. The global compact-\(\mathbf H^3\) split is exhaustive

Choose a minimal-visual-aperture facet \(E\).  Its closed star has all
normalized pairings below 14.  Cyclically write its side walls as
\(\delta_1,\ldots,\delta_N\).

If every \(\delta_i\perp E\), restriction to \(E^\perp\) is literal:
the same roots bound a compact arithmetic Coxeter polygon, and their
Lorentzian packet inherits the full ambient cycle field.  Fields of
arithmetic Coxeter polygons are finite by the exact genus-zero real-DM
orbifold-curve argument: pass to a maximal arithmetic Fuchsian over-lattice,
whose coarse quotient is still genus zero, and apply
Long--Maclachlan--Reid.  No bounded index is required; the field is a
commensurability invariant.

Otherwise select \(j\) with \(E\cdot\delta_j>0\) and take the edge
\(r=E\cap\delta_{j+1}\).  Compact simplicity makes
\(\delta_j,\delta_{j+2}\) its endpoint-exclusive walls.  Let their
normalized pairing be \(u\).

- If \(u=2\), the independent parabolic pair inherits the full field and
  forces \(k=\mathbf Q\).
- If \(u<2\), an inclusion-minimal positive principal packet is a Coxeter
  triangle, a rank-at-least-four Lanner packet, or contains a proper
  parabolic packet.  In the last case it is a noncompact quasi-Lanner
  simplex, whose \(k\)-rational cusp forces \(k=\mathbf Q\).  The first two
  field sets are finite.
- Otherwise \(2<u<14\).  If \(E\cdot\delta_{j+2}=0\), the three walls form
  the previously audited tail chart.  If it is positive, they form exactly
  the arbitrary-label fork audited in Sections 1--4.

These alternatives use every possible sign of the two relevant incidences
and every possible value of \(u\), so there is no uncovered compact
polyhedron branch.

## 6. Scope of the PASS

The compact-\(\mathbf H^3\) field theorem is a genuine integral-point
pushback on the exact quadric-arrangement chart, with a real algebraic
orbifold-curve branch.  Its terminal inputs are relative Fekete theory and
the Long--Maclachlan--Reid arithmetic genus-zero theorem.  It is therefore
AG-arrangement-mediated arithmetic, not a consequence of BAB or KSBA alone.

The executable audit is
[`audit_compact_h3_general_fork_arrangement_field_finiteness_independent_hostile.py`](../scripts/audit_compact_h3_general_fork_arrangement_field_finiteness_independent_hostile.py),
with independent tests in
[`test_compact_h3_general_fork_arrangement_field_finiteness_independent_hostile.py`](../tests/test_compact_h3_general_fork_arrangement_field_finiteness_independent_hostile.py).

The terminal owner bytes audited here have SHA-256

```text
4a971c64460a7d234114ecb12b793e613b080945850098da37b903d2c50876d5  results/compact_h3_general_fork_arrangement_field_finiteness.md
9426c6a4c1652217a5e7fe15ca2a317cd0503f0bcd5a4cb6405713eed397e219  scripts/certify_h3_general_fork_field_finiteness.py
06255c5e987f8444ecd8a80deb2ef0c9bf6efee67415370d4c73c7037a13ad22  tests/test_h3_general_fork_field_finiteness.py
```

The independent executable bytes have SHA-256

```text
d1d243406fe672105493cf23e2daf0fcd96a7bd0dc0535f431bff4105bee8b31  scripts/audit_compact_h3_general_fork_arrangement_field_finiteness_independent_hostile.py
a6cfe5f124fc2c88b0aa872d60dcd41cd0bf99536cb86fd5a05c3effd10413ed  tests/test_compact_h3_general_fork_arrangement_field_finiteness_independent_hostile.py
```

The combined owner and independent suites pass: **9 tests passed**.
