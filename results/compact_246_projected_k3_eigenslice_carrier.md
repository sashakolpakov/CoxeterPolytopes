# A projected-root K3 eigenslice carrier for the compact `Delta(2,4,6)` chamber

**Status (2026-08-11): proved by an exact finite root-window certificate.**

The compact arithmetic triangle

\[
                         P=\Delta(2,4,6)
\]

has no carrier in the current exact-return three-branch compiler, but it does
have a natural carrier after making precisely the escape allowed in that
no-go theorem: for the rational ground field, replace literal full-nef return
by a proper rational projector eigenslice in a larger K3 Picard lattice.

Concretely, there is an even lattice \(S\) of signature \((1,5)\), a
primitive rank-three sublattice \(T=(\mathbf Z^3,15G)\), and a rational
self-adjoint orthogonal projector \(e:S_{\mathbf Q}\to T_{\mathbf Q}\) such
that a projective complex K3 surface with \(\operatorname {NS}(X)=S\) has

\[
 \boxed{
 \mathbf P\bigl(\operatorname {Nef}(X)_{\mathbf R}
       \cap E_1(e)\bigr)=\Delta(2,4,6).}                 \tag{0.1}
\]

Each of the three restricted facets is the projection of a primitive
square-\((-2)\) class of \(S\), hence becomes a smooth rational curve after
Weyl framing.  An exact positive-definite enumeration proves that no other
ambient square-\((-2)\) mirror meets the closed restricted triangle.

This is a **positive fourth-category calibration**, not a repair of the
current theorem for all inputs.  The carrier has Picard rank six while its
returned eigenslice has rank three.  It does not give literal full-nef return
for \(k=\mathbf Q\), an integral action of the raw Coxeter group on all of
\(S\), a class-wide existence theorem, or finite-host control.  The K3 Picard
ceiling does give the carrier-relative bound \(n\le19\); see Section 5.

## 1. The rank-six Picard lattice

Use the negative wall-normal Gram

\[
 G=\begin{pmatrix}
 -12&6&0\\
 6&-6&3\\
 0&3&-2
 \end{pmatrix}.                                          \tag{1.1}
\]

It has signature \((1,2)\), determinant \(36\), and normalized edge labels
\((4,6,2)\).  Put \(s=15\), let

\[
                    T=(\mathbf Zr_1\oplus\mathbf Zr_2
                         \oplus\mathbf Zr_3,\,sG),       \tag{1.2}
\]

and adjoin basis vectors \(\delta_1,\delta_2,\delta_3\).  In the ordered
basis \((r_1,r_2,r_3,\delta_1,\delta_2,\delta_3)\), define

\[
 A=\begin{pmatrix}sG&G\\G&-2I_3\end{pmatrix}.            \tag{1.3}
\]

This is visibly integral and even.  Its Schur complement over \(sG\) is

\[
 K=-2I_3-{G\over15}
  =\begin{pmatrix}
  -6/5&-2/5&0\\
  -2/5&-8/5&-1/5\\
  0&-1/5&-28/15
  \end{pmatrix}.                                          \tag{1.4}
\]

The leading principal minors of \(-K\) are

\[
                         {6\over5},\quad {44\over25},
                         \quad {1214\over375},            \tag{1.5}
\]

so \(K\) is negative definite and \(S=(\mathbf Z^6,A)\) has signature
\((1,5)\).  Exact Smith reduction gives

\[
 \det A=-393336,
 \qquad \operatorname {SNF}(A)=
       \operatorname {diag}(1,1,1,1,18,21852).            \tag{1.6}
\]

In particular, \(\ell(A_S)=2\).  The first three basis vectors form a direct
summand of the underlying abelian group, so \(T\hookrightarrow S\) is
primitive.

The orthogonal projection is

\[
 e=\begin{pmatrix}I_3&I_3/15\\0&0\end{pmatrix}.           \tag{1.7}
\]

Direct multiplication gives \(e^2=e\) and \(e^tA=Ae\).  Thus its
eigenvalue-one space is exactly \(T_{\mathbf Q}\), while its kernel is the
negative-definite Schur-complement space.  Moreover

\[
                     e(\delta_i)={r_i\over15},
 \qquad \delta_i^2=-2.                                   \tag{1.8}
\]

The three physical roots therefore restrict to the three raw facet lines.

## 2. The complete visible-root reduction

Every \(\rho\in S\) has unique coordinates

\[
              \rho=x+\sum_{i=1}^3a_i\delta_i,
       \qquad x\in\mathbf Z^3,\quad a\in\mathbf Z^3.     \tag{2.1}
\]

Its projection is \(x+a/15\).  Set

\[
                         y=15x+a\in\mathbf Z^3.           \tag{2.2}
\]

The restricted mirror is \(G(-,y)=0\), and \(y\equiv a\pmod {15}\).
Orthogonal decomposition of (1.3) gives the exact root equation

\[
 \rho^2=-2
 \quad\Longleftrightarrow\quad
 G(y,y)=G(a,a)+30\bigl(\lVert a\rVert^2-1\bigr).         \tag{2.3}
\]

If the restricted mirror meets the open hyperbolic positive cone, its normal
has \(G(y,y)<0\): in signature \((1,2)\), the orthogonal of a timelike normal
is negative definite, while a null normal is only tangent to the boundary.
Hence

\[
                         a^t(30I_3+G)a<30.                \tag{2.4}
\]

The matrix in (2.4) is

\[
 M=\begin{pmatrix}18&6&0\\6&24&3\\0&3&28\end{pmatrix}, \tag{2.5}
\]

with leading principal minors \(18,396,10926\).  The dual Cauchy bounds are

\[
  a_1^2<{1105\over607},\qquad
  a_2^2<{840\over607},\qquad
  a_3^2<{660\over607}.                                  \tag{2.6}
\]

Thus every coordinate is in \(\{-1,0,1\}\).  Exhausting those 27 triples
leaves exactly

\[
                   a=0,\quad \pm e_1,\quad\pm e_2,
                         \quad\pm e_3.                    \tag{2.7}
\]

The case \(a=0\) is impossible for a root: then \(y=15x\), and (2.3) would
give \(15G(x,x)=-2\).  A zero projection is impossible as well.  If \(y=0\),
then \(a=15b\), and (2.3), divided by 15, has left side divisible by 15 and
right side 2.

For \(a=\pm e_i\), equation (2.3) becomes

\[
                 G(y,y)=G(e_i,e_i)\in\{-12,-6,-2\},
               \qquad y\equiv\pm e_i\pmod {15}.          \tag{2.8}
\]

It remains only to prove that a solution of (2.8) whose mirror meets `P`
cannot escape far along the indefinite quadric.

## 3. The finite hidden-root window

Take the rational interior point and three vertices

\[
 p=G^{-1}(1,1,1)^t=\left({11\over12},2,{5\over2}\right),
 \qquad v_i=G^{-1}e_i.                                   \tag{3.1}
\]

Then \(G(p,r_i)=1\), \(p^2=65/12\), and each \(v_i\) is the finite vertex
opposite the \(i\)-th wall.  The exact values of \(\cosh(2d(p,v_i))\) are

\[
                         {177\over65},\qquad {79\over65},
                         \qquad {17\over13}.             \tag{3.2}
\]

Hyperbolic distance from \(p\) is convex on geodesic segments, so every
point of `P` lies within the vertex radius \(R\) defined by

\[
                         \cosh(2R)={177\over65}.          \tag{3.3}
\]

Define the positive-definite ellipsoid

\[
 Q_p(y)={2G(p,y)^2\over G(p,p)}-G(y,y).                  \tag{3.4}
\]

Its matrix and leading principal minors are

\[
 Q_p={1\over65}
 \begin{pmatrix}
 804&-366&24\\
 -366&414&-171\\
 24&-171&154
 \end{pmatrix},qquad
 {804\over65},\ {612\over13},\ 36.                    \tag{3.5}
\]

If \(G(y,y)=-d\) and \(y^\perp\) meets `P`, the distance formula for a
spacelike hyperplane gives

\[
                         Q_p(y)\le d\cosh(2R).            \tag{3.6}
\]

Here \(d\le12\), so \(Q_p(y)\le2124/65\).  Applying Cauchy to the three
coordinate functionals using \(Q_p^{-1}\) gives

\[
 y_1^2\le {31329\over4225},\qquad
 y_2^2\le {111864\over4225},\qquad
 y_3^2\le {36108\over845}.                              \tag{3.7}
\]

Therefore

\[
                           |y_1|\le2,\qquad |y_2|\le5,
                           \qquad |y_3|\le6.              \tag{3.8}
\]

This box contains only \(5\cdot11\cdot13=715\) lattice points.  Exact
enumeration against the sharper per-norm inequality (3.6) leaves 22
candidate normals.  Imposing (2.8) modulo 15 leaves exactly

\[
                              y=\pm e_1,\quad\pm e_2,
                              \quad\pm e_3.               \tag{3.9}
\]

Thus every ambient square-\((-2)\) mirror meeting the **closed** restricted
triangle is one of its three original facet mirrors.  Conversely the roots
\(\delta_i\) expose all three.  The component of the restricted root
arrangement containing \(p\) is exactly `P`.

There is no ambiguity upstairs over a facet.  For restricted normal
\(y=e_i\), the ellipsoid and congruence calculation forces \(a=e_i\), hence
\(x=0\) and \(\rho=\delta_i\); the negative sign gives only
\(-\delta_i\).  A generic point of that facet therefore lies on the unique
ambient root hyperplane \(\delta_i^\perp\), and the segment from \(p\) to it
meets no earlier root hyperplane.  Thus \(\delta_i^\perp\) is an actual wall
of the ambient chamber containing \(p\), not merely a root mirror whose
projection happens to have the right line.

## 4. K3 realization and rational Coxeter correspondences

The K3 lattice has signature \((3,19)\).  A complement of \(S\) has rank 16,
signature \((2,14)\), and

\[
                         16\ge\ell(A_S)+2=4.              \tag{4.1}
\]

The stable primitive-embedding criterion therefore embeds \(S\) primitively
in the K3 lattice.  A generic period in its orthogonal has exact Picard
lattice \(S\).  There is no root with zero projection, and Section 3 shows
that \(p\) is regular for the ambient root arrangement.  Weyl framing sends
the ambient root chamber containing \(p\) to the ample chamber.  The three
\(\delta_i\) are then effective simple roots.  They are irreducible: if one
decomposed effectively, at a generic point of its nef wall every negative
component would supply another ambient root wall, contradicting the
uniqueness proved at the end of Section 3.  Adjunction therefore makes them
smooth rational \((-2)\)-curves, proving (0.1).

The raw reflection matrices preserve \(G\) integrally on \(T\).  Extending a
raw reflection \(s_i\) by the identity on \(T_{\mathbf Q}^{\perp}\) gives,
in the basis (1.3),

\[
       \widetilde s_i=
       \begin{pmatrix}s_i&(s_i-I_3)/15\\0&I_3\end{pmatrix}. \tag{4.2}
\]

Each is a rational involutive isometry of \(S_{\mathbf Q}\), commutes with
the eigenframe, and is self-adjoint for the intersection form.  The projector
and these operators are algebraic on \(\operatorname {NS}_{\mathbf Q}\): a
rational endomorphism of the divisor space is represented by a rational
linear combination of divisor-product correspondences using a dual divisor
basis.

This last statement is deliberately rational.  The matrices (4.2) need not
preserve \(S\) integrally, and the physical K3 reflections in the
\(\delta_i\) do not equal the raw Coxeter reflections on the eigenslice.
Those are real constraints on any attempt to promote this single carrier to
a new universal category.

## 5. Exact conclusion

The construction proves the nonempty escape clause left by the
[`Delta(2,4,6)` exact-compiler counterexample](compact_246_raw_maximal_property4_counterexample.md):

\[
 \boxed{
 \begin{gathered}
 \operatorname {NS}(X)=S,\quad \rho(X)=6,\quad
 T=(\mathbf Z^3,15G)\subset S\text{ primitive},\\
 \delta_i^2=-2,\quad e(\delta_i)=r_i/15,\quad
 \mathbf P(\operatorname {Nef}(X)\cap E_1(e))=\Delta(2,4,6).
 \end{gathered}}                                         \tag{5.1}
\]

So `Delta(2,4,6)` does not obstruct algebraic geometry itself.  It separates
two contracts sharply:

- the current rational-field **literal full-nef** Property-4 compiler is not
  exhaustive; but
- a **proper rational projector/eigen-nef** K3 category already recovers this
  counterexample with exact physical projected walls and no hidden cuts.

This fourth species still has an immediate carrier-relative AG ceiling.  If
an \(n\)-dimensional chamber is returned on the eigenvalue-one space of a
complex K3 Picard lattice, then

\[
                         n+1=\dim E_1(e)
                            \le\rho(X)\le20,
 \qquad\text{hence}\qquad n\le19.                        \tag{5.2}
\]

Thus no new dimension theorem is needed for the explicitly K3 version of
the category.  What remains open is class-wide nonemptiness and finite-host
control.  An arbitrary non-K3 rational eigenslice category would, of course,
need its own boundedness theorem.

## 6. Replay

```bash
python3 -m py_compile \
  scripts/certify_compact_246_projected_k3_eigenslice_carrier.py

PYTHONPATH=src:. python3 -m \
  scripts.certify_compact_246_projected_k3_eigenslice_carrier --check-archive

PYTHONPATH=src:. pytest -q \
  tests/test_compact_246_projected_k3_eigenslice_carrier.py
```

Frozen certificate:
[`compact_246_projected_k3_eigenslice_carrier.json`](../data/enumerations/compact_246_projected_k3_eigenslice_carrier.json).
