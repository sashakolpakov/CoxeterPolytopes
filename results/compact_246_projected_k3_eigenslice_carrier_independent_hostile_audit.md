# Independent hostile audit: the projected (2,4,6) K3 eigenslice carrier

**Verdict (2026-08-11): PASS at the stated projected-root scope.** The
rank-six lattice is an even hyperbolic K3 Picard lattice candidate, the
rank-three trace block is primitive, the rational projector is orthogonal
and self-adjoint, and the finite root calculation is exhaustive. Exactly
three ambient root hyperplanes meet the closed restricted triangle, with
unique roots \(\pm\delta_i\) above them. After primitive K3 embedding,
generic period choice, and Weyl framing, these are actual irreducible
smooth rational facet curves and

\[
\mathbf P\!\left(\operatorname{Nef}(X)_{\mathbf R}\cap E_1(e)\right)
=\Delta(2,4,6).
\]

This is one positive fourth-category carrier. It lies outside the original
\(k=\mathbf Q\) Property-2 requirement of literal full-nef return. It does
not prove a universal fourth branch or finite-host theorem, and its raw
Coxeter reflections extend only rationally, not integrally, to the full
Picard lattice.

The audit imports neither the owner certifier nor SymPy. It recomputes every
matrix over Python Fractions, obtains the Smith invariants from gcds of all
square minors, and rebuilds both finite ellipsoid enumerations. No owner byte
was edited.

## 1. Frozen owner packet

~~~text
f4c4a038610e7c47600af971154ca056d7ae143546812582619e83ea7286cd1b  results/compact_246_projected_k3_eigenslice_carrier.md
346ed6ed0395008e0a98888efd008a4d74be37e6551a4b2e377f21f42a6aac4a  scripts/certify_compact_246_projected_k3_eigenslice_carrier.py
6d94dbab080c4450bbe17bc7d3da091c550a0ac61965ee926de0eada42f537e1  data/enumerations/compact_246_projected_k3_eigenslice_carrier.json
315e144769eb8c86e1ca12551210e25a402c0c2b0ce863a54b40c3cfcdb8a485  tests/test_compact_246_projected_k3_eigenslice_carrier.py
~~~

## 2. Ambient lattice, primitive trace block, and projector

With

\[
G=\begin{pmatrix}-12&6&0\\6&-6&3\\0&3&-2\end{pmatrix},
\qquad
A=\begin{pmatrix}15G&G\\G&-2I_3\end{pmatrix},
\]

the exact LDL pivots of \(G\) are \((-12,-3,1)\). The Schur complement of
\(15G\) in \(A\) is

\[
K=-2I_3-\frac1{15}G,
\]

and the leading minors of \(-K\) are

\[
\frac65,\qquad\frac{44}{25},\qquad\frac{1214}{375}.
\]

Hence \(A\) has signature \((1,5)\). It is integral and even. Independent
minor-gcd reduction gives determinantal divisors

\[
1,1,1,1,18,393336
\]

and therefore

\[
\det A=-393336,\qquad
\operatorname{SNF}(A)=(1,1,1,1,18,21852),\qquad
\ell(A_S)=2.
\]

The first three coordinate vectors form a direct summand of
\(\mathbf Z^6\), so

\[
T=(\mathbf Z^3,15G)\hookrightarrow S=(\mathbf Z^6,A)
\]

is primitive. The matrix

\[
e=\begin{pmatrix}I_3&I_3/15\\0&0\end{pmatrix}
\]

satisfies \(e^2=e\) and \(e^tA=Ae\). Its eigenvalue-one space is precisely
\(T_{\mathbf Q}\), while its kernel has Gram \(K\) and is negative
definite. For the last three basis vectors,

\[
\delta_i^2=-2,\qquad e(\delta_i)=r_i/15.
\]

Thus the three proposed physical roots project to the raw wall lines.

## 3. Why every visible root enters a seven-point \(a\)-shell

Write an arbitrary ambient root as

\[
\rho=x+\sum_i a_i\delta_i,\qquad x,a\in\mathbf Z^3,
\qquad y=15x+a.
\]

Completing the square in \(A\), independently checked coefficientwise,
gives

\[
15\rho^2=G(y,y)-G(a,a)-30\lVert a\rVert^2.
\]

Consequently \(\rho^2=-2\) is equivalent to

\[
G(y,y)=G(a,a)+30(\lVert a\rVert^2-1).                 \tag{3.1}
\]

The sign step is essential. If \(G(y,y)>0\), then \(y^\perp\) is negative
definite and misses the hyperbolic positive cone. If \(G(y,y)=0\), its
hyperplane is tangent only at the ideal boundary. Every point of the closed
target triangle is an interior positive vector because the triangle is
compact. Therefore a nonzero restricted root mirror meeting the closed
triangle must satisfy

\[
G(y,y)<0.
\]

Combining this with (3.1) yields the positive-definite ellipsoid

\[
a^t(30I_3+G)a<30,\qquad
30I_3+G=
\begin{pmatrix}18&6&0\\6&24&3\\0&3&28\end{pmatrix}.    \tag{3.2}
\]

Its leading minors are \(18,396,10926\). Dual Cauchy bounds give

\[
a_1^2<\frac{1105}{607},\qquad
a_2^2<\frac{840}{607},\qquad
a_3^2<\frac{660}{607}.
\]

Thus every \(a_i\in\{-1,0,1\}\). Exhausting all \(27\) triples in (3.2)
leaves exactly

\[
a=0,\quad \pm e_1,\quad\pm e_2,\quad\pm e_3.
\]

The zero case cannot be a root: \(y=15x\) would make (3.1) imply
\(15G(x,x)=-2\). Nor can a root have zero projection: \(y=0\) forces
\(a=15b\), after which (3.1) equates a multiple of \(15\) with \(2\).
Every visible root therefore has \(a=\pm e_i\), norm
\[
G(y,y)\in\{-12,-6,-2\},
\qquad y\equiv a\pmod{15}.                              \tag{3.3}
\]

## 4. The \(Q_p\) window and exhaustive residue filter

Take

\[
p=G^{-1}(1,1,1)^t=\left(\frac{11}{12},2,\frac52\right),
\qquad p^2=\frac{65}{12}.
\]

Its wall slacks are all one. The three finite vertices \(v_i=G^{-1}e_i\)
give

\[
\cosh(2d(p,v_i))=
\frac{177}{65},\quad\frac{79}{65},\quad\frac{17}{13}.
\]

Hyperbolic distance from \(p\) is convex on geodesic segments, and the
triangle is the geodesic convex hull of its vertices. Hence every point of
\(P\) lies within the radius \(R\) satisfying

\[
\cosh(2R)=\frac{177}{65}.                               \tag{4.1}
\]

Define

\[
Q_p(y)=\frac{2G(p,y)^2}{G(p,p)}-G(y,y)
=\frac1{65}y^t
\begin{pmatrix}
804&-366&24\\
-366&414&-171\\
24&-171&154
\end{pmatrix}y.
\]

Its leading minors are \(804/65,612/13,36\), so it is positive definite.
For \(G(y,y)=-d\), the exact hyperplane-distance identity is

\[
Q_p(y)=d\cosh\!\left(2\,d(p,y^\perp)\right).            \tag{4.2}
\]

If \(y^\perp\) meets \(P\), equations (4.1)--(4.2), together with
\(d\le12\), give

\[
Q_p(y)\le\frac{2124}{65}.
\]

Dual Cauchy bounds are

\[
y_1^2\le\frac{31329}{4225},\qquad
y_2^2\le\frac{111864}{4225},\qquad
y_3^2\le\frac{36108}{845}.
\]

Therefore

\[
|y_1|\le2,\qquad |y_2|\le5,\qquad |y_3|\le6.            \tag{4.3}
\]

The independent replay checks all \(715\) points in this box against the
sharper norm-dependent version of (4.2). Exactly \(22\) remain. Applying
the three norms and six congruences in (3.3) leaves only

\[
(a,y)=(\pm e_i,\pm e_i),\qquad i=1,2,3.
\]

Indeed, all coordinate bounds in (4.3) are below the modulus \(15\), so no
other representative of any residue class is omitted. Then \(x=(y-a)/15=0\):
the unique ambient roots over the three facet lines are
\(\pm\delta_1,\pm\delta_2,\pm\delta_3\).

This also proves that the \(\delta_i^\perp\) are actual walls of the ambient
root chamber containing \(p\). At a generic point of a raw facet, only the
corresponding \(\delta_i\)-hyperplane occurs upstairs, while the segment
from \(p\) to that point has no earlier root crossing.

## 5. Primitive K3 realization and exact eigen-nef equality

The complement of \(S\) in the K3 lattice has desired signature \((2,14)\)
and rank \(16\). Since

\[
16\ge\ell(A_S)+2=4,
\]

the stable primitive-embedding theorem supplies a primitive embedding
\(S\hookrightarrow\Lambda_{\mathrm{K3}}\). A generic period in
\(S^\perp\) has exact Picard lattice \(S\). The vector

\[
12p=(11,24,30)\in T
\]

is integral, positive, and regular for the ambient root arrangement.
Weyl framing sends its ambient root chamber to the ample chamber.

The roots \(\delta_i\) are then effective simple roots. They are
irreducible, not merely effective classes: at a generic point of their nef
wall, every component of an effective decomposition is orthogonal to a
positive class and hence is itself a \((-2)\)-root. The unique-hyperplane
calculation above makes every such component proportional to \(\delta_i\);
primitivity forces one component with multiplicity one. Adjunction makes
each curve smooth rational.

There are two inclusions in the chamber equality. The connected interior of
\(P\) contains \(p\) and has no ambient root crossing, so it lies in the
ambient chamber slice. Conversely, the three actual \(\delta_i\)-walls
impose the three raw halfspaces, so that slice cannot extend through any
facet of \(P\). Taking closures gives the exact equality claimed above.

## 6. Rational correspondences, dimension, and scope

Every raw reflection \(s_i\) is integral on \(T\). Acting as \(s_i\) on
\(T_{\mathbf Q}\) and as the identity on \(T_{\mathbf Q}^{\perp}\) gives

\[
\widetilde s_i=
\begin{pmatrix}s_i&(s_i-I_3)/15\\0&I_3\end{pmatrix}.
\]

The audit verifies that these matrices are involutive \(A\)-isometries,
self-adjoint, and commute with \(e\). Their denominators also verify the
negative statement: they do not preserve \(S\) integrally. Separately, the
ambient K3 reflection in \(\delta_i\) neither preserves \(T\) nor projects
to the raw reflection \(s_i\).

The divisor-product algebraic-correspondence claim is exactly rational and
requires no automorphism assertion. If \(D_i^\vee\) is the
intersection-dual divisor basis, then any rational endomorphism \(f\) of
\(\operatorname{NS}(X)_{\mathbf Q}\) is induced there by

\[
\Gamma_f=\sum_i D_i^\vee\times f(D_i),
\qquad
(\Gamma_f)_*(\alpha)
=\sum_i(\alpha.D_i^\vee)f(D_i)=f(\alpha).
\]

The independent inverse-Gram calculation verifies this formula for \(e\)
and all three \(\widetilde s_i\).

For this K3 category only,

\[
n+1=\dim E_1(e)\le\rho(X)\le20,
\qquad n\le19.
\]

The proved and unproved claims are therefore sharply separated:

- proved: one Picard-rank-six projected-root K3 carrier for
  \(\Delta(2,4,6)\), with exact eigen-nef return;
- not proved: literal full-nef Property 2 for \(k=\mathbf Q\);
- not proved: integral full-Picard extensions of the raw reflections;
- not proved: equality of ambient root reflections with raw reflections;
- not proved: universal nonemptiness or finite-host control for a fourth
  projected K3 branch; and
- not proved: any bound for an arbitrary non-K3 rational eigenslice category.

## 7. Replay

The independent artifacts are:

- [auditor](../scripts/audit_compact_246_projected_k3_eigenslice_carrier_independent_hostile.py);
- [archive](../data/enumerations/compact_246_projected_k3_eigenslice_carrier_independent_hostile_audit.json); and
- [tests](../tests/test_compact_246_projected_k3_eigenslice_carrier_independent_hostile_audit.py).

Run:

~~~text
python3 -m scripts.audit_compact_246_projected_k3_eigenslice_carrier_independent_hostile --check-output
pytest -q tests/test_compact_246_projected_k3_eigenslice_carrier_independent_hostile_audit.py
~~~
