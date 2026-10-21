# Hostile audit of the global minimality-14 integral-arrangement field funnel

**Date:** 2026-08-12  
**Verdict:** **PASS after scoped repairs.**

The audited theorem is a genuine unconditional **ground-field finiteness**
result in its stated scope:

\[
 P\subset\mathbb H^n,quad n\geq6,quad
 P\text{ finite-volume arithmetic Coxeter of simplest type}
 \quad\Longrightarrow\quad
 k(P)\text{ lies in a fixed finite set}.
\]

The proof does not use fixed-field/fixed-rank reflective-lattice finiteness,
the Vinberg dimension cutoff, or ABSW.  It also does not prove finiteness of
chambers, groups, facets, lattices, or global arrangement heights.  Those
nonclaims are essential: this theorem is real progress toward the requested
pushback, but it is not the terminal all-dimensional Coxeter-finiteness
theorem.

The initial draft had one false overstatement and several compressed gates.
The owner file now repairs them:

1. a dependent singular Gram packet is no longer claimed to inherit the
   entire field;
2. the tail packet explicitly proves \(k=F_d(\alpha)\);
3. the fork argument specifies its fixed base field and bounded exceptional
   embeddings;
4. the \(u<2\) branch uses a minimal positive-index packet rather than the
   false assertion that every minimal indefinite packet is Lanner;
5. the quadrilateral restriction uses \(k\)-rational root lines and a
   \(W\)-stable lattice; and
6. the global word “effective” and the suggestion of pure AG boundedness
   have been removed.

## 1. Full-field inheritance: exact surviving scope

Let \(K_I\) be the cycle field of a principal packet \(A_I\).  If
\(K_I\subsetneq k\), total reality gives a nonidentity embedding of \(k\)
which fixes \(K_I\), and it extends to the Gram coefficient field.  Fixed
two-cycles make every nonzero edge change only by sign; fixed longer cycles
make the edge-sign cocycle have trivial holonomy.  Hence

\[
                         A_I^\sigma=D A_I D
\]

for a diagonal sign matrix \(D\), so inertia is retained.

At a nonselected place the ambient form is definite.  Therefore:

- if \(A_I\) has a positive eigenvalue, its conjugate cannot be a Gram
  matrix in that definite space;
- if \(A_I\) is parabolic and its normals are linearly independent, their
  conjugates remain independent, so their conjugate Gram must be definite
  rather than singular.

This proves full-field inheritance for every packet actually used: dotted
hyperbolic pairs, tail/fork packets, four-node transition packets,
quasi-Lanner packets, and the independent two-node affine packet at \(u=2\).

The discarded stronger statement was false as a proof principle.  An
oversized dependent collection in a definite space can have a singular
negative-semidefinite Gram matrix, so singularity alone gives no
contradiction and need not force the full field.

## 2. Tail curve: norm and Fekete gates

For

\[
 A_d(u)=
 \begin{pmatrix}
 -2&u&0\\u&-2&2\cos(\pi/d)\\0&2\cos(\pi/d)&-2
 \end{pmatrix},
 \quad
 \alpha=u^2,quad \ell_d=4\sin^2(\pi/d),
\]

the determinant is \(2(\alpha-\ell_d)\).  The positive dotted pair gives

\[
 k=\mathbf Q(\alpha,\ell_d)=F_d(\alpha),
 \qquad F_d=\mathbf Q(\ell_d).
\]

Thus the relative degree and norm calculations really take place inside the
full ambient field, not a visible proper subfield.  At every nongeometric
embedding,

\[
                     0<\alpha^\tau<\ell_d^\tau,
\]

while \(4<\alpha_0<196\).

The cyclotomic norm is exact.  If \(\zeta_d\) is primitive, then

\[
 \ell_d=(1-\zeta_d)(1-\zeta_d^{-1}),\qquad
 N_{F_d/\mathbf Q}(\ell_d)=\Phi_d(1)=
 \begin{cases}
 p,&d=p^a,\\
 1,&d\text{ is not a prime power}.
 \end{cases}
\]

Consequently the product of interval capacities is

\[
 {\gamma(d)\over4^{\varphi(d)/2}}<1
\]

for every \(d\geq3\).  The standard Fekete theorem with one exceptional
geometric embedding bounds \([F_d(\alpha):F_d]\) for fixed \(d\).
Once the degree is bounded, all conjugates are bounded and the relative
minimal-polynomial coefficients occupy a finite Minkowski box.

For the uniform label bound, put
\(\beta=\alpha(\alpha-\ell_d)\).  It is nonzero and integral, and

\[
 1\leq |N_{k/\mathbf Q}(\beta)|
 <
 196^2{4\over\ell_{d,0}^2}
 \left({\gamma(d)^2\over4^{\varphi(d)/2}}\right)^{[k:F_d]}.
\]

Using

\[
 \ell_{d,0}\geq16/d^2,qquad
 \gamma(d)\leq d,qquad
 \varphi(d)^2\geq d/2
\]

reduces \(d\geq20000\) to the impossible inequality

\[
 1<601d^6/2^{\sqrt{d/2}}<1.
\]

The endpoint comparison
\(2^{100}>601\cdot20000^6\) and the claimed monotonicity are exact.
No omitted exceptional prime-power case was found.

## 3. Fork curve: the capacity product really contracts

Write

\[
 a=2\cos(\pi/s),\qquad b=2\cos(\pi/k).
\]

The determinant of the fork packet is

\[
 2\bigl(u^2+abu+a^2+b^2-4\bigr).
\]

At a nongeometric embedding, negative definiteness puts \(u\) strictly
between the two roots

\[
 { -ab\pm\sqrt{(4-a^2)(4-b^2)}\over2}.
\]

Their distance is

\[
 \sqrt{(4-a^2)(4-b^2)}
       =4|\sin(\pi/s)\sin(\pi/k)|<4.
\]

For \(s,k\in\{3,4,5\}\), let
\(E=\mathbf Q(a,b)\).  Under every embedding of this fixed totally real
field, \(|a|,|b|<2\), so every interval has capacity strictly below one and
their finite product is strictly below one.

Let \(K\) be the full cycle field and \(L=E(u)\).  The squares of the
three entries and their triangle product leave only vertex-sign ambiguity,
so the crude uniform bound \([L:K]\leq8\) is safe.  Hence at most eight
embeddings of \(L\) lie over the selected embedding of \(K\); on them
\(|u|<14\).  All other embeddings lie in the contracting determinant
intervals.  Fekete's theorem with a bounded number of exceptional intervals
therefore bounds \([L:E]\).  Since every conjugate of \(u\) is then
bounded, only finitely many \(u\), fields \(L\), and subfields \(K\)
occur.

This closes the most compressed arithmetic step in the draft.  The relevant
form of Fekete's theorem explicitly permits a fixed finite number of
exceptional embeddings; it is the version stated as Theorems 1.2.1--1.2.2
in Nikulin's primary account
[arXiv:0708.3991](https://arxiv.org/abs/0708.3991).

## 4. Five transition families

The frozen seven-support archive gives the following routes:

| Family | selected auxiliary | route |
|---|---:|---|
| \(\Gamma_1\) | node 3 | fork, both labels in \(\{3,4,5\}\) |
| \(\Gamma_2\) | node 3 | fork, both labels in \(\{3,4,5\}\) |
| \(\Gamma_3\) | node 3 | tail |
| \(\Gamma_4\) | node 4 | tail |
| \(\Gamma_5\) | node 3 | tail |

Every selected three-node packet contains the dotted pair \(u>2\), so it
has a positive eigenvalue and inherits the full field.  The tail label is
bounded by \(d<20000\); the two fork coefficient fields come from a finite
list.  Thus the fields of all five four-node transition families are finite.

The executable route check is not, by itself, a proof of the parameter
grammar: its declarations of which labels are bounded rely on the frozen
archive.  The independent archive reconstruction supplies that missing
mathematical evidence, so this is a test-coverage limitation rather than a
gap in the theorem.

## 5. Global edge and dimension-ten argument

For a noncompact simplest-type arithmetic chamber, nonuniformity makes its
defining \(k\)-form isotropic.  A nonzero isotropic \(k\)-vector remains
isotropic under every real embedding, which is incompatible with a
nonidentity definite place.  Hence \(k=\mathbf Q\).

In the compact case, vertex links are spherical Coxeter simplices and the
chamber is simple.  The minimal-visual-aperture cone calculation gives a
facet \(E\) whose closed star has normalized entries below \(14\).

For an edge with endpoint walls \(e_1,e_2\):

- \(e_1e_2=2\) is an independent parabolic pair and forces \(k=\mathbf Q\);
- if \(e_1e_2<2\), an inclusion-minimal positive-index packet inherits the
  full field; it is Lanner/triangle if all proper packets are elliptic, and
  quasi-Lanner otherwise.  The quasi-Lanner simplex is noncompact
  arithmetic and again has field \(\mathbf Q\).

Thus fields outside the rational, Lanner, and triangle branches have
\(e_1e_2>2\) on every edge.  If a finite vertex graph along \(E\) had a
component of size at least three, adjoining the opposite endpoint wall to a
connected three-node subgraph would produce one of the five transition
supports.  It cannot be \(\Gamma_5\), because deleting an endpoint in
\(\Gamma_5\) leaves only a two-node endpoint component.  Avoiding the
first four transition fields therefore makes every vertex graph a matching
plus isolated nodes.

The \(h\)-vector/Dehn--Sommerville bound for the average \(a\) of vertices
of a two-face of the \((n-1)\)-facet is

\[
 a<
 \begin{cases}
 4+4/(n-2),&n\text{ even},\\
 4+4/(n-3),&n\text{ odd}.
 \end{cases}
\]

Counting nonright corners gives

\[
 a\left(\binom{n-1}{2}+\lfloor(n-1)/2\rfloor\right)
 \geq5\binom{n-1}{2}.
\]

At \(n=10\) and \(11\), the lower bound equals \(4.5\) while the upper
bound is strict below \(4.5\); the contradiction strengthens thereafter.
This part is exact and agrees with the primary high-dimensional reduction.

## 6. Dimension six through the quadrilateral branch

The same strict average bound gives \(a<5\) for every \(n\geq6\), so the
narrow facet has a triangular or quadrilateral two-face.

A triangular face contradicts \(e_1e_2>2\) on one of its edges, because the
two endpoint-exclusive walls meet at the third vertex.  Hence its field is
already in the Lanner/triangle list.

For a quadrilateral \(M_4\), let \(Q\) be the \(n-2\) walls containing
its plane and cyclically label its four cutting walls
\(\delta_1,\ldots,\delta_4\).  A nonright corner makes
\(\{\delta_1,\delta_2\}\) a complete two-node component of the finite
vertex graph, so \(\delta_1,\delta_2\perp Q\).

If \(\delta_3\not\perp q\) for some \(q\in Q\), order the four nodes as

\[
                       (\delta_1,\delta_3;\delta_2,q).
\]

The first two are the dotted endpoint pair of the edge cut out by
\(Q\cup\{\delta_2\}\).  The remaining finite support is

\[
 \{13,23,24\}\sim\{13,14,23\}=\Gamma_4
 \quad\text{if }\delta_2\delta_3>0,
\]

or

\[
                         \{13,24\}=\Gamma_5
 \quad\text{if }\delta_2\delta_3=0.
\]

This is exactly the claimed exhaustive dichotomy.  Avoiding the five
transition fields forces \(\delta_3\perp Q\), and the cyclic argument
does the same for \(\delta_4\).

The orthogonal complement

\[
                         U=(\operatorname{span}_kQ)^\perp
\]

is a \(k\)-defined three-space of signature \((1,2)\) at the geometric
place and definite at every other place.  The \(k\)-definition is valid
because a \(k\)-rational reflection has a \(k\)-rational minus-one
eigenline.  Since every \(\delta_i\perp Q\), intrinsic projection does
not change any root, norm, angle, distance, or cyclic product.

Choose a full \(W\)-stable \(\mathcal O_k\)-lattice \(L\), obtained
from a reference lattice by taking a finite-index normal core and summing
finitely many translates.  Then

\[
                         L_U=L\cap U
\]

is a full rank-three torsion-free \(\mathcal O_k\)-module and hence
projective, and all four restricted reflections preserve it.  The
four-normal Gram spans \(U\), has rank three and a positive eigenvalue, so
its cycle field is the full ambient \(k\).  The compact quadrilateral gives
a finite-coarea reflection lattice in \(O(U)\); containment in
\(O(L_U)\) makes it arithmetic.  Thus this is an exact arithmetic Coxeter
quadrilateral, not merely a quasi-arithmetic projected face.

No rank or module obstruction remains in this step.

## 7. Polygon finiteness input

Every arithmetic Coxeter triangle or quadrilateral reflection lattice lies
in a maximal discrete reflection overgroup.  Indeed, a strict overgroup has
integer index at least two, so its chamber area drops by at least a factor of
two; the positive universal lower bound for finite-area hyperbolic
two-orbifolds makes every ascending chain terminate.  Finite extensions of
an arithmetic lattice remain arithmetic.

The independently audited genus-zero orbifold-curve theorem leaves only
finitely many maximal arithmetic reflection polygons.  For a fixed maximal
group, a triangle subgroup has index below

\[
                       \pi/\operatorname{Area}(P_{\max}),
\]

and a quadrilateral subgroup has index below

\[
                       2\pi/\operatorname{Area}(P_{\max}).
\]

A finitely generated group has finitely many subgroups of bounded index.
Therefore only finitely many triangle and quadrilateral groups, and hence
fields, occur.  The primary terminal theorem is Long--Maclachlan--Reid,
[Arithmetic Fuchsian Groups of Genus Zero](https://web.ma.utexas.edu/users/areid/genus0_final.pdf),
Theorem 1.1.

This argument proves finiteness.  The current repository does not enumerate
the resulting polygon list, so the owner theorem now correctly refrains from
calling the entire union effective.

## 8. AG attribution and remaining debt

The tail and fork inequalities are determinant equations in integral cycle
coordinates on exact quadric-arrangement charts.  In that precise sense the
result is arrangement-mediated and the bounded chart remembers the full
field.

The terminal tools are nevertheless arithmetic: Fekete for transition
integral points and arithmetic genus-zero orbifold-curve theory for polygons.
No BAB, KSBA, K-moduli, or general AG boundedness theorem proves the finite
field set.  More importantly,

\[
 \text{finite ground fields}
 \not\Longrightarrow
 \text{finite chambers or reflection groups}.
\]

The theorem leaves the dimension ceiling, remote facets and dotted heights,
native lattices, glue, and maximality-to-global-height implication
uncontrolled.  The owner report states this accurately.  It does not close
the project's remaining \(\texttt{MAX-ARR}\) gate.

## 9. Executable audit

The following completed successfully:

\[
\begin{array}{ll}
\texttt{PYTHONPATH=.:src python3 scripts/certify\_global\_minimality14\_field\_funnel.py}
 & \text{PASS},\\
\texttt{PYTHONPATH=.:src pytest -q tests/test\_global\_minimality14\_field\_funnel.py}\
\qquad\texttt{tests/test\_nikulin\_four\_node\_transition\_families.py}
 & 12\text{ passed},\\
\texttt{pandoc results/global\_minimality14\_integral\_arrangement\_field\_funnel.md}
 & \text{PASS}.
\end{array}
\]

The certificate checks the exact totient and cutoff arithmetic, the five
tail/fork routes, the dimension-ten matching contradiction, and the
dimension-six small-face inequality.  It does not prove Fekete, the
minimal-aperture geometry, or the rank-three lattice restriction; those were
audited above at the prose level.

Pinned SHA-256 values:

\[
\begin{array}{ll}
\texttt{05acee48bd7a6c03e92febc32434cc38bcd3581c0e9044dafec3f07f9dbab678}
 & \text{audited owner report},\\
\texttt{7296f74d01a7114de284b9a03c5141d37e769d4b6df97b9d4783e183adc93489}
 & \text{certificate script},\\
\texttt{2900482f5e8bb5bb0e1cfb35055545dd854a9e6276fc48ca6acbed931e01e194}
 & \text{certificate tests},\\
\texttt{f32af7c7b6ecfc00bcc6a2ac4fc0169468f46e715750a93f27370243cf161bb8}
 & \text{transition archive},\\
\texttt{68ab8d85002106fde4493540ee6c6409247211cdda5364a3181000ae0135b659}
 & \text{two-dimensional owner theorem}.
\end{array}
\]

These are the terminal audited pins for the theorem certificate.
