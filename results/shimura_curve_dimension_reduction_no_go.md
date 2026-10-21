# Genus-zero curve extraction from higher-dimensional arithmetic Coxeter chambers

## Conditional descent, an absolutely maximal counterexample, and the missing transverse bound

**Status (2026-08-12).** The successful real-DM orbifold-curve theorem in
dimension two does not extend to all dimensions through any of the proposed
automatic extractions: rank-three principal subdiagrams, minimal indefinite
subdiagrams, Coxeter faces, vertex links, or reflection-mirror centralizers.
There is a precise positive theorem when an arithmetic Coxeter two-face is
given as extra marked data. It inherits the ambient ground field and produces
the same exact genus-zero real-DM curve as in dimension two. Neither the
existence nor the arithmeticity of such a face is automatic.

The compact arithmetic tetrahedron \(H_6\) in Felikson's notation gives a
sharp control inside the target class. Its reflection group is absolutely
maximal, yet:

1. every rank-three principal Coxeter subdiagram is spherical;
2. the unique minimal indefinite principal packet has rank four;
3. none of its four intrinsic facets is a Coxeter polygon;
4. every vertex link is spherical; and
5. every reflection-mirror centralizer has compactified oriented coarse curve
   of genus one, not genus zero.

There is a second, independent obstruction to pushback. The absolutely
maximal arithmetic \(4\)-simplex \([5,3,3,3]\) and the absolutely maximal
arithmetic straight \(5\)-prism based on \([5,3,3,3,3]\) contain the same
arithmetic Coxeter triangle \((2,3,10)\), hence the same real-DM orbifold
curve \(\mathbf P^1_{(2,3,10)}\), over the same field. The curve therefore
does not bound the ambient dimension or recover the transverse form.

This is a no-go for **dimension reduction to genus-zero arithmetic curves
alone**. It is not a claim that no richer algebro-geometric carrier can ever
work.

## 1. The exact conditional face-curve theorem

Let \(P\subset\mathbb H^n\) be a finite-volume arithmetic Coxeter polytope
with ground field \(k\). Suppose a marked two-dimensional face
\(S\subset P\) satisfies both additional conditions:

- with its intrinsic metric and induced walls, \(S\) is a finite-area Coxeter
  polygon; and
- its reflection group \(W_S<\operatorname {Isom}(\mathbb H^2)\) is
  arithmetic, not merely quasi-arithmetic.

Write \(W_S^+=W_S\cap\operatorname {Isom}^+(\mathbb H^2)\). Doubling the
mirrored disk \(S\) gives the real-framed orbifold log curve

\[
 {\cal C}(S)=
 \left(
   \overline{[\mathbb H^2/W_S^+]}_{\mathbf R},
   D_{\rm cusp},c,\mathcal H
 \right).                                                   \tag{1.1}
\]

Here \(c\) is induced by the orientation-reversing coset and \(\mathcal H\)
is the chosen real half with its cyclic side marking. Its compactified coarse
curve is \(\mathbf P^1\). The real quotient, cone orders, punctures, and
cyclic marking recover the complete polygon \(S\), exactly as in the
two-dimensional theorem.

> **Theorem 1.1 (marked arithmetic Coxeter-face descent).** Under the two
> hypotheses above:
>
> 1. \({\cal C}(S)\) is an exact, functorial genus-zero arithmetic real-DM
>    curve for the marked pair \((P,S)\);
> 2. the ground field of \(S\) is exactly the ambient field \(k\);
> 3. the maximal arithmetic Fuchsian overgroups of \(W_S^+\) belong to the
>    finite Long--Maclachlan--Reid list; hence only finitely many fields and
>    Fuchsian commensurability classes can occur; and
> 4. if \(W_S\) is itself maximal among discrete reflection groups, then only
>    finitely many marked-isometry classes of \(S\) occur.

**Proof.** Exactness and genus zero follow by doubling the polygon, as in the
dimension-two construction. Bogachev--Kolpakov, Theorem 1.4, says that every
Coxeter face of a quasi-arithmetic Coxeter polytope is quasi-arithmetic with
the **same ground field**. Arithmetic implies quasi-arithmetic, so
\(k(S)=k(P)=k\).

Choose a maximal arithmetic Fuchsian lattice \(\Lambda^+\) containing
\(W_S^+\). The finite orbifold map

\[
       [\mathbb H^2/W_S^+]\longrightarrow
       [\mathbb H^2/\Lambda^+]                              \tag{1.2}
\]

induces a finite map from a coarse \(\mathbf P^1\); therefore the target
coarse curve also has genus zero. Maximal arithmetic Fuchsian groups are
congruence. Long--Maclachlan--Reid's genus-zero theorem leaves only finitely
many conjugacy classes of such \(\Lambda^+\), proving the field and
commensurability assertions. If \(W_S\) is reflection-maximal, the exact
two-dimensional theorem applies directly to \(S\). \(\square\)

The reflection-maximal hypothesis in the last clause matters. A fixed
arithmetic triangle group can have many finite-index reflection subgroups, so
the finite list of maximal Fuchsian overgroups does not by itself give a
finite list of nonmaximal face polygons.

### 1.1 Why the hypotheses are real restrictions

Bogachev--Kolpakov prove only quasi-arithmetic descent in general. Their
Theorem 1.6 gives a useful sufficient condition for an arithmetic facet: all
adjacent dihedral angles have the form \(\pi/(2m)\). They also compute a
Bugaenko example in which Coxeter faces of an arithmetic polytope are
properly quasi-arithmetic. Thus

\[
  P\text{ arithmetic},\quad S\text{ a Coxeter face}
  \centernot\Longrightarrow S\text{ arithmetic}.            \tag{1.3}
\]

The Long--Maclachlan--Reid theorem is a theorem about arithmetic Fuchsian
groups. It cannot be applied to a properly quasi-arithmetic face.

For an arbitrary \(k\)-rational totally geodesic \(2\)-subform the problem is
worse. It is not selected by the Coxeter chamber, its minimal trace field may
be a proper subfield of \(k\), its quotient genus is uncontrolled, and it
does not carry a return map to the wall system. The same-ground-field theorem
above is special to a Coxeter **face**, not to an arbitrary ternary subform.

## 2. An absolutely maximal arithmetic tetrahedron with no genus-zero Coxeter extraction

Put

\[
 k=\mathbf Q(\sqrt5),\qquad
 \varphi={1+\sqrt5\over2},
\]

and consider the norm-two Gram matrix

\[
G=
\begin{pmatrix}
2&-\varphi&0&-1\\
-\varphi&2&-1&0\\
0&-1&2&-1\\
-1&0&-1&2
\end{pmatrix}.                                               \tag{2.1}
\]

Its Coxeter diagram is a four-cycle with one edge labelled \(5\) and the
other three labelled \(3\). This is the compact tetrahedron \(H_6\) in
Felikson's Table 1.

### 2.1 Exact hyperbolicity, compactness, and arithmeticity

At the geometric embedding, the leading principal determinants are

\[
 2,\quad {5-\sqrt5\over2},\quad 3-\sqrt5,\quad
 -{1+5\sqrt5\over2}.                                        \tag{2.2}
\]

Hence \(G\) has signature \((3,1)\). Under the nonidentity embedding of
\(k\), they become

\[
 2,\quad {5+\sqrt5\over2},\quad 3+\sqrt5,\quad
 {-1+5\sqrt5\over2},                                        \tag{2.3}
\]

all positive. The four rank-three principal determinants are

\[
       3-\sqrt5,\quad 3-\sqrt5,\quad 4,\quad 4.              \tag{2.4}
\]

Their full Sylvester chains are positive. Thus every proper principal
subdiagram is elliptic and the simplex is compact. The entries are algebraic
integers, the \(5\)-edge cyclic product generates \(k\), and (2.3) is
positive definite. Vinberg's arithmeticity criterion therefore makes the
reflection group arithmetic over \(k\).

Felikson classified all Coxeter decompositions of bounded hyperbolic
tetrahedra. Her complete tables contain no nontrivial decomposition of
\(H_6\). A proper reflection overgroup would cut \(H_6\) into copies of its
smaller Coxeter chamber and hence give precisely such a decomposition.
Therefore the reflection group of \(H_6\) is maximal among all discrete
reflection groups. This classification input is used only to ensure that the
counterexample lies in the requested maximal class; it is not being used as
a finiteness theorem.

### 2.2 Rank-three and minimal-indefinite extraction fail

Equation (2.4), together with the positive lower principal minors, shows that
all four rank-three principal packets are spherical. There is no Lorentzian
rank-three Coxeter subdiagram from which to make a Fuchsian triangle group.
The only minimal indefinite principal subdiagram is the full rank-four
diagram (2.1), so the minimal-indefinite operation returns a hyperbolic
tetrahedron, not a curve.

Every vertex link is also spherical: its Coxeter diagram is one of those
rank-three elliptic subdiagrams. Link recursion therefore never produces a
hyperbolic curve on this example.

### 2.3 None of the intrinsic facets is a Coxeter polygon

Deleting a diagram node is not the correct intrinsic facet Gram matrix when
the deleted wall is not orthogonal to all neighboring walls. If node \(i\)
is fixed and \(J\) is the set of remaining nodes, the projected facet-normal
Gram matrix is the Schur complement

\[
                 G^{(i)}=G_{JJ}-G_{Ji}G_{ii}^{-1}G_{iJ}.     \tag{2.5}
\]

For a pair of projected normals put

\[
 c_{ab}^{(i)}=-{G^{(i)}_{ab}\over
        \sqrt{G^{(i)}_{aa}G^{(i)}_{bb}}},\qquad
 \mu_{ab}^{(i)}={\pi\over\arccos c_{ab}^{(i)}}.              \tag{2.6}
\]

Exact cosine comparisons give the following open intervals for the three
orders of the four triangular facets:

| fixed node | exact order intervals |
|---|---|
| \(0\) | \((5,6),(4,5),(3,4)\) |
| \(1\) | \((4,5),(5,6),(3,4)\) |
| \(2\) | \((8,9),(3,4),(2,3)\) |
| \(3\) | \((8,9),(2,3),(3,4)\) |

Every value lies strictly between consecutive integers. In particular no
intrinsic facet angle is \(\pi/m\) with integral \(m\), so none of the four
facets is a Coxeter polygon. This also illustrates why a rank-three
principal subdiagram cannot be silently identified with a geometric facet.

### 2.4 The mirror-centralizer curve has genus one

All four finite Coxeter labels in the diagram are odd, while the two missing
edges have label \(2\). Therefore the odd Coxeter diagram is \(C_4\). It is
connected, so all simple reflections, and hence all reflections, are
conjugate. Its first Betti number is

\[
                  b_1(C_4)=4-4+1=1.                         \tag{2.7}
\]

Brink's theorem in Allcock's explicit topological form gives, for any
reflection \(s\),

\[
 C_W(s)/\langle s\rangle
       =W_\Omega\rtimes\Gamma_\Omega,\qquad
 \Gamma_\Omega\simeq\pi_1(C_4)\simeq\mathbf Z.              \tag{2.8}
\]

In Allcock's chamber model, the effective quotient of the mirror by the full
centralizer is obtained from the tetrahedral boundary by deleting the two
opposite even ridges and the incident zero-strata. Its natural bordered
compactification is an annulus. Doubling the annulus across its two mirrored
boundary components gives

\[
           g_{\rm double}=2g+b-1=2\cdot0+2-1=1.             \tag{2.9}
\]

The deck group \(\Gamma_\Omega\) is orientation-preserving on the mirror:
each odd-edge transport is a power of a product of two ambient reflections,
and it fixes the selected normal half-line. Hence the index-two
orientation-preserving quotient is exactly this double, rather than a
further unoriented quotient.

Because the odd diagram is connected, (2.9) applies to every reflection
mirror. The effective orientation-preserving stabilizer is an arithmetic
Fuchsian lattice: it is finite index in the arithmetic stabilizer of the
\(k\)-rational mirror. It therefore defines an arithmetic algebraic orbifold
curve, but its compactified coarse genus is one. The automatic genus-zero
premise needed by the Long--Maclachlan--Reid reduction is absent. A larger
Fuchsian overgroup may happen to have genus zero in a particular case, but
genus one does not force that; no genus-zero quotient follows from reflection
maximality in \(\operatorname {Isom}(\mathbb H^3)\).

Combining Sections 2.2--2.4 proves the promised scoped no-go:

> **Theorem 2.1 (failure of the natural universal curve extractions).**
> There exists an absolutely reflection-maximal compact arithmetic Coxeter
> polytope for which rank-three principal extraction, minimal-indefinite
> extraction, Coxeter-facet recursion, vertex-link recursion, and
> reflection-mirror centralizers produce no forced genus-zero arithmetic
> Coxeter/Fuchsian curve. The polytope is \(H_6\) with Gram matrix (2.1).

This theorem does not formalize the word "canonical" so broadly that a rule
may choose an unrelated arithmetic curve by fiat. It disproves the intrinsic
constructions that retain a geometric or group-theoretic return map to the
chamber.

## 3. Even a successful genus-zero face does not bound the ambient form

Consider the compact arithmetic \(4\)-simplex \(P_4=[5,3,3,3]\). With the
same norm-two convention its path Gram matrix has off-diagonal entries
\(-\varphi,-1,-1,-1\). Intersect the first two marked facets. The Schur
complement on the remaining three walls is

\[
C=
\begin{pmatrix}
1-\sqrt5/5&-1&0\\
-1&2&-1\\
0&-1&2
\end{pmatrix}.                                               \tag{3.1}
\]

After diagonal normalization to norm two, the three off-diagonal values are

\[
             -2\cos(\pi/10),\qquad -1,\qquad 0.              \tag{3.2}
\]

Thus the intrinsic two-face is the Coxeter triangle
\(\Delta(2,3,10)\), and its exact real-DM curve is

\[
             \mathbf P^1_{(2,3,10)}\quad\text{over}
             \quad\mathbf Q(\sqrt5),                         \tag{3.3}
\]

with its chosen real half.

The compact straight Coxeter \(5\)-prism \(P_5\) based on
\([5,3,3,3,3]\) has an isometric \([5,3,3,3]\) base facet. It therefore
contains exactly the same marked triangle and exactly the same curve (3.3).

Both ambient reflection groups are absolutely reflection-maximal. The
orientation subgroup of \(P_4\) realizes the unique minimum-covolume compact
orientable arithmetic \(4\)-orbifold; the corresponding statement for
\(P_5\) is the Emery--Kellerhals minimum theorem in dimension five. If either
reflection group had a proper discrete reflection overgroup, its orientation
subgroup would be a proper arithmetic overgroup of smaller covolume, a
contradiction.

Consequently the genus-zero curve, its cone orders, real half, and ground
field do not determine even the ambient dimension:

\[
 (P_4,\Delta(2,3,10))\longmapsto\mathbf P^1_{(2,3,10)}
 \longmapsfrom(P_5,\Delta(2,3,10)).                          \tag{3.4}
\]

Adding the normal representation would distinguish these two objects, but
then the missing quantity is precisely its rank. The genus-zero curve theorem
supplies no bound on that rank. Likewise, a packet of all extracted curves
plus incidence data might eventually recover the diagram, but its
cardinality, multiplicities, normal ranks, and transverse gluings must first
be bounded. Assuming such a bound is assuming the higher-dimensional part
of the desired theorem.

## 4. Audit of the proposed reductions

| extraction | what is inherited | exact obstruction |
|---|---|---|
| Lorentzian rank-three principal subdiagram | an arithmetic triangle when it exists | \(H_6\) has none |
| minimal indefinite principal subdiagram | a hyperbolic Coxeter packet | the minimal packet of \(H_6\) has rank four |
| intrinsic two-dimensional Coxeter face | same ground field by Bogachev--Kolpakov | need not exist; if it exists it may be properly quasi-arithmetic |
| facet recursion | quasi-arithmeticity for those faces that remain Coxeter | facets need not be Coxeter; all four \(H_6\) facets fail |
| vertex/link recursion | local spherical or Euclidean Coxeter data | compact \(H_6\) links are spherical, not Fuchsian |
| reflection-mirror centralizer | a natural arithmetic Fuchsian stabilizer | its curve need not have genus zero; \(H_6\) gives genus one |
| arbitrary ternary \(k\)-subform | an arithmetic Fuchsian group after integral choices | no functorial choice, no Coxeter return, genus and minimal trace field uncontrolled |
| one successful genus-zero face | ambient field and one Fuchsian commensurability class | \(P_4\) and \(P_5\) share the same \((2,3,10)\) curve |

The Long--Maclachlan--Reid theorem can bound the maximal Fuchsian group of a
successful arithmetic genus-zero face. It does not bound:

1. the ambient Lorentzian rank \(n+1\);
2. the number of ambient walls;
3. the normal representation along the curve;
4. the remaining Gram entries; or
5. the incidence/gluing of several curve packets.

To turn curve finiteness into finiteness of the original Coxeter polytopes one
would need a genuinely new theorem giving all three of the following:

\[
\begin{array}{ll}
\text{(E)}&\text{existence of a functorial arithmetic genus-zero two-section;}\\
\text{(N)}&\text{a uniform bound on its normal rank and on the number of curve orbits;}\\
\text{(G)}&\text{finiteness of transverse Gram extensions/gluings.}
\end{array}                                                  \tag{4.1}
\]

The \(H_6\) control disproves the natural versions of (E), and (3.4) shows
that the curve theorem alone gives neither (N) nor (G). Hence the known
two-dimensional AG theorem does not yield all-dimensional reflection-group
finiteness without an independent higher-dimensional input.

## 5. Reproducibility and primary sources

The exact determinant, Schur-complement, cosine-interval, odd-graph, and
\((2,3,10)\)-face checks are implemented in
scripts/certify_shimura_curve_dimension_reduction_no_go.py. The regression
tests are in tests/test_shimura_curve_dimension_reduction_no_go.py.
Absolute maximality is explicitly labelled as a literature input rather than
being smuggled into the symbolic certificate.

Primary sources:

- N. Bogachev and A. Kolpakov,
  [*On faces of quasi-arithmetic Coxeter polytopes*](https://arxiv.org/abs/2002.11445),
  especially Theorems 1.4 and 1.6 and the Bugaenko example.
- D. Long, C. Maclachlan, and A. Reid,
  [*Arithmetic Fuchsian Groups of Genus Zero*](https://web.ma.utexas.edu/users/areid/genus0_final.pdf),
  Theorem 1.1.
- A. Felikson,
  [*Coxeter Decompositions of Hyperbolic Tetrahedra*](https://arxiv.org/abs/math/0212010),
  Table 1 and the complete bounded-decomposition tables.
- D. Allcock,
  [*Reflection Centralizers in Coxeter Groups*](https://arxiv.org/abs/1106.3822),
  refining Brink's centralizer theorem and identifying
  \(\Gamma_\Omega\) with the fundamental group of the relevant odd-diagram
  component.
- M. Belolipetsky,
  [*On Volumes of Arithmetic Quotients of \(SO(1,n)\)*](https://arxiv.org/abs/math/0306423)
  and its [addendum](https://arxiv.org/abs/math/0610177), for the unique
  compact arithmetic minimum in even dimension; N. Bredon and R. Kellerhals,
  [*Hyperbolic Coxeter Groups and Minimal Growth Rates in Dimensions Four and Five*](https://arxiv.org/abs/2008.10961),
  for the identification with \([5,3,3,3]\) and the straight \(5\)-prism.
- V. Emery and R. Kellerhals,
  [*The Three Smallest Compact Arithmetic Hyperbolic 5-Orbifolds*](https://arxiv.org/abs/1205.2984),
  for the dimension-five minimum and its Coxeter realization.
