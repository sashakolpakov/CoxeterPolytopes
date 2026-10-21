# The anti-symplectic fixed chamber and the Bugaenko facet ladder

> **Subsequent decision.**  The nontrivial diagram reversal of the complete
> Bugaenko facet extends to a symmetry `A` of this complete H9 chamber.  With
> the intrinsic Pell-fold symmetry `g`, it satisfies
> `tr(A*g)=8+8*phi>10`; hence the stabilizer is infinite and this H9 chamber,
> as well as its parent H10 chamber, has infinite volume.  See
> [`og10_golden_h9_infinite_stabilizer.md`](og10_golden_h9_infinite_stabilizer.md).

## Result

Let

\[
H=\begin{pmatrix}-2\phi&1\\1&2\end{pmatrix}
  \oplus E_8\oplus\langle2\rangle
\]

over \(\mathcal O=\mathbb Z[\phi]\), and let \(C\subset\mathbb H^{10}\)
be the complete norm-two chamber constructed in the golden OG10 carrier.
The complete-chamber involution is

\[
h=\begin{pmatrix}1&1-\phi\\0&-1\end{pmatrix}\oplus I_9.
\]

It is the reflection in

\[
\nu=(\phi-1,2,0^9),\qquad \nu^2=6+2\phi,
\]

and \(hC=C\).  Its fixed field lattice is exactly

\[
L_9=H^h=\langle-2\phi\rangle\oplus E_8\oplus\langle2\rangle,
\tag{1}
\]

of signatures \(((9,1),(10,0))\).  The fixed slice

\[
P_9=C\cap\mathbb H(L_9)
\tag{2}
\]

is a chamber of the **complete** norm-two reflection arrangement of \(L_9\).
No projected or folded walls have to be added: every nonfixed simple-wall
orbit is globally ultraparallel to the fixed hyperplane and restricts to a
redundant inequality.

There is a second exact reduction.  The \(u=0\) facet of \(P_9\) is precisely
Bugaenko's compact eleven-wall polytope \(P_8\subset\mathbb H^8\).  Thus the
OG10 construction contains a successive right-angled facet ladder

\[
                 P_8\ \subset\ P_9\ \subset\ \widehat P_{10},
\tag{3}
\]

where

\[
\widehat P_{10}=C\cap\{(x,\nu)\leq0\}
\]

is a fundamental chamber for \(\langle W_2,h\rangle\).  Here \(P_9\) is the
\(\nu\)-facet of \(\widehat P_{10}\), and \(P_8\) is the final
\(\langle2\rangle\)-coordinate facet of \(P_9\).  The walls that meet either
facet meet it orthogonally; all remaining simple walls are parallel or
ultraparallel to it.

This is a structural theorem, not a new finite-volume polytope.  Subsequent
complete-facet symmetry analysis proves that \(P_9\) and \(C\) have infinite
volume and are noncompact; finite-sidedness remains open.

## 1. Why nonfixed wall orbits disappear

Write an ambient norm-two root as

\[
r=(x,y,z,u),\qquad z\in E_8\otimes\mathcal O.
\]

Direct calculation gives

\[
h(r)-r=-y\nu,
\qquad
(r,h(r))=2-(3+\phi)y^2,
\tag{4}
\]

and therefore

\[
(r+h(r))^2=8-2(3+\phi)y^2.
\tag{5}
\]

Suppose that \(r\) is a nonfixed simple root of the complete chamber.  Since
\(hC=C\), the distinct root \(h(r)\) is another oriented simple root.  The
simple system is acute.  If the two walls were not ultraparallel, then

\[
-2\leq(r,h(r))\leq0,
\]

so (4) would imply

\[
\frac2{3+\phi}\leq y^2\leq\frac4{3+\phi}.
\tag{6}
\]

At the conjugate embedding \(\phi'=1-\phi\), the form is positive definite.
Minimizing its binary block over \(x'\), with \(y'\) fixed, gives

\[
2=(r',r')\geq\frac{4-\phi}{2}(y')^2,
\qquad
(y')^2\leq\frac4{4-\phi}.
\tag{7}
\]

Write \(y=a+b\phi\), with \(a,b\in\mathbb Z\).  The exact inequalities

\[
\frac4{3+\phi}<1,
\qquad
\sqrt{\frac4{3+\phi}}+
\sqrt{\frac4{4-\phi}}<\sqrt5
\tag{8}
\]

force

\[
|b|\sqrt5=|y-y'|<\sqrt5.
\]

Hence \(b=0\).  But then (6) requires a nonzero integer \(a\) with
\(|a|<1\), a contradiction.  The square-root comparison in (8) is exact:
after moving the positive terms and squaring, its margin is \(25/121>0\).

It follows that

\[
                        (r,h(r))<-2                 \tag{9}
\]

for every nonfixed simple root, including every wall beyond the computed
prefixes.  For a fixed vector \(v\),

\[
(v,r)=(v,h(r))=\tfrac12(v,r+h(r)).
\]

By (5) and (9), \(r+h(r)\) is timelike.  Its orthogonal hyperplane contains
no point of \(\mathbb H(L_9)\).  Since the fixed interior point lies on the
chosen side, both inequalities in this orbit are redundant on the entire
fixed hyperbolic component.

The only walls remaining in (2) are therefore the \(h\)-fixed ambient simple
roots.  They are exactly the simple walls of the norm-two chamber of (1):
every norm-two root of \(L_9\) is an ambient fixed root, while the component
of the ambient root complement containing the fixed interior point cannot be
cut internally by any such mirror.  This proves the complete, rather than
prefix-level, fixed-slice assertion.

## 2. Bugaenko's compact polytope is a complete facet

Use coordinates \(r=(x,z,u)\) on \(L_9\).  The height-zero root

\[
s_u=(0,0,-1)
\]

is the wall normal for the final \(\langle2\rangle\) coordinate.  For every
other oriented simple root, acuteness gives

\[
(s_u,r)=-2u\leq0,
\qquad u\geq0.                                      \tag{10}
\]

Conjugate positive definiteness gives \(|u'|\leq1\).  If \(0<u<1\) and
\(u=a+b\phi\), then

\[
|b|\sqrt5=|u-u'|<2,
\]

so \(b=0\); no integer \(a\) lies strictly between zero and one.  Therefore

\[
u=0\quad\hbox{or}\quad u\geq1.                    \tag{11}
\]

In the second case, (10) is at most \(-2\): the wall is parallel or
ultraparallel to the \(s_u\)-facet and cannot cut its hyperbolic interior.
The cutting walls are exactly those with \(u=0\), which are the roots of

\[
L_8=\langle-2\phi\rangle\oplus E_8.
\]

The independently completed Vinberg transcript for this lattice gives
Bugaenko's compact \(\mathbb H^8\) chamber: 11 facets, 41 vertices, and

\[
(f_0,\ldots,f_8)=(41,164,316,374,294,156,54,11,1).
\]

The eleven \(u=0\) roots found inside the present fixed slice have, after an
explicit permutation, exactly the archived Bugaenko wall Gram matrix.  The
first nine roots in this order have coordinate determinant
\(2\phi-3\), whose field norm is \(-1\).  They therefore form an
\(\mathcal O\)-basis of \(L_8\), so the comparison is an integral lattice
identification, not merely a real diagram isometry.  Thus this is not merely
a facet isometric to an unknown subcone: it is the complete Bugaenko chamber.

## 3. The right-angled ladder and the dimension-reduction gate

Because \(h\) is reflection in \(\nu\) and normalizes \(W_2\), cutting the
\(W_2\)-chamber \(C\) by one side of \(\nu^\perp\) gives a fundamental
chamber \(\widehat P_{10}\) for \(\langle W_2,h\rangle\).  If a simple root
is fixed by \(h\), then it is orthogonal to \(\nu\).  If it is not fixed,
(9) is equivalent to

\[
\frac{(r,\nu)^2}{r^2\nu^2}>1,
\]

so its wall is ultraparallel to the \(\nu\)-facet.  This proves the first
right-angled step in (3).  Equations (10)--(11) prove the second.

The ladder gives a useful one-way decision theorem.  The ambient arithmetic
group is uniform, so any finite-volume fundamental polytope for a subgroup of
finite covolume is compact.  Consequently

\[
\operatorname{vol}(C)<\infty
\Longrightarrow \widehat P_{10}\text{ compact}
\Longrightarrow P_9\text{ compact}.                \tag{12}
\]

Therefore a proof that the complete \(L_9\) chamber is nonreflective,
infinite-volume, or noncompact immediately proves that the OG10
\(\mathbb H^{10}\) chamber has infinite volume.  This makes the smaller
fixed-lattice Vinberg problem the correct next primal computation.  The
converse is false: compactness of \(P_9\) is necessary, not sufficient, for
compactness of \(\widehat P_{10}\) or \(C\).

## 4. The 61-wall fixed prefix still escapes

Filtering the exact 1,144-wall ambient enumeration through control coordinate
\(4+5\phi\) gives 61 fixed roots.  The earlier involution-closed 640-root set
had 44 fixed roots and 298 nontrivial pairs.  Its complete exact cone census
was

| finite rays | cusps | hyperideal rays | opposite-sheet rays |
|---:|---:|---:|---:|
| 784 | 0 | 3113 | 0 |

with the concise hyperideal witness

\[
(1,2\phi,4\phi,6\phi,8,6,4,2,3\phi,5-3\phi),
\qquad v^2=178-110\phi>0.
\]

One of the seventeen new fixed walls cuts that ray.  The exact replacement
for the 61-wall cone has active facets

\[
(0,1,2,4,5,6,12,48,49)
\]

and ray

\[
v_{61}=\left(
1,2\phi,4\phi,6\phi,
\frac{64}{13}\phi,\frac{48}{13}\phi,
\frac{32}{13}\phi,\frac{16}{13}\phi,
\frac{40}{13}\phi,\frac1{13}\phi
\right),
\]

with

\[
v_{61}^2=\frac{210-128\phi}{169}>0.       \tag{13}
\]

All 61 wall pairings are nonpositive and exactly the nine listed facets are
active.  Thus the present 61-wall outer cone is not finite-volume.  Equation
(13) does **not** prove that the complete chamber escapes: later fixed simple
roots may cut this ray.

## 5. The algebro-geometric involution

The fixed slice has an honest OG10 realization, not just a lattice-theoretic
interpretation.  Extend \(h\) by \(-I\) on the period complement
\(\langle2\rangle\oplus\langle66\rangle\), as in the complete-chamber
involution theorem.  Replaying the order-44 discriminant gluing and taking
the integral kernels of \(h\mp I\) in the full Beauville lattice gives

\[
H^2(X,\mathbb Z)^h
 \cong E_8(-1)^2\oplus\langle2\rangle\oplus\langle-2\rangle^3,
\tag{14}
\]

and

\[
H^2(X,\mathbb Z)_h
 \cong\langle2\rangle^2\oplus\langle-2\rangle\oplus\langle-6\rangle.
\tag{15}
\]

The exact invariant-lattice audit gives rank 20, signature \((1,19)\),
determinant \(-16\), and Smith diagonal
\(1^{16},2^4\).  The coinvariant lattice has rank 4, signature \((2,2)\),
determinant 48, and Smith diagonal \((2,2,2,6)\).  Their orthogonal direct sum
has index 16 in the OG10 Beauville lattice.  These data, including the
nontrivial discriminant action, match row 3 of Billi--Grossi's Table A.7.

It remains to ensure that the birational involution can be made regular
without abandoning the fixed slice.  The visible OG10 Kähler-wall norm
orbits are represented by

\[
2,\qquad4,\qquad6-2\phi.
\]

If such a wall hyperplane contained \(\operatorname{Fix}(h)\), its normal
would be proportional over \(\mathbb Q(\sqrt5)\) to \(\nu\).  The ratio of
its norm to \(\nu^2=6+2\phi\) would therefore be a square in the field.  The
rational norms of the three ratios are

\[
\frac1{11},\qquad\frac4{11},\qquad\frac5{11},
\tag{16}
\]

none of which is a square in \(\mathbb Q\).  Since the norm of a square in a
quadratic field is a rational square, (16) excludes every containing Kähler
wall.

Choose a generic point of \(P_9\), avoiding the remaining locally finite
wall restrictions.  It lies in a unique Kähler-type chamber.  Because \(h\)
fixes the point and permutes the wall arrangement, that chamber is
\(h\)-stable.  Markman's Kähler-type chamber theorem realizes it as the
Kähler cone of a same-period birational model.  The positive Hodge-monodromy
extension of \(h\) preserves that Kähler cone, so global Torelli realizes it
as a biregular automorphism.  The representation
\(\operatorname{Aut}(X)\to O(H^2(X,\mathbb Z))\) is faithful for OG10-type
manifolds by Mongardi--Wandel; since the cohomological action has order two,
the automorphism is an involution.  It acts by \(-1\) on \(H^{2,0}\), hence
is anti-symplectic.

Thus \(P_9\) is the golden identity-embedding section of an invariant
eigen-movable cone for an actual anti-symplectic OG10 involution.  It is not
the full positive cone of (14), which has hyperbolic dimension 19, and no
claim about a smooth quotient or its nef cone is needed.

## 6. Scope and next step

The proved statements are:

1. the complete fixed slice is exactly the norm-two chamber of
   \(\langle-2\phi\rangle\oplus E_8\oplus\langle2\rangle\);
2. every nonfixed ambient wall orbit is globally redundant on that slice;
3. Bugaenko's compact \(\mathbb H^8\) polytope is a complete right-angled
   facet of it;
4. the current 61-wall outer cone still has a certified hyperideal escape;
5. the involution regularizes on a same-period OG10 model and has the exact
   invariant/coinvariant lattices (14)--(15); and
6. noncompactness of the complete fixed chamber would settle the original
   \(\mathbb H^{10}\) chamber negatively; the subsequent infinite-stabilizer
   theorem supplies exactly this conclusion.

Finite-sidedness remains unproved, and no new compact Coxeter polytope is
claimed.  Subsequent symmetry analysis proves infinite volume and
noncompactness of both \(P_9\) and \(C\).  A direct rank-10 Vinberg recursion
has also replayed the 61 roots and exhausted the next control coordinate
\(3+6\phi\), adding ten roots; every one is strict on the unchanged
hyperideal ray.  See
[`og10_golden_h9_vinberg_through_3_plus_6phi.md`](og10_golden_h9_vinberg_through_3_plus_6phi.md)
for that bounded primal corroboration.  Further shells can study
finite-sidedness but are no longer a finite-volume decision gate for this
carrier.

## Sources and reproducibility

The complete ambient symmetry and its Hodge-monodromy extension are proved in
the [complete-chamber involution report](og10_golden_complete_chamber_involution.md).
The OG10 wall types and Kähler-chamber realization are in the
[complete wall/movable-cone report](og10_golden_complete_wall_chambers.md).
Bugaenko completion is independently archived in
[`data/known_examples/bugaenko_vinberg_transcript.json`](../data/known_examples/bugaenko_vinberg_transcript.json)
and [`data/known_examples/bugaenko_compact_h8.json`](../data/known_examples/bugaenko_compact_h8.json).

The external geometric inputs are Markman's
[Kähler-type chamber and Torelli results](https://arxiv.org/abs/1101.4606),
Mongardi--Onorati's
[OG10 wall classification](https://arxiv.org/abs/2010.12511), and
Billi--Grossi's
[classification of non-symplectic prime-order OG10 automorphisms, Table A.7](https://doi.org/10.1093/imrn/rnaf159).
Faithfulness of the action of OG10 automorphisms on \(H^2\) is due to
Mongardi--Wandel,
[Theorem 3.1](https://arxiv.org/abs/1411.0759).

The exact generator, frozen certificate, and regression test are:

- [`scripts/certify_og10_golden_h9_fixed_slice.py`](../scripts/certify_og10_golden_h9_fixed_slice.py),
- [`data/candidates/og10_golden_h9_fixed_slice.json`](../data/candidates/og10_golden_h9_fixed_slice.json), and
- [`tests/test_og10_golden_h9_fixed_slice.py`](../tests/test_og10_golden_h9_fixed_slice.py).

Run

```sh
PYTHONPATH=.:src python3 scripts/certify_og10_golden_h9_fixed_slice.py
PYTHONPATH=.:src pytest -q tests/test_og10_golden_h9_fixed_slice.py
```

The optional historical full-cone replay takes longer:

```sh
PYTHONPATH=.:src python3 scripts/certify_og10_golden_h9_fixed_slice.py --full-census
```
