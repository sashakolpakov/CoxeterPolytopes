# Multi-ideal triangle obstructions and a bound on simplex-face size

Two exact finite-volume Coxeter square pyramids show that triangle-port
transitivity can fail for triangles with exactly two ideal vertices, and
also for triangles with three ideal vertices. In both examples the chosen
cusp is simple and its whole link is a single simplex factor. Thus choosing
a simplex-factor face at that cusp does not repair the argument: the other
endpoint of the shared edge can still be a nonsimple ideal vertex.

A separate original-wall-star argument does give a justified bound on
aggregate size: **any actual simplex face of a finite-volume hyperbolic
Coxeter polytope has dimension at most9**. This does not assume that the
induced face angles are Coxeter angles. It sharpens the one-ideal cusp
block budget, without bounding the ambient dimension or covering the
remaining multi-ideal triangles.

The [standalone checker](../scripts/fresh_multi_ideal_triangle_scope_20261017.py)
and [certificate](fresh_multi_ideal_triangle_scope_20261017.json) use exact
Fraction pairs in Q(sqrt2). They reconstruct both pyramids, all vertices,
active sets, face incidences, Gram gates, and ideal-support semidefiniteness.
They also check the block maximum below by dynamic programming through
256 ports. Its proof for all port counts is given separately here.

## Two literal Coxeter counterexamples

Use the Lorentz form diag(-1,1,1,1). Put

\[
 b=(0,0,0,-1),\quad
 s_1=(a_1,1,0,a_1),\quad s_2=(a_2,0,1,a_2),\quad
 s_3=(a_3,-1,0,a_3),\quad s_4=(a_4,0,-1,a_4).
\]

All five normals have square1. Their products are

\[
 b\cdot s_i=-a_i,\qquad
 s_i\cdot s_{i+1}=0,\qquad s_1\cdot s_3=s_2\cdot s_4=-1.
\]

Indices on the successive sides are cyclic. Choose a_i=cos(pi/m_i) with
m_i in {3,4}. Then all native contacts are right, pi/3, pi/4 or parallel.
The common future center H=(2,0,0,1) has square-3 and pairs strictly
negatively with every normal.

In the chart t=1 these inequalities give the rectangular pyramid

\[
 0\le z\le1,\qquad
 -a_3(1-z)\le x\le a_1(1-z),\qquad
 -a_4(1-z)\le y\le a_2(1-z).
\]

Its apex is nu=(1,0,0,1). Its base vertices are the four combinations of
x in {a_1,-a_3} and y in {a_2,-a_4}, with z=0. Both choices below put all
vertices in the closed ball, with bounded Euclidean links at the ideal
vertices, hence give finite-volume Coxeter chambers.

**Two-ideal triangles.** Set (m_1,m_2,m_3,m_4)=(4,4,3,3). The base vertex

\[
 v=(1,1/\sqrt2,1/\sqrt2,0)
\]

is ideal; the other three base vertices are finite. The faces on s_1 and
s_2 are triangles, each with precisely the ideal vertices v and nu and
one finite vertex. They share edge v--nu. The third two-face through the
three ports at v is the quadrilateral base. Thus the two relations at v
do not complete to a third triangular relation.

**Three-ideal triangles.** Set every m_i=4. All five vertices are ideal.
The same two side triangles and base square give the same failure, now
with all three vertices of each triangle ideal.

In both cases v is incident with exactly three facets. Its link is the
affine (4,4,2) triangle and is combinatorially Delta2, so the selected
factor fiber is the entire link. The apex nu is incident with four facets;
its link is Delta1 x Delta1. Positive null relations include

\[
 \sqrt2\,b+s_1+s_2=\sqrt2\,v,\qquad
 s_1+s_3=(a_1+a_3)\nu,\qquad
 s_2+s_4=(a_2+a_4)\nu.
\]

The shared ideal edge therefore has a nonsimple endpoint nu. This is
exactly where the degree-three/Balinski argument fails. No inference from
a numerical near-ideal point or an unverified Gram signature is involved.

These examples rule out unrestricted transitivity and a disjoint-port
capacity1 for multi-ideal simplex blocks even within a simplex cusp fiber.
They do not rule out a different incidence inequality using additional
endpoint data.

## A simplex-face bound from the original Coxeter star

Let A be an actual k-dimensional simplex face of an n-dimensional
finite-volume Coxeter polytope P, with k >= 2. Every vertex of A is finite
or ideal. Let W be the n-k original ambient facet normals supporting A.
They are independent and span a positive-definite subspace.

Each facet of A is positive-dimensional, of dimension k-1. Simplicity at
edges gives exactly n-k+1 ambient facets containing it. Hence each facet
of A supplies **one** additional original normal b_i beyond W. These
k+1 black normals are distinct. This produces the original star

\[
 \Sigma=W\cup\{b_0,\ldots,b_k\},\qquad |\Sigma|=n+1.
\]

This does not say that an ideal vertex of P has only n incident facets.
The n normals obtained by deleting b_i are merely a selected subset of
the walls through the corresponding vertex of A.

Project the black normals orthogonally off W. The projected normals are
positive multiples of the geometric facet normals of A. Since A is a
nondegenerate projective k-simplex, its k+1 vertex rays are independent;
its facet normals form the corresponding dual basis. Their span is the
full Lorentzian (k+1)-space. Thus the Gram matrix of Sigma has signature
(n,1) and in particular is indefinite. No assertion about the projected
angles being Coxeter angles is used.

For every i, the star Sigma minus b_i consists of normals through an
actual vertex of A. It is positive definite at a finite vertex and
positive semidefinite at an ideal vertex, since the latter normals lie in
the orthogonal complement of a null vector. Therefore any indefinite
principal subdiagram of Sigma must contain **all k+1 black nodes**.

Choose an inclusion-minimal indefinite principal subdiagram L. It is
connected and every proper principal submatrix is positive semidefinite.
If every proper submatrix is positive definite, L is Lanner. Otherwise a
proper connected affine component C has a positive null vector. An
outside node adjacent to C makes C plus that node indefinite by its
nonzero coupling to that null vector. Minimality forces these to be all
of L. This also excludes a proper affine component accompanied by extra
orthogonal spherical or affine components. Thus L has precisely the
Lanner/quasi-Lanner scope of the standard classification.

The classification input is Prokhorov's
[§2, Proposition1, printed p.407](https://www.mathnet.ru/links/1e830e3298f5df6e148b8dc52fd64e1c/im1484_eng.pdf):
the largest quasi-Lanner diagrams have order10. The primary
[article record](https://www.mathnet.ru/eng/im1484) identifies the original
paper and translation. Its indexed primary PDF text was retrieved during
this check; a fresh full-PDF fetch timed out. The same exact input and the
minimal-indefinite reduction were already explicitly reviewed in the
[frozen compact-facet proof](fresh_nikulin_narrow_20261015_compact_facet.md).
This package does not independently reclassify those diagrams.

Since L contains all black nodes,

\[
 k+1\le |L|\le10,\qquad\boxed{k\le9.}
\tag{1}
\]

For k=1 the same conclusion is trivial. The local face-incidence facts
used to obtain W and the black normals are stated in
[Potyagailo–Vinberg, printed p.64](https://ems.press/content/serial-article-files/43032).

If the **ambient P is compact**, it has no affine subdiagram, so L is
Lanner and the familiar order-five bound instead gives k<=4. Compactness
of A alone does not justify that improvement: deleting a white node from
L can still leave an affine subdiagram containing all black nodes, even
when all black deletions are positive definite. The report makes no
order-five claim in that case.

## Sharper one-ideal cusp budget

The [one-ideal factor-fiber theorem](fresh_cusp_one_ideal_simplex_capacity_20261017.md)
gives disjoint blocks of k ports in a fiber with m=r_i+1 ports. By (1),
every nontrivial block has 2 <= k <= 9. Define

\[
 F_9(m)=36\lfloor m/9\rfloor+\binom{m\bmod9}{2}.
\tag{2}
\]

Then sum binom(k,2) <= F9(m). To prove this for every m, move a port from
a smaller nonempty block to a largest block of size below9. The objective
does not decrease, and increases when both blocks have size at least2.
A leftover block of size1 is an unused port and contributes zero. Repeat
until all blocks except the last have size9; the remainder has size m mod9.
This proves (2) without a rank cutoff. The exact dynamic-programming check
is an independent finite confirmation of this elementary optimization.

For a cusp with factor ranks r_i and V=product(r_i+1), every triangle with
exactly one ideal vertex belongs to one such block. Hence

\[
 \boxed{\quad
 2N_{3,1}(v)\le V\sum_i\frac{2F_9(r_i+1)}{r_i+1}.
 \quad}
\tag{3}
\]

For r_i<=8 the summand equals r_i. For larger factors it is at most8,
so in particular

\[
 2N_{3,1}(v)\le V\sum_i\min(r_i,8)
 \le\min((n-1)V,8s_vV).
\tag{4}
\]

This improves the previous worst-case bound 2E_v=(n-1)V when the cusp
has large factors. It does not pay for triangles with two or three ideal
vertices, or by itself provide the alpha_n E_v budget required by the
original-face inequality. The exact counterexamples above show why that
missing scope cannot be filled by simply repeating the transitivity proof.

The complementary
[finite-ended-port theorem](fresh_finite_ended_simplex_ports_20261017.md)
groups triangles having at least two finite vertices into maximal simplex
faces, which may have several ideal vertices. Those blocks are not the
same as the maximal one-ideal blocks in (3), so their incidence variables
must not be identified without an additional allocation argument.

**FINAL/FROZEN.** Generation and a fresh saved replay passed, as did source
hash, syntax and report-link checks. Exact_hub_counts independently read
the original-star/minimal-indefinite argument and the block-filling proof;
the primary classification was not independently re-enumerated. This is
an internally checked face-size lemma and two exact local counterexamples,
not a new ambient dimension ceiling, external review or novelty claim.
