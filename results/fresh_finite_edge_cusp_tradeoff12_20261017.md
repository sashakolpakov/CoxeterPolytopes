# Eight-simplex faces and finite-edge triangle coefficient12

Every actual eight-dimensional simplex face of a finite-sided,
finite-volume acute Coxeter polytope has at most three ideal vertices.
This statement concerns the face as a geometric simplex; its induced
angles are not assumed to be Coxeter angles. The ambient chamber may
have nonsimple ideal vertices and arbitrarily many other supporting
facets.

Together with the frozen nine-simplex result, it proves

\[
 2(T_0+T_1)\le12E_{ff}
 =6n v_f-6E_{fi}\le6n v_f.                       \tag{1}
\]

Here T_s counts triangular two-faces with exactly s ideal vertices;
E_ff counts finite-to-finite edges, E_fi finite-to-ideal edges, and v_f
finite vertices. This is a complete internally reviewed lemma, with
exact arithmetic checks described below. It is not an ambient
dimension bound or a sharpness claim for coefficient12.

## Geometric inputs and the34 cases

The [original-star theorem](fresh_multi_ideal_triangle_scope_20261017.md)
gives dimension at most9 for every actual simplex face. An eight-simplex
S has nine original black facet normals B and n−8 positive-definite
supporting white normals W. The star W∪B is indefinite. Every black
deletion is a selected Gram matrix at an actual vertex and is positive
semidefinite; it is positive definite exactly when that vertex is finite.
This uses simplicity at positive-dimensional faces, not simplicity of
every ambient ideal vertex.

A minimal indefinite subdiagram L must contain every black root. Its
order is9 or10. The [frozen exact census](fresh_high_rank_quasilanner_census_20261017.json)
contains four order-nine and three order-ten types, in agreement with
[Chein's primary tables, printed pp.15–16](https://www.numdam.org/item/M2AN_1969__3_3_3_0.pdf).
The full census generated146,262 affine attachments, retained84 labeled
attachments, and found26 types of orders6–10; its source and saved data
are bound here. The present checker does not replay that entire census.
It independently rebuilds the seven required high-rank matrices and all
their deletion tests from literal edge lists using Fraction arithmetic
over Q(sqrt2).

For order9, L=B. For order10, there is one internal white root w∈L\B.
There are therefore four plus3·10=34 cases. The result of every case is
saved individually, including every diagonal threshold and its minimizing
ties. Census indices used below are literal data identifiers:

| L | Internal white choices | Number of cases | Upper bound on ideal vertices |
|---|---|---:|---:|
| 19 | none | 1 | 2 |
| 20 | none | 1 | 2 |
| 21 | none | 1 | 3 |
| 22 | none | 1 | 2 |
| 23 | 0,…,9 | 10 | 2 |
| 24 | 0,1,9 | 3 | 2 |
| 24 | 2,…,8 | 7 | 3 |
| 25 | 0,…,8 | 9 | 2 |
| 25 | 9 | 1 | 2 |

These are upper bounds, not assertions that every boundary value is
realized in every ambient completion.

## External whites when a black affine deletion exists

Let X=W\L be the external whites. If L\{b} is affine for a black b,
its kernel has strictly positive coordinates on all its nodes. The
actual vertex Gram on W∪(B\{b}) is positive semidefinite. A zero
quadratic vector of a positive-semidefinite matrix lies in its kernel,
so for every x∈X its pairing with the positive affine combination is
zero. Since every summand is nonpositive, x is orthogonal to every node
of L\{b}.

If two black deletions are affine, their supports cover L. Every external
white is then orthogonal to L. The external positive-definite block
cannot change any vertex's finite/ideal status. This handles19 cases.

If exactly one black deletion b is affine, external whites can contact
only b among the roots of L. Eliminating their positive-definite Gram
decreases only the diagonal at b:

\[
 G_L\longmapsto G_L-\tau e_be_b^t,\qquad\tau\ge0. \tag{2}
\]

This includes orthogonality to the internal white, if present. For every
other BLACK deletion j, its matrix H_j=G_{L\setminus\{j\}} is positive
definite, and the vertex condition is exactly

\[
 \tau\le\tau_j=1/(H_j^{-1})_{bb}.                 \tag{3}
\]

Only black deletions constrain actual vertices. In particular the
checker does not impose an extra vertex condition at an internal white
deletion. It evaluates all applicable thresholds, their minimum, and
every equality tie. At the minimum, the number of ideal vertices is
one plus the number of tied deletions; below it there is only the
original ideal vertex. In these14 cases there are at most two ties,
and the sole two-tie case is order-nine L21. This proves the claimed
bound of three for all33 cases with a black affine deletion.

The only remaining case is L25 with internal white9. Its nine BLACK
roots themselves form affine E8. It cannot be handled by (2), since
none of its black deletions is initially affine.

## Complete attachment analysis for the affine-E8 black block

Label the black affine-E8 roots0,…,8, all with squared norm2 and edges

\[
 0\! -\!1,\quad0\! -\!2\! -\!3,\quad
 0\! -\!4\! -\!5\! -\!6\! -\!7\! -\!8,
\]

all labeled3. Normalize every white root to squared norm2. A nonzero
contact magnitude d_j=−(u,b_j) of a white root with a black root is at
least1. At each black-deletion vertex, its contact vector obeys
d^tG^{-1}d≤2. The finite Gram inverse is entrywise nonnegative.

Deleting black8 gives the following inverse diagonal, in order0,…,7:

\[
 (30,8,14,4,20,12,6,2).                           \tag{4}
\]

Thus every white-to-black contact except those at7 and8 vanishes. A
contact at7 has magnitude at most1 and therefore exactly1. Deleting
black3 gives inverse principal block on nodes7,8

\[
 \begin{pmatrix}2&1\\1&1\end{pmatrix}.           \tag{5}
\]

Both contacts cannot be nonzero: even their minimum cost is5>2. A sole
contact at8 has magnitude at most sqrt2. The complete possibilities for
a white root touching B are consequently

| Profile | Nonzero black contact |
|---|---|
| A | node8, magnitude1, label3 |
| B | node8, magnitude sqrt2, label4 |
| C | node7, magnitude1, label3 |

This proof uses all Coxeter labels at once. It imposes no finite label
cutoff and does not rely on sampling or on enumerating white diagrams.
It also excludes parallel/dotted black contacts. A white root right to
all B belongs to the remaining set Z below.

Let U be all white roots directly touching B. It is nonempty because the
minimal diagram L includes the original touching white. At the vertex
obtained by deleting black3, the projection-cost Gram for the three
profiles is

\[
 K=\begin{pmatrix}
 1&\sqrt2&1\\
 \sqrt2&2&\sqrt2\\
 1&\sqrt2&2
 \end{pmatrix}.                                  \tag{6}
\]

Every entry is strictly positive. A B or C root has residual squared
norm zero. If there were any second direct white root, their residual
cross product would be strictly negative: its native product is
nonpositive and a positive entry of (6) is subtracted. A PSD matrix
cannot have a nonzero entry in a zero-diagonal row. Thus a B or C root
must be the only direct white root.

For two A roots the residual diagonal entries are1. Their residual
product is their native product minus1, hence at most−1. PSD requires it
to be at least−1. Thus the native contact is zero, and their residual
Gram is

\[
 H_A=\begin{pmatrix}1&-1\\-1&1\end{pmatrix}.     \tag{7}
\]

Three A roots would be pairwise orthogonal before projection and have
residual Gram2I−J, with a negative eigenvalue. Hence U has exactly one
of the four forms A, B, C, or two mutually orthogonal A roots.

## All remaining whites, including arbitrarily many indirect contacts

Put Z=W\U. By definition Z is right to B. Its Gram G_Z is positive
definite because W is. Write R for the Z-to-U product matrix. Eliminating
Z changes the direct-white Gram by subtracting

\[
 D=R^tG_Z^{-1}R.                                  \tag{8}
\]

The matrix D is PSD and entrywise nonnegative. The latter follows from
R≤0 and the nonnegative inverse of the positive-definite Stieltjes
matrix G_Z. No bound on |Z|, attachment lengths, number of components,
or contact labels is used.

For a single B or C, the residual scalar at black-deletion3 was zero;
the PSD condition after (8) forces D=0. For two A roots, (7) has positive
null vector(1,1). PSD of H_A−D forces
(1,1)D(1,1)^t=0. Entrywise nonnegativity forces every entry of D to be
zero. In both situations the remaining whites have no effect on any
vertex status. Equivalently, D=0 and positive definiteness of G_Z imply
R=0.

For a single A root the same vertex condition gives0≤D≤1. The certificate
computes all nine thresholds for decreasing the A-root diagonal:

| Deleted black | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| Threshold | 7/6 | 10/9 | 9/8 | 1 | 6/5 | 5/4 | 4/3 | 3/2 | 2 |

Only deletion3 can become ideal, and only at D=1. The complete status
checks for the four direct-white configurations are:

| Direct whites | Ideal black deletions after the allowed correction |
|---|---|
| A,0≤D<1 | none |
| A,D=1 | 3 |
| B,D=0 | 3 |
| C,D=0 | 3,8 |
| two A,D=0 | 3 |

In particular the exceptional case has at most two ideal vertices.
Together with the other33 cases, this proves the eight-simplex theorem.

## Payment bound and its remaining scope

The [finite-ended-port identity](fresh_finite_ended_simplex_ports_20261017.md)
partitions finite edges among maximal actual simplex blocks with at least
two finite vertices. For a block with f finite and i ideal vertices,
its twice-triangle payment divided by its number of finite edges is

\[
 2i+2(f-2)/3.
\]

For block dimension at most7, f+i≤8, so this is at most12. For dimension8,
the theorem above gives i≤3 and f+i=9, so it is at most26/3. For
dimension9, the
[frozen previous theorem](fresh_finite_edge_cusp_tradeoff_20261017.md)
gives exactly the possible counts(9,1),(8,2),(7,3), whose largest
coefficient is28/3. This proves (1), using
2E_ff+E_fi=n v_f at finite simple vertices.

The only formal equality type left by this relaxation is a seven-simplex
with(f,i)=(2,6). This package neither realizes nor excludes it. The
classification data for a minimal diagram must not be substituted for
the actual face counts without controlling all supporting whites.

For comparison with the
[exact cusp/truncation budget](fresh_finite_volume_simplex_truncation_20261017.md),
let Q denote the number of quadrilateral two-faces, and retain that
report's α_n,b_n and cusp-link edge counts E(L_v). The necessary
unpaid-face budget improves to

\[
 2(T_2+T_3)+Q\ge
 (b_n-6n)v_f+6E_{fi}+\alpha_n\sum_v E(L_v).        \tag{9}
\]

This follows from b_n v_f+α_nΣE(L_v)≤D_5(P) and
D_5(P)≤12E_ff+2(T_2+T_3)+Q. It leaves triangles with two or three ideal
vertices and the quadrilateral term uncontrolled. In particular, it
does not imply a dimension exclusion, including in dimensions23–26.
The primary literature distinctions in the previous report remain
unchanged.

## Exact checks and internal review

The [standalone checker](../scripts/fresh_finite_edge_cusp_tradeoff12_20261017.py)
and [saved data](fresh_finite_edge_cusp_tradeoff12_20261017.json) use no
project arithmetic, SymPy, or optimizer imports. They independently
reconstruct seven order9/10 diagrams from the frozen literal census,
verify their inertias and all positive affine kernels, and evaluate
all34 black choices. All14 rank-one cases retain their full threshold
lists and equality sets. The checker also reconstructs the literal
affine E8 Gram, all nine finite inverses, the complete attachment bounds,
the four direct-white configurations, their vertex inertias, and the
single-A thresholds. Finally it checks every integer(f,i) in the stated
coefficient domains.

Root derived the eight-simplex refinement and independently generated
and freshly replayed the full high-rank census. Frontier_literature's
hostile mathematical check found no gap and independently recomputed
the34 cases and affine-E8 matrices; its additional inverse-diagonal
argument above proves profile completeness directly. This fresh checker
records that independent small arithmetic. The geometric Schur and
positive-kernel steps were internally reviewed as written proofs.

Generation and fresh saved replay PASS; syntax, bound source hashes and
local report links PASS. The classification ceiling, completeness of the
frozen high-rank census, and geometric predecessor theorems remain
explicit inherited dependencies. No external peer review, formal
proof-assistant verification, or ambient dimension bound is asserted.
All predecessor files are unchanged.
