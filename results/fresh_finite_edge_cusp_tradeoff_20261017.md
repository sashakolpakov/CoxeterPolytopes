# A finite-edge triangle bound and the all-ideal literature boundary

Let P be a finite-sided, finite-volume acute Coxeter n-polytope. Write T_s
for its number of triangular two-faces with exactly s ideal vertices,
v_f for its finite vertices, E_ff for edges with two finite endpoints, and
E_fi for edges with one finite and one ideal endpoint. Then

\[
 2(T_0+T_1)\le14E_{ff}
 =7n v_f-7E_{fi}\le7n v_f.                       \tag{1}
\]

This improves the coefficient16 obtained from the dimension-nine
simplex-face bound alone. The improvement uses the three order-ten
minimal indefinite diagrams, together with a restriction on additional
supporting white roots. It does not identify an induced simplex face with
a Coxeter simplex. It does not provide an ambient dimension bound.

## Finite-edge blocks

The [finite-ended-port theorem](fresh_finite_ended_simplex_ports_20261017.md)
assigns every finite edge to a unique maximal actual simplex face S with
at least two finite vertices. Include isolated finite edges as
dimension-one blocks. If f=f(S) and i=i(S), it proves

\[
 \sum_S\binom f2=E_{ff},\qquad
 T_0+T_1=\sum_S\left[\binom f3+i\binom f2\right]. \tag{2}
\]

The [original-star argument](fresh_multi_ideal_triangle_scope_20261017.md)
proves dim S≤9 for every actual simplex face of a finite-volume Coxeter
polytope. In particular f+i≤10. The payment per finite edge in one block
is

\[
 c(f,i)=2i+\frac{2(f-2)}3.                         \tag{3}
\]

Without further information, (3) is at most16, with formal extremizer
(f,i)=(2,8). Below we prove that every actual nine-dimensional simplex
face instead has

\[
 (f,i)\in\{(9,1),(8,2),(7,3)\}.                 \tag{4}
\]

Its coefficient is consequently at most28/3. For dimensions at most8,
f+i≤9 and f≥2 give c(f,i)≤14, with equality possible in this numerical
relaxation only at (f,i)=(2,7). Summing (3) proves the first inequality in
(1). Finite vertices are simple, so counting their incident edges gives
2E_ff+E_fi=n v_f, proving the remaining statements.

The pair (2,7) is not asserted to be realizable, or impossible. Thus14
is a proved uniform coefficient, not a sharpness claim.

## The primary order-ten classification input

Chein's complete table has exactly three order-ten diagrams. Use node
labels0,…,9 and label3 on every unspecified edge:

| Name | Edges | Affine node deletions |
|---|---|---|
| A | Path0–1–…–7, plus2–8 and6–9 | 0,7,9 |
| B | Path0–1–…–8, plus6–9 | 0 |
| C | B, with edge0–1 relabeled4 | 0,8 |

The source is [Chein1969, appendix, printed p.16](https://www.numdam.org/item/M2AN_1969__3_3_3_0.pdf),
with the matrix definition and classification scope on printed p.4.
A [byte-bound source copy](../data/references/fresh_finite_edge_cusp_tradeoff_20261017_chein.pdf)
is included. The program reconstructs these literal graphs; it does not
independently enumerate all hyperbolic Coxeter diagrams.

All roots have squared norm2, except that for rational arithmetic in C
the first root is multiplied by sqrt2. That positive diagonal congruence
changes its squared norm to4 and the label4 product to−2; the Cartan
two-cycle remains2. Hence inertia and all zero/nonzero kernel statements
are unchanged. The program checks full inertia(9,1,0), every principal
deletion's inertia, and a strictly positive kernel at each affine
deletion. Every other deletion is positive definite. For an isolated
Coxeter simplex these give the counts in the table, but an arbitrary
face requires the additional argument below.

## Additional supporting roots cannot create eight ideal vertices

Let S be an actual nine-simplex face. Its original supporting normals
are W∪B, with |W|=n−9 and ten black roots B, one for each facet of S.
For dimension at least2, simplicity at positive-dimensional faces
ensures precisely this star description even when some vertices of S
are nonsimple ideal vertices in P. The archived argument extracts a
minimal indefinite diagram L containing all black roots. Since |L|≤10,
we have L=B. Thus B is one of A,B,C above, and W is positive definite.

For each black root b, the selected normals W∪(B\{b}) belong to the
corresponding actual vertex of S. Their Gram matrix is positive
semidefinite. It is positive definite exactly at a finite vertex. At an
ideal vertex these n selected normals lie in the null hyperplane of that
vertex, whose positive index is only n−1, and cannot have positive
definite Gram.

Suppose B\{b} is affine and let its positive kernel coefficients be h_j.
In a positive-semidefinite block matrix every vector with quadratic form
zero is in the full kernel. Therefore, for any w∈W,

\[
 0=\sum_{j\ne b}h_j(w,b_j).
\]

Every h_j is positive and every product(w,b_j) is nonpositive. Hence w
is orthogonal to every node of B\{b}. If two different black deletions
are affine, their unions cover B, so W is orthogonal to all of B. This
applies to diagrams A and C. Their finite/ideal counts remain(7,3) and
(8,2), respectively, regardless of the size or structure of W.

Diagram B has only one affine deletion, at node0. The same argument
forces all W-to-B contacts except those at0 to vanish. Eliminating the
positive-definite Gram G_W decreases only the diagonal of node0:

\[
 G_B\longmapsto G_B-\tau e_0e_0^t,\qquad
 \tau=c^tG_W^{-1}c\ge0.                           \tag{5}
\]

For a finite deletion j, let H_j=G_{B\setminus\{j\}}. By a rank-one
Schur complement, its modified matrix is positive semidefinite exactly
when

\[
 \tau\le\tau_j:=1/(H_j^{-1})_{00}.
\]

The certificate calculates these exact thresholds:

| Deleted j | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|
| τ_j | 2 | 3/2 | 4/3 | 5/4 | 6/5 | 7/6 | 9/8 | 1 | 10/9 |

All vertex conditions require τ≤1. For τ<1 only the original affine
deletion0 is ideal. At τ=1 precisely deletion8 becomes additionally
ideal. This proves (4). Each of its three pairs occurs already among
the three classified Coxeter simplices themselves, so the list is exact
for possible vertex counts of actual nine-simplex faces.

No corresponding assertion about all eight-simplex faces follows.
For example, a minimal order-ten diagram can contain nine black roots
and one internal white root. If the unique affine deletion is the white
root, none of the black vertex deletions is affine before adding further
whites. The orthogonality argument then gives no control. Reading the
face counts directly from the isolated ten-node diagram would be an
invalid shortcut.

## What the current primary literature check establishes

The following scopes are different:

| Scope | Primary result located |
|---|---|
| All vertices ideal, and P simple | n≤8: Felikson–Tumarkin, Theorem1 |
| All vertices ideal, and P right-angled | n≤6: Kolpakov, Corollary7 |
| General finite-volume right-angled, few finite vertices | Egorov's2026 preprint claims v_f≥4 in n7 and larger lower bounds in n8–12 |
| All vertices ideal, unrestricted Coxeter angles and cusp products | No additional dimension theorem covering this whole class was identified in this bounded source search |

The [Felikson–Tumarkin proof](https://arxiv.org/pdf/math/0502413) uses
simple ideal vertices: their affine diagrams are connected and have
exactly n nodes. It proves absence of triangular two-faces when n>5 and
restricts quadrilaterals, then applies the simple-polytope average-face
inequality. Product cusp links with several affine factors do not meet
that hypothesis. In particular, for n≥9 an all-ideal Coxeter chamber,
if one exists, must have at least one nonsimple ideal vertex.

The [Kolpakov proof, pp.1958–1959](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-p.pdf)
uses the minimum of24 facets for an ideal right-angled four-polytope and
the Nikulin–Khovanskii face-average inequality. Right-angled faces remain
right-angled, so every four-face supplies this lower bound. General
induced Coxeter faces do not have that closure property.

The [Egorov preprint, Theorem1.2 and Table2](https://arxiv.org/pdf/2603.29925)
claims finite-vertex lower bounds4,41,1312 in dimensions7,8,9. Its argument
uses parity for three-dimensional right-angled faces and orthogonal
gluing, then double counting. This is a recent preprint result, recorded
as a primary-source claim, not independently verified here and not
transferred to arbitrary Coxeter angles. The finite-vertex count table
is not a restriction on the number of finite vertices of general
Coxeter polytopes.

This was a targeted search, not an exhaustive current-best survey.
Statements about noncompact finite-volume examples in dimension19 or21
do not, by themselves, assert that those examples have all vertices
ideal. The words noncompact and all-ideal must not be substituted for
one another.

## The exact missing budget

Equation(1) controls precisely triangles with a finite edge. It vanishes
in the all-ideal case and gives no estimate for T_2,T_3 or quadrilaterals.
Using the notation of the
[finite-volume truncation report](fresh_finite_volume_simplex_truncation_20261017.md),
put Q=four-sided two-face count, d=n−1,
E(L_v)=edge count of the product-of-simplices cusp link, and

\[
 \alpha_n=5/\eta_n-1,\quad
 b_n=\alpha_n\binom n2,
\qquad
 \eta_n=\begin{cases}4(n-1)/(n-2)&n\text{ even},\\4n/(n-1)&n\text{ odd}.
 \end{cases}
\]

The archived face-average inequality and the elementary polygon identity
give

\[
 b_n v_f+\alpha_n\sum_v E(L_v)
 \le D_5(P)
 \le14E_{ff}+2(T_2+T_3)+Q.
\]

Consequently any such chamber necessarily satisfies

\[
 2(T_2+T_3)+Q
 \ge (b_n-7n)v_f+7E_{fi}
       +\alpha_n\sum_v E(L_v).                   \tag{6}
\]

The dropped terms from polygons of at least six sides are nonnegative
and could be retained for a stronger exact accounting. In the all-ideal
case (6) becomes a purely cusp/multi-ideal-triangle/quadrilateral budget.
At n23–26 the coefficient b_n−7n is negative, so this finite-edge estimate
alone is far from a dimension exclusion. A new bound for the remaining
two-/three-ideal triangles and quadrilaterals, compatible with general
product cusp links, is still required. The previously certified
multi-ideal pyramid counterexamples explain why the finite-ended
transitivity proof cannot simply be extended to all cusp ports.

## Certificate status

[Standalone Fraction source](../scripts/fresh_finite_edge_cusp_tradeoff_20261017.py)
and [saved exact data](fresh_finite_edge_cusp_tradeoff_20261017.json)
check the three literal diagrams, all30 vertex deletions, positive affine
kernels, the nine rank-one thresholds, and every integer(f,i) in the
stated coefficient domains. Generation and fresh saved replay PASS.
The source binds the primary PDF and the frozen geometric predecessor
proofs. Classification completeness and geometric arguments are written
proof inputs; no optimizer, finite-label extrapolation, external peer
review, formal proof-assistant verification, or new dimension ceiling is
claimed. Predecessor files were not edited.
