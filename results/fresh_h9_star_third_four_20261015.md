# Single label-4 third contact with high anchor blocks

No compact hyperbolic 9-polytope realizes the following label family on the saved star incidence. The four root blocks are
\[
A=(3,4),\quad B=(6,7),\quad C=(9,10),\quad D=(12,13),
\]
with internal labels at least 5, 7, 7 and 7, respectively. Third leaf 14 has own contacts `(2,4)`. Anchors 11 and 17 each have one right and one label-3 contact to block D; their nonright contacts may use the same root or opposite roots. All other labels and dotted distances are unrestricted. Owner-root and last-root reversals are included explicitly.

This covers the first outside state at the frozen **60,015-visit cursor**, without asserting that the whole star incidence is excluded. No search continuation is part of this package.

The [checker](../scripts/fresh_h9_star_third_four_20261015.py) and [certificate](fresh_h9_star_third_four_20261015.json) bind the complete literal incidence, saved batch 014, and the frozen [opposite-single-third proof and independent review](fresh_h9_star_single_third_opposite_20261015_review.md). The new proof uses that package's guard recipes as input, but reconstructs their matrices, positive kernels, conditional label dependencies and actual supporting vertices at the weaker third label 4. It also rebuilds the shared-anchor versions. It does not import the earlier producer's arithmetic.

## Domains at the equality boundary

Normalize every root to square 1. Each of the 24 cross-block pairs is forced right by an actual relaxed affine path with labels 4–3–4. The eight roots consequently have a positive block-diagonal Gram matrix.

On the actual support obtained by deleting owner root 3 and retaining leaf 14, the leaf's own projection cost is now **exactly (1/2)**. A nonright contact to any foreign block costs at least (1/2), because every inverse diagonal of a root block with label at least 4 is at least 2. Their sum is at least 1, contradicting positive definiteness of that actual support. Equality is forbidden as well. Thus leaf 14 remains right to all foreign blocks.

Each anchor has an actual support with its right owner root deleted. Its retained own contact costs (1/4), so it has at most one foreign contact, necessarily label 3. The high-middle paths force this contact into block A. A foreign contact forces the owner internal label below 6, hence **exactly 5** under the present hypotheses.

The other literal affine paths, stars and forks require only third label at least 4. They give the same conservative contact rules in both anchor modes:

- mutual anchor label at most 5;
- at most one anchor has an owner-foreign contact;
- each third-anchor contact is right or label 3, and both cannot be nonright;
- a nonright mutual contact excludes both third-anchor and foreign contacts;
- a foreign contact at the third's nonright owner root forces both third-anchor contacts right.

For shared anchors some rules could be strengthened, but no such strengthening is used. The script checks **108 positive-kernel occurrences and six actual projection budgets** across the two modes. There are fourteen literal contact patterns per mode, reduced by naming the foreign anchor first to ten formulas. Every needed role is checked on actual faces; no global anchor-swap symmetry is assumed.

## The cases without a foreign contact

Put
\[
c=\cos(\pi/m_A),\quad t=\cos(\pi/m_D),\quad
q=1-c^2,\quad D=(1-t^2)^{-1},\quad \beta=1/\sqrt2.
\]
Let (mu) be the mutual anchor cosine, and (delta_x,delta_y) the third-anchor cosines. The latter have product zero. Both anchor residual diagonals are
\[
a=1-D/4<0,
\]
and the third residual diagonal is
\[
v=1-\frac{1}{2q}<0
\]
because (m_A\ge5). The determinant (M) of the two-anchor residual block is
\[
M=\begin{cases}
1-D/2-\mu D/2-\mu^2,&\text{shared anchors},\\
1-7D/16-\mu tD/2-\mu^2,&\text{opposite anchors}.
\end{cases}
\]
Both are strictly negative: (D>100/19>4) and (mu\ge0). The three-normal residual determinant is therefore
\[
\det R=vM-a(\delta_x^2+\delta_y^2)>0.
\]
This treats all six no-foreign patterns in each mode for every finite owner label at least 5 and every finite last label at least 7. It requires no coefficient bank or algebraicity restriction on the unspecified normals.

## The four foreign-contact formulas

A foreign contact pins (c=\cos(\pi/5)). Define
\[
F=\frac1{4q},\qquad
T=\begin{cases}\beta c/(2q),&\text{contact at root 3},\\
\beta/(2q),&\text{contact at root 4}.
\end{cases}
\]
Write (h=1) for shared anchors and (h=t) for opposite anchors. The residual matrix is
\[
R=\begin{pmatrix}
1-D/4-F&-hD/4&-\delta_x-T\\
-hD/4&1-D/4&-\delta_y\\
-\delta_x-T&-\delta_y&1-\beta^2/q
\end{pmatrix}.
\]
For contact at root 3 the permitted pairs (2(\delta_x,\delta_y)) are `(0,0)`, `(0,1)`, `(1,0)`; for contact at root 4 only `(0,0)` remains.

The checker forms (P=q(1-t^2)\det R) directly. All four shared-anchor polynomials and three opposite-anchor polynomials have strictly positive degree-two Bernstein coefficients on ([9/10,1]).

The remaining opposite-anchor case, with foreign contact at root 3 and ((\delta_x,\delta_y)=(0,1/2)), is
\[
P(t)=\frac{3(1+\sqrt5)}{32}t^2
-\frac{\sqrt{10}+\sqrt2}{64}t
-\frac{5+3\sqrt5}{64}.
\]
This is exactly the exceptional polynomial in the preceding theorem, with the owner and third label-4/5 roles exchanged. It is negative for last label 7 and positive for every last label at least 8. The exact enclosing intervals are
\[
9/10<\cos(\pi/7)<901/1000,\qquad
\cos(\pi/m)>23/25\quad(m\ge8).
\]
The source rechecks the cubic/Chebyshev argument and the quadratic-field comparison establishing these bounds. It keeps the integer split; positivity on the full continuous interval would be false.

Thus eight foreign formulas give **nine sign regions, with 24 strictly positive and three strictly negative exact Bernstein coefficients**. In every permitted case the residual determinant is nonzero. Three residual normals cannot have a nonzero Gram determinant in the two-dimensional complement of the actual PD8 support in Lorentz dimension ten. This proves the compact-family exclusion.

## Verification and scope

Primary generation and saved replay pass. The [independent review](fresh_h9_star_third_four_20261015_review.md) also passes generation and saved replay. Its separate `Fraction` implementation checks all 108 literal guards, six budgets, both domains, the no-foreign identity, all nine foreign regions and their 24 positive/three negative coefficients, eight orientations and the saved 60,015 cursor. No new predicates or resumed search are included in this package.

```sh
python scripts/fresh_h9_star_third_four_20261015.py --verify
```

The certificate checks eight literal owner/last-root orientations against every saved vertex and minimal nonface. It concerns the stated core family only, with no finite cutoff on any high label. The frozen 397-predicate search and its nine pending frames remain unchanged.
