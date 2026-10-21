# Exact metric gate for the first star assignment

The first mandatory-label survivor of star case 0 has **no hyperbolic metric completion**. This excludes one saved assignment, not the star incidence. The [exact checker](../scripts/fresh_h9_star_cube_metric_20261014.py) and its [certificate](fresh_h9_star_cube_metric_20261014.json) pass generation and saved replay.

The eight doubled-negative normals, at indices `[3,4,6,7,9,10,12,13]`, form four normalized B2 blocks. All 24 cross-block contacts are forced right: introducing label three gives a non-positive-definite actual-vertex Gram, and increasing that contact or any omitted contact cannot restore positive definiteness. The four own positive-leaf/cut-leaf pairs `(2,14),(5,15),(8,16),(11,17)` are forced right by the same test. Literal vertex indices and exact nonpositive LDL pivots are archived for all 28 claims.

The positive eight-dimensional block has inverse blocks

\[
\begin{pmatrix}2&\sqrt2\\\sqrt2&2\end{pmatrix}.
\]

Each of the eight leaves has contacts `−1/2,−1/2` with its own B2 block. For each leaf the incidence has an actual spherical support containing that leaf, one own B2 normal, and all six foreign B2 normals. Its Schur cost is at least

\[
\frac14+2\sum_{\text{foreign }j}\cos^2(\pi/m_j).
\]

The inverse entries and acute signs make every omitted cross term nonnegative. Two foreign nonright contacts already cost at least `5/4`; a single foreign label at least four costs at least `5/4`. Therefore there are exactly seven possible projection patterns for each leaf: no foreign contact, or one of six foreign contacts at label three. Every actual support used for this argument is explicitly checked.

Eliminating the positive block leaves a Lorentzian space of dimension at most two. A leaf's residual square is either `−√2/2` or `−(1+√2)/2`. Choose the final positive leaf `P3=11` and cut leaf `C3=17` as anchors. Their mutual original pairing is zero. In all 49 projection-pattern pairs their residual Gram determinant is strictly negative, so they are independent and no degenerate anchor branch is omitted.

The first cut leaf `C0=14` meets both anchors. Separate actual-vertex failures at label four force each of these two angles to label two or three; the argument excludes every larger integer label by monotonicity. Thus all possibilities for the three residual vectors are covered by

\[
49\times7\times2^2=1,372
\]

exact `Q(√2)` calculations. Every residual three-by-three determinant is nonzero: 150 are negative and 1,222 positive. Adding back the positive eight-dimensional block gives an eleven-normal principal Gram matrix of rank eleven, whereas an H9 realization has rank at most ten. No dotted distance or remaining central normal can repair this obstruction.

This is a finite exhaustive rank argument, not a numerical fit. It does not rely on the numerical contact-probe file. The remaining mandatory assignments of this same incidence are still available to the continuing label search.

The [separate review](fresh_h9_star_cube_metric_20261014_review.md) independently reconstructs all 30 domain witnesses, eight Schur supports, 49 anchor determinants, and 1,372 residual determinants using algebraic-field matrix arithmetic. Its generation and saved replay pass.

```sh
python scripts/fresh_h9_star_cube_metric_20261014.py --verify
```
