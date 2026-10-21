# Star case 0: all own contacts π/3, arbitrary internal labels

**No compact hyperbolic 9-dimensional Coxeter polytope realizes this label family.** In the fixed 18-facet star incidence, suppose each of its eight leaves meets both roots of its own two-node block at angle π/3. The four internal block labels may be arbitrary finite integers; every other meeting label and all 15 dotted distances remain free.

The exact [producer](../scripts/fresh_h9_star_zero_20261014_own_three.py) and [certificate](fresh_h9_star_zero_20261014_own_three.json) prove this for the literal case0 incidence, with no identification of different incidence types. Let
\[
S=(3,4;6,7;9,10;12,13),\quad P=(2,5,8,11),\quad C=(14,15,16,17).
\]
The two leaves of each block are its corresponding entries of P and C. Root squares are normalized to one.

Each own triple is a minimal nonface and must be Lannér. Its labels are (3,3,m), so m≥4; m=3 is affine and m=2 spherical. If any internal label is at least6, the already reviewed [two-dense block theorem](fresh_h9_dotted_central_20261014_additional_gates.md) applies: choose that block and any appropriate other block with label≥4. The certificate supplies a literal placement for each of the four choices of high block, obtained from the [complete target table](fresh_h9_dotted_central_20261014_next_targets.md). Thus all remaining internal labels lie in {4,5}, leaving exactly16 choices. This is a proof for unbounded integer labels, not a cutoff assumption.

The thirty actual-vertex witnesses in the [first-star independent review](fresh_h9_star_cube_metric_20261014_review.md) remain valid when any internal label increases from4 to5. The normalized off-diagonal entries only decrease; the least eigenvalue of a symmetric matrix with nonpositive off-diagonal entries is nonincreasing under this change, by its nonnegative Perron eigenvector. The witnesses therefore force every cross-block pair and every Pi–Ci pair to be right. They also restrict the two contacts C0–P3 and C0–C3 to labels2 or3. Omitted contacts in those witnesses were already relaxed to right, so no additional labels are being assumed.

Consequently S is a positive direct sum of four dihedral blocks. For internal cosine c, its inverse block is
\[
\frac1{1-c^2}\begin{pmatrix}1&c\\c&1\end{pmatrix}.
\]
Each inverse diagonal is at least2. For each leaf the certificate identifies an actual spherical support containing that leaf, one own root and all six foreign roots. The retained own contact contributes1/4 to its projection cost. Two foreign nonright contacts would contribute at least1, while one foreign label≥4 also contributes at least1. Both violate positive definiteness. Hence a leaf has either no foreign contact or exactly one foreign π/3 contact: seven conservative projection patterns. The argument covers both internal labels4 and5; some patterns may be impossible, but retaining them cannot invalidate an exclusion.

Take the three additional normals P3=11, C3=17 and C0=14. For each of the16 internal-label choices, the checker exhausts their7³ projection choices and the two independent contact labels2/3. It reconstructs the three-dimensional Schur Gram over \(\mathbb Q(\sqrt2,\sqrt5)\). **All21,952 determinants are nonzero.** The784 two-anchor determinants are negative as an additional check; there are no zero residual cases requiring a further vertex gate. Since S has rank8, every resulting eleven-normal principal Gram has rank11. A hyperbolic9 realization has ambient dimension10, giving the contradiction.

This closes the necessary assignment found after21 visits by the [new resumable search](fresh_h9_star_zero_20261014_search_batch_001.json), whose internal labels were (4,4,4,5), as part of the whole all-own-π/3 family. It does not exclude the full star incidence: changed own contacts remain to be searched. The saved21-visit cursor is preserved; continuation uses a separately bound family gate rather than restarting the completed prefix.

Generation and the saved exact replay pass. The [independent review](fresh_h9_star_zero_20261014_own_three_review.md) reconstructs all21,952 determinants,784 anchor determinants,30 actual-vertex witnesses and128 budget checks with separate Fraction arithmetic and no producer imports. The first-state domain proofs and two-dense theorem remain explicitly bound dependencies.
