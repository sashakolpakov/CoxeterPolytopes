# Independent review of the opposite-single core theorem

**PASS.** The [independent checker](../scripts/fresh_h9_star_opposite_core_20261014_review.py) and [saved certificate](fresh_h9_star_opposite_core_20261014_review.json) validate the new domains and all 3,612 replacement vertex witnesses in the [primary theorem](fresh_h9_star_opposite_core_20261014.md). No issue was found. Generation and saved replay passed.

The scope is the literal star case 0, with finite internal labels at least \((4,4,4,7)\), leaf 14 having own labels \((3,3)\), and leaves 11 and 17 having opposite single own labels \((2,3)\), \((3,2)\). Both last-root orientations are covered. Every noncore contact is unrestricted. This excludes that family, not the whole incidence, and performs no new label search.

## Fresh arithmetic and literal domains

The review imports only previously independent Fraction arithmetic and Laplace/Cramer helpers, never the producer's arithmetic. It checks the frozen source manifests recursively and reconstructs the literal core inside all 378 vertex incidences.

* All 24 relaxed cross-block paths have positive affine null vectors. Both third-contact upper paths have a positive vector of quadratic value \((2-\sqrt5)/2<0\). The remaining path has a positive affine null vector and rules out simultaneous third labels \((4,4)\). These 27 witnesses use only the stated core labels, lower bounds and trial contacts.
* A separate rational congruence verifies the mutual-anchor label-6 obstruction: its congruent matrix is the weighted path Laplacian with edge weights \(1/4,3/4\). Thus the mutual label is at most 5.
* All three retained-own eight-root budgets and the high-last triangle are actual faces. Own cost \(1/4\), foreign inverse diagonal at least 2, and nonnegative inverse entries permit at most one foreign contact, necessarily label 3. If a foreign internal label is at least 6 its inverse diagonal is at least 4, which rules that contact out.

The local permutation \((11\ 17)(12\ 13)\) preserves all nine maximal induced core faces and the fixed core data. Consequently one can orient the third-contact domain as \(\{2,3\}\times\{2,3,4\}\). The review also explicitly maps every selected fallback support to another actual face. This is a local face-list symmetry; no full-incidence automorphism is assumed. Separately, the last-root swap \((12\ 13)\) is verified on every full vertex and minimal nonface to cover the other own-contact orientation.

## High-first continuous argument

For first internal label at least 6, the residuals of leaf 14 and the sum of the two anchors have squared norms \(-(K+H)\), \(-(L+F)\), with \(K>2,L>1,H,F\ge0\). The own-block inverse formulas give these strict bounds; all middle-block contributions are positive-semidefinite projection costs. Their mutual pairing has magnitude at most \(\sqrt2+\sqrt{HF}\). The checker independently expands

\[
(2+H)(1+F)-(\sqrt2+\sqrt{HF})^2
=(\sqrt{2F}-\sqrt H)^2.
\]

The actual determinant is strictly greater than this nonnegative quantity, including when the square vanishes. Its two negative diagonal entries therefore give a negative-definite plane, impossible in a Lorentz space of index one. This part has no finite label cutoff.

## Inherited rank bank and replacement failures

The 47,040 polynomial cases are inherited through the [frozen independent opposite-high review](fresh_h9_star_opposite_high_20261014_review.md); their coefficients are **not recomputed here**. The new checker verifies every literal parameter key, the bound source hashes, and the saved Bernstein signs: 43,412 coarse positive cases, 16 refined positive cases, and 3,612 cases requiring vertex contradictions.

For every one of those 3,612 cases, the review reconstructs its relaxed core matrix and actual placement. It independently computes leading determinants by Laplace expansion, then uses Cramer's rule to construct a nonnegative vector with nonpositive quadratic value. All saved failure pivots agree. There are 240 distinct selected full matrices; the producer tried 262 matrix types. Because the witness vectors are nonnegative, increasing internal labels or adding further acute contacts preserves the contradiction. Thus no formerly required noncore label enters these replacement failures.

The result is an independent check of the new domain reduction, all new matrix witnesses and their monotone transfer. It is not a second implementation of the inherited 47,040 polynomial calculations or a new exhaustive incidence search.

```sh
python scripts/fresh_h9_star_opposite_core_20261014_review.py --verify
```
