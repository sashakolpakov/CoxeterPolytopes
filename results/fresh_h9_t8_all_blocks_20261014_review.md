# Independent review of the remaining T8 internal-label extension

**PASS.** The [primary all-block certificate](fresh_h9_t8_sparse_search_20261014_all_blocks.json) is independently verified for its stated sparse compact-H9 family. The condition `(4,5)=4` may be replaced by `(4,5)>=4`, with no finite upper cutoff. The other conditions remain

```text
(1,2)>=7, (7,8)>=4, (12,13)>=4;
(1,14)=2, (2,14)=(4,14)=(5,14)=3.
```

All other labels may vary subject to the literal incidence and compact Coxeter conditions. In particular all nine nonempty anchor-own support patterns remain covered. This is not an exclusion of the whole T8 incidence.

The [independent algebra](fresh_h9_t8_all_blocks_20261014_review_algebra.json) reconstructs both two-node inverse Grams, each foreign-contact projection and all third-wall projections with separate Fraction arithmetic in Q(sqrt(2),sqrt(5)). It imports only previously independent polynomial operations, with no producer arithmetic or SymPy. It expands the residual determinant and two-anchor determinant directly and verifies their symbolic identities.

At changed label five the third-wall cost gives

```text
Q = 25/19 + (1+sqrt(5))/2,
```

and a foreign contact to root 4 or 5 contributes `(3+sqrt(5))/4` to the magnitude of the third residual pairing. The exact changed bank has 256 cases when the other middle block has label four, and 128 additional mixed-foreign cases when it has label five. All 5,472 Bernstein coefficients in these 384 cases are strictly positive. Their minimum is `261/2000-31sqrt(5)/1000>0`.

The [application checker](fresh_h9_t8_all_blocks_20261014_review.json) compares every power coefficient, Bernstein coefficient, interval and projection cost with the primary artifact. It also verifies the exact partition into 384 changed and 432 inherited rows across the two displayed middle-label choices. These row counts do not constitute a discrete cutoff; the remaining labels follow the argument below.

The [independently reviewed sparse domains](fresh_h9_t8_sparse_frontier_20261014_review.md) continue to hold when label `(4,5)` increases. Their positive null-vector witnesses remain nonpositive under decreasing off-diagonal Gram entries; their nonnegative inverse-Gram cost lower bounds only increase. For either middle label at least six, a foreign nonright contact costs at least one before adding a positive retained own-contact cost, and is impossible at the corresponding actual finite vertex.

With no foreign contact into `(4,5)`, the anchor determinant and residual off-diagonal terms are unchanged. Only Q increases, and the independently certified anchor numerator satisfies `A<0`. The identity `dN/dQ=-A>0` therefore preserves strict positivity. The label-five cases without a contact to the other middle block reuse the new label-four bank literally; cases without a contact to `(4,5)` reuse the preceding middle-block bank at larger Q. These observations cover every remaining finite-label range and justify every inherited row.

The positive rank-eight support leaves at most two residual dimensions in ambient dimension ten. The positive three-residual determinant contradicts that rank. This uses the complete literal domain argument, not just a comparison of numeric coefficients.

```sh
python3 scripts/fresh_h9_t8_all_blocks_20261014_review_algebra.py --verify
python3 scripts/fresh_h9_t8_all_blocks_20261014_review.py --verify
```

Generation and saved replay pass for the [algebra script](../scripts/fresh_h9_t8_all_blocks_20261014_review_algebra.py) and [application script](../scripts/fresh_h9_t8_all_blocks_20261014_review.py). The previous independent domain and coefficient certificates remain frozen dependencies.
