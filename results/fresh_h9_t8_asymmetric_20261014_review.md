# Independent review of all T8 own-contact patterns

**PASS: independent generation and saved replay.** The [literal-domain checker](../scripts/fresh_h9_t8_asymmetric_20261014_review.py) and [certificate](fresh_h9_t8_asymmetric_20261014_review.json) verify the entire stated [primary family](fresh_h9_t8_asymmetric_20261014.md), including all nine nonempty support patterns for the two anchors. The first internal label ranges over every finite integer at least 7, the last over every finite integer at least 4, and all four anchor-own labels are unrestricted finite integers at least 2. Other mandatory labels retain their saved values. This is not an exclusion of the whole T8 incidence.

The separate [polynomial reconstruction](../scripts/fresh_h9_t8_block_polynomials_20261014_review.py) and [algebra certificate](fresh_h9_t8_block_polynomials_20261014_review.json) use Fraction arithmetic in `Q(sqrt2)`, with no producer arithmetic or SymPy imports. They reconstruct the foreign-contact projections from the two normalized B2 inverses, expand the scaled residual 3×3 and anchor 2×2 determinants directly, and verify their denominator-clearing identities. Every multivariate power coefficient and Bernstein coefficient is independently regenerated and matched to the primary certificate.

The six continuous boxes contain 408 cases and 11,736 strictly positive Bernstein coefficients. The smallest coefficient is exactly `7/50`. They prove the residual numerator positive and the anchor numerator negative throughout the boxes. Increasing the first internal cosine-square increases the residual determinant, so the inequalities extend to every allowed first label without a cutoff.

The application to the literal incidence is checked separately. The review verifies 24 actual-vertex cross-contact guards, all four retained-own-root anchor budgets, both third-anchor budgets, nine own-support patterns, and eighteen freshly reconstructed failures at third-contact label 4. Shared own roots force the anchor mutual contact right; otherwise its full finite cosine interval is retained. Both empty own supports are excluded by the saved minimal missing triples. For last labels 4, 5 or 6, a single own label 3 would make its missing triple positive semidefinite, which is incompatible with compactness. Thus those single contacts have label at least 4 and cannot have a foreign contact. This justifies the two low-label boxes.

No incidence symmetry is assumed. Renaming the own-block roots or exchanging the two anchor roles only relabels the residual algebra; every original support pattern has its own literal vertex witnesses.

The positive rank-eight support leaves a residual space of dimension two in a compact H9 realization. Three residual normals must therefore have zero determinant. The proved positive determinant gives the contradiction. No dotted distance was assigned or bounded in this argument.

```
python3 scripts/fresh_h9_t8_block_polynomials_20261014_review.py --verify
python3 scripts/fresh_h9_t8_asymmetric_20261014_review.py --verify
```

Both saved replays pass. The following middle-label continuation is outside this report's fixed middle-edge scope and requires its own certificate.
