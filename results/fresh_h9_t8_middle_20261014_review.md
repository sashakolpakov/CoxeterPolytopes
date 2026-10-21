# Independent review of the unrestricted middle-label extension

**PASS: independent generation and saved replay.** The [domain and coefficient checker](../scripts/fresh_h9_t8_middle_20261014_review.py) and [certificate](fresh_h9_t8_middle_20261014_review.json) extend the [independently checked own-pattern theorem](fresh_h9_t8_asymmetric_20261014_review.md) to every finite label at least 4 on the additional middle edge `(7,8)`. The first internal label remains at least 7; the last internal label is at least 4; all four anchor-own labels remain arbitrary finite labels at least 2. Other mandatory labels stay fixed.

Increasing the middle label from 4 preserves the old literal actual-vertex exclusions. Each relaxed Gram is a symmetric matrix with nonpositive off-diagonal entries; making an off-diagonal entry more negative cannot restore positive definiteness. The middle inverse Gram entries `1/(1-c²)` and `c/(1-c²)` both increase for `c>=0`, so all lower Schur-cost bounds on nonnegative contact magnitudes persist as well.

For a middle label at least 6, a nonright foreign contact has projection cost at least `1/[4(1-3/4)]=1`, before the positive own-contact cost. The inherited actual-vertex budget forbids it. The 152 old cases with no foreign contact in that middle block therefore suffice. Middle label 4 uses the complete old bank.

At middle label 5, exactly 256 cases change. The separate [quartic algebra checker](../scripts/fresh_h9_t8_middle_20261014_review_algebra.py) and [certificate](fresh_h9_t8_middle_20261014_review_algebra.json) reconstruct the normalized inverse Gram in `Q(sqrt2,sqrt5)`, then rebuild every residual determinant and multivariate Bernstein coefficient using Fraction arithmetic. It reuses only the earlier independent sparse-polynomial operations, replacing their coefficient operations in its own process; no producer arithmetic or SymPy is imported.

All 4,320 new Bernstein coefficients are strictly positive. Their minimum is exactly

`7/50 - 31 sqrt(5)/2000 > 0`.

Every coefficient matches the primary certificate, including cases where cancellation changes which variables occur. The inherited 152 cases are checked as the exact complementary subset of the 408-case base. The same rank contradiction and monotonic extension in the first label then apply.

```
python3 scripts/fresh_h9_t8_middle_20261014_review_algebra.py --verify
python3 scripts/fresh_h9_t8_middle_20261014_review.py --verify
```

Both saved replays pass. This excludes the stated seven-label family; it does not exclude the whole T8 incidence or establish a compact H9 example.
