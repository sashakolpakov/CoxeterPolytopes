# Independent review of the single label-4 third-contact family

**PASS.** The [independent checker](../scripts/fresh_h9_star_third_four_20261015_review.py) and [certificate](fresh_h9_star_third_four_20261015_review.json) verify the [primary theorem](fresh_h9_star_third_four_20261015.md) using separate Fraction arithmetic in `Q(sqrt(2),sqrt(5))`, sparse permutation determinants and exact polynomial division. No producer arithmetic or SymPy is imported. Definitions from earlier frozen independent checkers are reused; their source bindings are checked recursively.

The precise family has owner block `(3,4)` of finite label at least 5, the other three internal blocks of finite label at least 7, and third leaf 14 with own contacts `(2,4)`. Anchors 11 and 17 have single label-3 own contacts at either the same last-block root or opposite roots. All noncore contacts remain free. This is a compact H9 family exclusion, not an exclusion of the whole incidence.

The checker independently reconstructs all **108 positive-kernel witnesses and six actual-support budgets**, including the changed shared-anchor supports and their available-edge dependencies. The third's retained own cost is exactly one half. Adding any foreign nonright contact costs at least another half, so positive definiteness excludes the equality boundary as well. The anchor budgets and literal affine paths leave fourteen contact patterns per mode, reduced by naming the foreign anchor first to ten formulas. Both anchor roles are checked on actual faces; no global anchor interchange is assumed.

For the six patterns without a foreign contact in each mode, a separate six-variable determinant expansion verifies the identity

\[
\det R=vM-a(\delta_x^2+\delta_y^2),
\]

where the anchor residual diagonal `a`, the third residual diagonal `v`, and the anchor minor `M` are strictly negative. The contact guards give `delta_x delta_y=0`. The determinant is therefore strictly positive for every permitted owner and last label; no finite label cutoff is used.

For the four foreign-contact patterns in each mode, the owner label is forced to 5. The checker reconstructs every residual matrix and its cleared quadratic determinant directly. It independently obtains the **nine polynomial regions and all 24 positive and three negative Bernstein coefficients**. The single opposite-anchor exception is verified coefficientwise against the earlier independent formula. Its last-label-7 interval is negative, while the interval for every last label at least 8 is positive. The exact cubic bounds for cosine of pi/7 and the radical bound for cosine of pi/8 are checked again. The intervening continuous sign change is not discarded or treated as a valid positive region.

In every allowed case three residual vectors would have nonzero Gram determinant in the two-dimensional orthogonal complement of an actual positive-definite eight-root support. This supplies the contradiction. The eight literal owner/last-root reversals are checked against the complete saved vertex set and missing-face clutter. The covered 60,015-visit cursor is checked explicitly.

Generation and saved replay pass:

```sh
python scripts/fresh_h9_star_third_four_20261015_review.py --verify
```

No previous large polynomial bank or label search is repeated. This review makes no claim about a subsequent transfer binder, resumed search, or finite-volume noncompact realization.
