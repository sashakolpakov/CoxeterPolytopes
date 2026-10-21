# Independent review: opposite-root T8 anchors with four unrestricted labels

**PASS.** The [independent checker](../scripts/fresh_h9_t8_continuation_20261014_opposite_review.py) verifies the [primary family exclusion](fresh_h9_t8_continuation_20261014_opposite.md), using separate Fraction arithmetic in `Q(√2,√5)` and polynomial operations. It imports neither the producer nor SymPy.

The domain is the specified opposite-root assignment of compact-H9 case 8: both internal labels are arbitrary finite integers at least 7, and the two anchors' nonright own labels are arbitrary finite integers at least 3. Other saved mandatory labels remain fixed. This is not an exclusion of the entire incidence.

The review checks all 24 literal cross-block guards, the changed opposite-root contact budgets, the mutual-anchor label bound, the third-root budgets and both finite-vertex upper-contact failures. In particular, an anchor with a foreign contact must retain own label 3; only foreign-free anchors have a continuously varying own contact.

For the minimal own contacts, all 400 cases are reconstructed by direct polynomial determinant expansion. Every residual numerator has three positive quadratic Bernstein coefficients on `9/10 ≤ c ≤ 1`, and every anchor-minor numerator has three negative coefficients. The independent roots-of-unity cubic argument proves `cos²(π/m) > 81/100` and `cos(π/n) > 9/10` for all finite `m,n ≥ 7`.

The general determinant and its derivative identities are separately verified as nine-variable polynomials with rational coefficients. For one varying own contact and one foreign contact, completing a square bounds half the derivative below by `7/32`. For two foreign-free anchors with nonright mutual contact, both derivatives are positive. With right mutual contact, a direct bound gives residual numerator greater than `2/5`. The anchor minor decreases as either allowed own contact increases, so monotonicity in the first internal label also persists.

Thus all four label ranges are handled without a finite cutoff. The [review certificate](fresh_h9_t8_continuation_20261014_opposite_review.json) binds the source manifests and records every checked case. Generation and saved-result replay pass:

```sh
python3 scripts/fresh_h9_t8_continuation_20261014_opposite_review.py --verify
```

The [next saved necessary assignment](fresh_h9_t8_continuation_20261014_next_case_08.json) has one anchor with own contacts `(2,3)` and the other with `(3,3)`. That asymmetric regime is outside this proof and remains unresolved.
