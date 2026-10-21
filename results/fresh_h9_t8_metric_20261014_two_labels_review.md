# Independent review of the two-label T8 extension

**PASS.** The [extension certificate](fresh_h9_t8_metric_20261014_two_labels.json) excludes the saved T8 compact-H9 metric family when the two internal labels are allowed to vary independently over **all finite `m >= 7` at `(1,2)` and `n >= 4` at `(12,13)`**. Other mandatory assignments and the whole incidence remain outside this conclusion.

The [separate checker](../scripts/fresh_h9_t8_metric_20261014_two_labels_review.py) and [review certificate](fresh_h9_t8_metric_20261014_two_labels_review.json) inherit the frozen [independent one-label domain review](fresh_h9_t8_metric_20261014_review.md), verify its literal support persistence, and reconstruct all **100 full bivariate determinant polynomials** using Fraction square-root-basis arithmetic and a separate permutation determinant. No primary producer or primary arithmetic is imported.

The changed inverse matrices do not invalidate the contact domains:

- Every nonpositive actual-vertex witness persists when the internal label `(12,13)` increases from 4: the corresponding nonnegative cosine weight only increases. This includes the cross-block and upper-contact guards.
- Each of the four anchor budget supports omits one of roots `12,13`. Consequently the changed edge is absent from these matrices and the original five-pattern anchor budgets are unchanged.
- The two third-leaf budget supports retain the entire varied block, but their required bound only improves. For `s=cos(pi/n)` in `[sqrt(2)/2,1)`, the inverse diagonal is `1/(1-s^2) >= 2`, with nonnegative off-diagonal entries. The known contacts already cost `1/2`; any additional foreign nonright contact costs at least `1/2`, contradicting strict positive definiteness.
- The mutual-anchor affine-triangle witness uses root 12 and does not include root 13, so it is unchanged.

The positive eight-root block has determinant `(1-t)(1-s^2)/4 > 0`, where `t=cos(pi/m)^2`. The equal anchor contacts to the varied block have projection norm

`k = (1+s)/[2(1-s^2)] = 1/[2(1-s)] >= k0 = 1+sqrt(2)/2`.

For each pair of foreign patterns, the anchor determinant is affine in `k`. Its derivative `a1` is strictly negative and its value at `k0` is negative, so it stays negative for all `k >= k0`. The separate coefficient comparison establishes, as a polynomial identity rather than a finite sample,

`dD/dk = a1*q14 + (ac-bc)^2 > 0`,

because `a1 < 0` and `q14 = -sqrt(2)/2 - 1/[4(1-t)] < 0`. Also `dD/dq14` equals the negative-valued anchor determinant, while `dq14/dt < 0`; hence `dD/dt > 0`.

Every determinant therefore has its lower boundary value at `(k,t)=(k0,3/4)`. All 100 values agree exactly with the already independently certified one-label values, whose minimum is `3/4+5sqrt(2)/8 > 0`. The PD-eight-block Schur decomposition forces rank eleven throughout the full two-parameter region. There is no finite angle cutoff and no residual zero case.

Generation and saved-result replay pass. Source manifests, syntax, and report links pass. The older one-label and star-family files are unchanged.

```sh
python3 scripts/fresh_h9_t8_metric_20261014_two_labels_review.py --verify
```
