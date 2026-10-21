# Independent review of the two-high-label T8 family

**PASS.** The [primary certificate](fresh_h9_t8_metric_20261014_two_high.json) excludes compact H9 realizations of the new saved T8 assignment with both internal labels **`m,n >= 7`** and both final anchors having own contacts `(2,3)`. The complete T8 incidence and other mandatory assignments remain outside the conclusion.

The [review script](../scripts/fresh_h9_t8_metric_20261014_two_high_review.py) and [certificate](fresh_h9_t8_metric_20261014_two_high_review.json) use separate Fraction square-root-basis arithmetic and multivariate polynomial determinants. They import no primary producer or primary arithmetic. All 100 determinant identities are reconstructed directly, rather than obtaining the result by shifting the primary table from the preceding family.

The decreased anchor contacts require fresh domain checks, and these pass:

- The mutual-anchor right-contact witness uses the actual triangle **`[11,13,17]`**. Root 13 still has label-3 contacts to both anchors; the old triangle using root 12 would not justify the claim.
- Each anchor budget now **retains root 13 and deletes root 12**. Its known contact therefore still costs `1/4`. The first high-label block permits no nonright anchor contact, and the two foreign B2 blocks permit at most one contact, necessarily label 3. These are the five projection patterns used in the census.
- Both fresh label-4 upper-contact failures are independently verified by exact leading determinants. Their actual vertices contain neither whole parameter pair, so the matrices are independent of both high labels.
- All 24 cross-block witnesses are checked on their literal supports. The twelve relaxed affine paths remain valid when an internal label is increased. The two third-leaf supports also remain valid: the last high-label block has inverse diagonal greater than 5, strengthening the lower bound 2 used to rule out foreign contacts.

The anchor projection norm is now

`k = 1/[4(1-cos(pi/n)^2)]`,

because only one of its two own pairings is nonzero. It is not the equal-contact expression from the preceding family.

The review independently proves **`k > 5/4`** without numerical angle bounds. Put `c=cos(2pi/7)`. The seventh-root identity gives `8c^3+4c^2-4c-1=0`; `c>1/2`, and this polynomial is strictly increasing there. Its value at `3/5` is `-29/125`, so `c>3/5` and `cos(pi/7)^2>4/5`. The same strict bound holds for every finite `n>=7`. Substituting `c=2u-1` independently reproduces the primary cosine-square cubic and its rational checks.

For the first label, the weaker bound `t=cos(pi/m)^2>3/4` suffices. The positive eight-root block has determinant `(1-t)(1-u)/4`, where `u=cos(pi/n)^2<1`. The 25 anchor determinants are affine decreasing in `k` and already strictly negative at `k=5/4`.

For each of the 100 contact cases, the full bivariate determinant satisfies

`dD/dk = a1*q14 + (ac-bc)^2 > 0`,

with `a1<0` and `q14=-sqrt(2)/2-1/[4(1-t)]<0`. Its derivative with respect to `q14` equals the anchor determinant, which is negative; since `dq14/dt<0`, the determinant also increases with `t`. All cases are strictly positive at `(k,t)=(5/4,3/4)`, with exact minimum **`(1+sqrt(2))/8`**. Thus the PD-eight-block Schur decomposition forces rank eleven on the entire region, contradicting ambient dimension ten.

Generation and saved-result replay pass, as do source manifests, syntax and links. There is no finite label cutoff, no numerical acceptance predicate and no residual zero case. Earlier T8 review files remain unchanged.

```sh
python3 scripts/fresh_h9_t8_metric_20261014_two_high_review.py --verify
```
