# Freeing all three internal labels in the sparse T8 gate

The sparse T8 obstruction remains valid when `(4,5)=4` is weakened to **any finite label at least four**. The exact hypotheses are now

```text
(1,2)>=7;  (4,5),(7,8),(12,13)>=4;
(1,14)=2;  (2,14)=(4,14)=(5,14)=3.
```

All other meeting labels and dotted distances remain free. The [primary checker](../scripts/fresh_h9_t8_sparse_search_20261014_all_blocks.py) and [certificate](fresh_h9_t8_sparse_search_20261014_all_blocks.json) pass generation and saved replay. The [independent Fraction review](fresh_h9_t8_all_blocks_20261014_review.md) separately reconstructs every changed inverse, projection and polynomial coefficient, checks the inherited cases and the continuous argument, and also passes both runs.

The underlying incidence and [sparse domain proof](fresh_h9_t8_sparse_frontier_20261014.md) are unchanged. Increasing the first middle label preserves every relaxed nonpositive actual-vertex witness and increases the inverse-block entries used for lower Schur-cost bounds. Thus the nine anchor-own patterns, finite third-anchor labels `{2,3}`, and at-most-one foreign contact per anchor remain complete domains.

Write `p` for the label on `(4,5)`, `q` for the label on `(7,8)`, and `t=cos(pi/m)^2` for `(1,2)`. Wall 14 has contacts of label three with root 2 and both roots 4 and 5. Its negative residual square is

`Q = 1/[4(1-t)] + 1/[2(1-cos(pi/p))] - 1`.

This replaces the former fixed-middle formula. For `p=4`, it is exactly the frozen sparse theorem.

For `p>=6`, every nonright anchor contact into `(4,5)` would cost at least one before the positive own-contact cost. Both anchors are therefore right to that block. Their residual Gram and mutual pairings with 14 are then unchanged from the appropriate frozen middle-label case, while Q only increases. Its residual determinant increases with Q because the two-anchor determinant numerator A is strictly negative. The whole unbounded `p>=6` range is consequently inherited, not sampled.

At `p=5`, put `phi=(1+sqrt(5))/2`. The exact quantities are

```text
Q at the first-label lower bound = 25/19 + phi;
foreign self cost              = (5+sqrt(5))/10;
foreign opposite-root overlap  = (5+3sqrt(5))/20;
foreign overlap with wall14    = (3+sqrt(5))/4.
```

The checker derives the last three from the normalized inverse two-node Gram. The first follows from `cos(pi/7)^2>81/100`, already certified in the base theorem.

There are only two changed finite blocks to check:

- For `q=4`, 256 cases have at least one foreign contact at 4 or 5. Without such a contact, only Q increases, so the old certificate applies.
- For `q=5`, only the 128 cases with one foreign contact at 4 or 5 and the other at 7 or 8 are additionally different. The other cases agree with the preceding `q=4` block or inherit the frozen `q=5` certificate with larger Q.

For `q>=6`, contacts into `(7,8)` are forbidden by the same cost-at-least-one argument, leaving a subset already covered. These splits exhaust every finite label; no finite cutoff is imposed.

The resulting 384 changed cases use the same six continuous own-contact boxes as the earlier theorem. Their 5,472 tensor Bernstein coefficients, over `Q(sqrt(2),sqrt(5))`, are strictly positive. The smallest is

`261/2000 - 31sqrt(5)/1000 > 0`.

The certificate records the 432 complementary inherited occurrences in the two 408-case tables as well. These are occurrences, not a claim of 432 different polynomials. Every changed case stores its variable order, exact intervals, power coefficients and Bernstein coefficients, including cancellations that alter its variable list.

The positive residual determinant and negative anchor determinant extend monotonically to every larger first internal label. As before, the positive eight-wall support leaves residual rank at most two in H9, contradicting the selected three-wall residual determinant. This proves the stated all-label family exclusion without choosing any dotted distance.

```sh
python scripts/fresh_h9_t8_sparse_search_20261014_all_blocks.py --verify
```

This theorem excludes the first assignment found by the new sparse continuation. It does not classify the whole case-8 incidence or supply a compact H9 example.
