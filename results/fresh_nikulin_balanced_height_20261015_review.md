# Independent check: equal heights do not permit unit charge

The [primary obstruction](fresh_nikulin_balanced_height_20261015.md) passes an independent check by [determinant-cofactor arithmetic](../scripts/fresh_nikulin_balanced_height_20261015_review.py), with [saved exact data](fresh_nikulin_balanced_height_20261015_review.json). The reviewer uses SymPy principal determinants to compute inverse diagonals, rather than the producer's Fraction inversion, and imports no project code.

All 36 distinct E8 omission pairs and their 324 common unit-contact rows are reconstructed. There are 179 rows respecting both individual endpoint costs at most two. Every omission pair has such a row with total cost below 16/9. Since `p=q=e_k`, the two heights coincide, `M=1`, and the energy is zero. Thus every pair violates the proposed unit-charge local lemma. Identical omissions fail trivially with the common row on the omitted node itself.

The smallest necessary coefficient furnished by these witnesses is exactly **8/7**, uniquely for omissions `(1,6)`. Their common row at node 8 has endpoint costs 8/9 and 2/3. Ordinary E7 also has an equal-height counterexample: omissions `(3,6)`, common row at node 4, and costs 1/2 and 1 give total 3/2 below the target 14/9.

These are counterexamples to a local estimate under its stated cost and contact hypotheses. They are not assertions that the rows occur in a complete hyperbolic chamber. No eight-sevenths upper estimate has been proved, and no new dimension bound follows. The reviewed candidate remains 169. Any further improvement must use a stronger local estimate or additional global information.

```sh
python scripts/fresh_nikulin_balanced_height_20261015_review.py --verify
```
