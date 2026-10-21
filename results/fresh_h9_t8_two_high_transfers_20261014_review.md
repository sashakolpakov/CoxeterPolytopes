# Independent literal transfers of the T8 negative-plane obstruction

**PASS.** The [separate checker](../scripts/fresh_h9_t8_two_high_transfers_20261014_review.py) and [certificate](fresh_h9_t8_two_high_transfers_20261014_review.json) verify all four [primary proof transfers](fresh_h9_t8_sparse_search_20261014_bridge_variants.json). Generation and saved replay pass.

There are two distinct sets of sufficient conditions for excluding a compact H9 realization of the literal T8 incidence:

```text
(1,2)>=7, (4,5)>=7, (7,8)>=4, (12,13)>=4; or
(1,2)>=7, (4,5)>=4, (7,8)>=7, (12,13)>=4.
```

Every other meeting label and divergent distance is free. The first line uses bridge 14 and the second bridge 15. Either anchor 11 or anchor 17 supplies the needed Lanner triple on block `(12,13)`.

The checker constructs these roles directly, without assuming any global incidence symmetry. For each of the four bridge/anchor choices, it independently locates the eight-root positive support, all 45 selected meeting pairs, 24 cross-block witnesses, four third-wall budgets, two anchor budgets, four affine stars and three missing triples in the literal source incidence. It also checks every primary support against its saved containing vertex and matches the two condition sets.

The scalar argument is reused from the [independently verified two-high-block theorem](fresh_h9_t8_sparse_search_20261014_two_high_blocks_review.md): the bridge has negative residual square less than minus three halves, the anchor less than minus one quarter, and their projections use disjoint blocks. The residual mutual pairing is zero or minus one half. Its determinant is strictly greater than one eighth, so these residuals span a forbidden negative-definite plane.

All supports needed for that argument are present under each new assignment. No additional angle is fixed, no polynomial bank is repeated and no dotted distance is evaluated. The whole T8 incidence remains open.

```sh
python3 scripts/fresh_h9_t8_two_high_transfers_20261014_review.py --verify
```
