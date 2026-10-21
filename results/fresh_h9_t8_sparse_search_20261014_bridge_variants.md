# Two literal bridge placements in T8

The [two-high-block negative-plane theorem](fresh_h9_t8_sparse_search_20261014_two_high_blocks.md) applies in two placements. In each row, all displayed labels are finite integers and every other label and dotted distance is free:

| High blocks, labels at least seven | Other blocks, labels at least four | Joining wall | Independent anchor |
|---|---|---|---|
| `(1,2)` and `(4,5)` | `(7,8)` and `(12,13)` | 14 | 11 or 17 |
| `(1,2)` and `(7,8)` | `(4,5)` and `(12,13)` | 15 | 11 or 17 |

Neither placement admits a compact H9 Coxeter realization of the literal case-8 incidence. This does not exclude the whole incidence.

The [checker](../scripts/fresh_h9_t8_sparse_search_20261014_bridge_variants.py) and [certificate](fresh_h9_t8_sparse_search_20261014_bridge_variants.json) pass generation and saved replay. They bind the primary theorem and its [independent arithmetic and proof review](fresh_h9_t8_sparse_search_20261014_two_high_blocks_review.md).

The proof transfer does not assert a global incidence symmetry. For each of four local renamings, the checker rechecks all 45 meeting pairs, 24 cross-block null-witness supports, two anchor budgets, four joining-wall budgets, four upper-angle stars, the positive eight-wall support and three missing triples in the actual saved incidence. The Gram matrices and scalar inequalities are identical after these variable renamings. These are all the incidence premises of the reviewed proof; no outside wall label occurs in its equations.

For the second row, the joining wall's missing triples are `{1,2,15}` and `{7,8,15}`. The anchor's missing triple is `{11,12,13}` or `{12,13,17}`. Thus the anchor and joining wall again have disjoint support projections after the checked Schur-cost guards. Their residuals span a negative-definite plane, contradicting Lorentz negative index one.

```sh
python scripts/fresh_h9_t8_sparse_search_20261014_bridge_variants.py --verify
```

The new computation checks literal applicability of a previously independently verified theorem. It does not independently repeat its field arithmetic or enumerate new label assignments.

The [separate transfer review](fresh_h9_t8_two_high_transfers_20261014_review.md) independently rebuilds all four placements from the saved actual vertices and checks every mapped premise. Its generation and saved replay pass. The final continuation binds both this review and the original negative-plane review.
