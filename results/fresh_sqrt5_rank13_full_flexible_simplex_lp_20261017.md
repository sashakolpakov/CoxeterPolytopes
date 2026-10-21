# Enlarged flexible-star simplex payment model

The complete saved selected-endpoint triangle catalogue gives a new exact obstruction to the rank-13 target of 12, even after maximal-simplex aggregation and removal of all fixed-discriminant stars. This is a restricted payment-model result, not a dimension bound or a realization of the local stars in one compact polytope.

The [source](../scripts/fresh_sqrt5_rank13_full_flexible_simplex_lp_20261017.py) reconstructs its entire domain from the frozen [197,957-occurrence separation](fresh_sqrt5_rank13_full_triangle_simplex_separation_20261017.json) and the preceding [class-filtered model](fresh_sqrt5_rank13_discriminant_simplex_lp_20261017.json). The exact coefficients and certificates are in the [result](fresh_sqrt5_rank13_full_flexible_simplex_lp_20261017.json).

## Exact result

| Quantity | Value |
| --- | ---: |
| Maximal simplex stars and distinct aggregate inequalities | 10,185 |
| Whole rank-13 vertex forests | 1,870 |
| Retained old forests / newly added forests | 605 / 1,265 |
| Nonnegative corner-orbit variables | 67,193 |
| Exact lower bound | 25250122565220 / 2100000009539 |
| Exact feasible upper bound | 34353908270 / 2857142857 |

The lower bound is strictly greater than 12; the upper bound is strictly less than 14. The numerical optimum was approximately 12.0238678942. The two rational bounds bracket the optimum; no equality of the rational bounds is asserted.

Every vertex budget counts all 78 unordered root pairs. For a maximal simplex with `m` finite vertices, its aggregate corner demand is `2 * binomial(m, 3)`, and its literal corner occurrences are reconstructed from every positive-definite one-root deletion. No facet mark is removed. Each included star has an odd-sized connected component, so it remains in the conservative flexible class of the earlier determinant test. Flexibility is only a necessary compatibility treatment: it does not prove shared Hasse invariants, one common lattice, or simultaneous global realization.

The optimizer ran once with a 300-second limit. Exact primal recovery, an unminimized exact dual, and all literal vertex and face inequalities are checked using rational arithmetic. The saved replay reconstructs the full domain and verifies both certificates without invoking an optimizer.

## Obstruction support

The positive dual support contains 346 aggregate rows: 103 triangles, 102 tetrahedra, and 141 four-simplices. There are 203 positive whole-forest multipliers. The largest forest multipliers occur at:

- `A1 + A3 + H3 + H4 + I2(5)`;
- `2 A1 + A4 + H3 + H4`;
- `2 A1 + H3 + 2 H4`;
- `A1 + A5 + H3 + H4`.

The largest single face contribution is a four-simplex star with component sizes `[1, 3, 4, 6]`. Its finite vertex forests are `2 A1 + A4 + H3 + H4`, `2 A1 + H3 + 2 H4`, `A1 + A2 + 2 H3 + H4`, `A1 + A3 + H3 + H4 + I2(5)`, and `A1 + A5 + H3 + H4`. Thus the obstruction persists among exceptional-component compensation types, rather than being explained only by the separate A13 zero-payment counterexample.

## Scope and next implication

This model contains all flexible stars from the two-selected-endpoint catalogue, but omits fixed-discriminant stars, triangles outside that catalogue, other maximal simplex stars, quadrilateral aggregates, and vertex forests absent from these inputs. A rule satisfying these inequalities is therefore only a partial rule. Conversely, the exact dual rules out target 12 for a nonnegative whole-forest corner payment satisfying this particular set of aggregate demands.

The upper bound below 14 motivates a separate rank-14 investigation; it does not supply a rank-14 rule. That investigation must reconstruct rank-14 vertex budgets and face domains. In addition, the ambient quadratic form would have odd rank 15, so its determinant squareclass can change under a common totally positive scaling; the even-rank discriminant partition used here does not transfer unchanged.

Generation and a fresh saved replay passed. No new dimension bound follows from this package.
