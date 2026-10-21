# Exact arithmetic parameters for joint A13/D13 quadrilaterals

The complete finite-anchor calculation reduces the possible all-selected quadrilateral corners at a D13 vertex from 78 to **35**, and those at an A13 vertex from 78 to **21**. Here every other corner may belong to any of the frozen 215 high-baseline rank-13 ADE forests. These are necessary local arithmetic possibilities, not assertions that the retained configurations occur in a compact polytope.

The [source](../scripts/fresh_sqrt5_rank13_irreducible_quad_parameters_20261017.py) and [saved exact output](fresh_sqrt5_rank13_irreducible_quad_parameters_20261017.json) bind the [complete finite-pattern census](fresh_sqrt5_rank13_irreducible_joint_quads_20261017.md), the [finite-anchor theorem](fresh_sqrt5_compact_quad_finite_anchors_20261017.md), its [independent review](fresh_sqrt5_compact_quad_finite_anchors_review_20261017.md), and the [refined triangle census](fresh_sqrt5_rank13_triangle_third_costs_20261017.md).

| Initial vertex | Retained literal maps | Possible quadrilateral corners | Refined triangle corners |
|---|---:|---:|---:|
| D13 | 95 | 35 | 5 |
| A13 | 53 | 21 | 6 |

There are 263,839 anchor tests across the 9,895 colored records. They yield 165 integral parameter pairs before the sign gates, counted separately in each record, and 136 surviving solutions in 112 colored records. The saved data identify every retained literal map, its free-product squares, any field-valued products or cycle product, and the opposite-pair and edge determinants.

The finite-anchor theorem says at least one of the two opposite full normal products is strictly greater than minus six. Both choices of the small member are tested. When a fixed 3-edge path makes that product a ground-field element, the exact eight-element direct-product list applies; otherwise its square is in the exact 59-element list. These reduced lists use the ADE rational-edge exclusion and must not be applied unchanged to non-ADE fixed contacts. The parameter calculation also distinguishes the case where neither product individually lies in the field but their mutual product belongs to an integral cycle. It solves the resulting linear or quadratic equation over Q(sqrt(5)), testing all exact square roots and both branches. An identically zero specialization is excluded by the separately proved positive-null-vector argument.

Both opposite residual pairs must have negative determinant. A positive definite opposite pair would give an additional finite vertex: eliminate the positive definite common-wall and opposite-wall block, and project a strict timelike center onto its common kernel. Its Stieltjes inverse makes the added normal coefficients nonnegative, so all remaining wall inequalities stay strict. This contradicts the stipulated opposite incidence. A singular pair instead supplies a nonzero positive null combination, feasible against every wall, contrary to compactness.

After these tests, the two old-edge residual 3-blocks must have negative determinant at the physical embedding and positive determinant at the conjugate embedding. Together with the already positive definite adjacent corner blocks and the zero full determinant, this gives physical/conjugate residual inertias (2,1,1)/(3,0,1). The field calculations use exact rational coefficients throughout.

Generation and a fresh complete saved replay passed:

```sh
python scripts/fresh_sqrt5_rank13_irreducible_quad_parameters_20261017.py
python scripts/fresh_sqrt5_rank13_irreducible_quad_parameters_20261017.py --verify
```

A separate root source/proof read found no issue in the square-root enumeration, both anchor axes, the direct/squared/mutual-product substitutions, or the determinant sign argument. That read is not an independent arithmetic replay of this parameter package. The inherited finite-anchor theorem has a separate independent arithmetic review.

The refined triangle corner sets are subsets of these quadrilateral sets. Giving a selected corner at most one-half for a two-selected triangle or one-quarter for an all-selected quadrilateral consequently gives vertex totals 10 for D13 and 27/4 for A13. The [separate incidence composition](fresh_sqrt5_rank13_irreducible_joint_incidence_20261017.md) states the resulting frequency inequality. Other selected types remain uncontrolled by this calculation; no dimension theorem is claimed.
