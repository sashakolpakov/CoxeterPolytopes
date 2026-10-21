# Joint quadrilateral patterns at A13 and D13 vertices

This census allows all four vertices of a quadrilateral to range over the [215 high-baseline ADE types](fresh_sqrt5_rank13_high_ade_endpoints_20261017.md), while requiring the initial vertex to be A13 or D13. It is complete for these two initial types and retains free opposite wall products.

| Initial forest | Tested attachment/bond choices | Literal finite-corner patterns |
|---|---:|---:|
| D13 | 11,780 | 7,219 |
| A13 | 20,850 | 11,806 |

Together these give9,895 colored patterns. Each saved literal map identifies both endpoint attachments, their destination types, the opposite destination type, the new–new finite product, and its exact colored isomorphism. Every root pair at both initial types retains at least one pattern with a genuinely coupled determinant factor. Thus rational-factor exclusion alone leaves all78 possible quadrilateral corners.

The checker uses the complete joint endpoint lists. For each pair of released roots, it tests every pair of endpoint attachments and both new–new products, zero and minus one. It then requires the fourth corner to belong to the joint set. For each colored pattern it independently computes the rational four-root Schur matrix, its bivariate determinant, and exact factorization. Rational univariate factors are excluded by the [rational-edge arithmetic gate](fresh_sqrt5_mixed_d4_d5_quad_factor_gate_20261017.md); coupled factors remain unresolved. No nonlinear univariate factor occurs in this data.

The fixed-edge graph forces both opposite products into the ground field in8,044 patterns, one product in1,633 patterns, and neither product separately in218 patterns. In the last group their squares are still in the field, and a cycle may additionally force their mutual product into it. A later parameter checker must preserve this distinction.

Using the original complete two-special-triangle corner sets gives the coarse per-vertex payments87/4 for D13 and22 for A13. The later [third-corner refinement](fresh_sqrt5_rank13_triangle_third_costs_20261017.md) lowers the triangle counts to5 and6 respectively, but with all78 quadrilateral pairs still present the corresponding payments remain83/4 and21. These are well above the target12. This is an unresolved arithmetic pattern list and a diagnosis of a weak estimate, not an impossibility result for refined joint incidence.

The [source](../scripts/fresh_sqrt5_rank13_irreducible_joint_quads_20261017.py) and [saved output](fresh_sqrt5_rank13_irreducible_joint_quads_20261017.json) pass generation and fresh complete saved replay:

```sh
python scripts/fresh_sqrt5_rank13_irreducible_joint_quads_20261017.py
python scripts/fresh_sqrt5_rank13_irreducible_joint_quads_20261017.py --verify
```

The [finite-anchor theorem](fresh_sqrt5_compact_quad_finite_anchors_20261017.md) supplies the next complete arithmetic reduction. It has not been applied inside this census. Neither the other213 selected initial types nor any unselected type is controlled here, and no dimension bound is asserted.
