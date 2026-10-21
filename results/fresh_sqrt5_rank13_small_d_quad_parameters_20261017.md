# Exact arithmetic quadrilateral parameters at the small D-component forests

Applying the complete finite-anchor theorem to the small-component census leaves the following necessary local possibilities when all other corners range over the full 215-type ADE set:

| Initial forest | Retained literal maps | Possible quad pairs | Refined triangle pairs | Corner-payment upper bound |
|---|---:|---:|---:|---:|
| D5+2D4 | 1,907 | 72 | 11 | 83/4 |
| 2D5+A3 | 762 | 63 | 13 | 19 |
| 3D4+A1 | 6,489 | 75 | 21 | 24 |

All three bounds exceed the dimension-13 h-vector threshold of 12. This is a limitation of the uniform corner payments on the enlarged set, not an impossibility theorem for better payments or for the vertex types.

The [source](../scripts/fresh_sqrt5_rank13_small_d_quad_parameters_20261017.py) and [saved exact parameters](fresh_sqrt5_rank13_small_d_quad_parameters_20261017.json) bind the [complete finite-pattern census](fresh_sqrt5_rank13_small_d_joint_quads_20261017.md), the [finite-anchor theorem](fresh_sqrt5_compact_quad_finite_anchors_20261017.md), its [independent review](fresh_sqrt5_compact_quad_finite_anchors_review_20261017.md), and the [triangle third-cost census](fresh_sqrt5_rank13_triangle_third_costs_20261017.md).

Across the 11,331 colored records, there are 496,221 anchor tests and 530 integral parameter pairs before the sign gates, counted separately per record. Exactly 504 solutions survive in 398 records. The method is the same complete both-axis calculation as the [irreducible parameter certificate](fresh_sqrt5_rank13_irreducible_quad_parameters_20261017.md): one opposite full product must exceed minus six; the correct eight-element direct or 59-element squared anchor list is tested on either axis. Products forced into the ground field, products known only through their squares, and a possible integral mutual cycle product are treated separately. All exact square roots and polynomial branches are tested. An explicit assertion checks that each specialized equation has degree at most two.

Identity specializations are excluded by the positive-null-vector compactness lemma. Both opposite residual determinants must be negative. The physical and conjugate edge determinant signs, together with the already spherical adjacent corners and full determinant zero, give residual inertias (2,1,1) and (3,0,1). The lists are necessary arithmetic local data; further integrality or global incidence conditions could remove retained cases.

Generation and a fresh complete saved replay passed:

```sh
python scripts/fresh_sqrt5_rank13_small_d_quad_parameters_20261017.py
python scripts/fresh_sqrt5_rank13_small_d_quad_parameters_20261017.py --verify
```

The large corner counts depend strongly on migration to other forest types. Restricting every corner to the five checked types D13, A13, D5+2D4, 2D5+A3 and 3D4+A1 reduces the three quad-pair counts to 27, 17 and 18. Requiring a single repeated initial forest reduces them further to 18, 15 and 9. The [separate five-type incidence certificate](fresh_sqrt5_rank13_five_type_incidence_20261017.md) uses the literal restricted maps and matching triangle occurrences. The non-ADE third forests are recorded in a [separate extraction](fresh_sqrt5_rank13_small_d_third_types_20261017.md). No dimension theorem is asserted here.
