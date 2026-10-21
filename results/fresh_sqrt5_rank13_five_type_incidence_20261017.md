# A joint frequency bound for five rank-13 vertex types

For a compact arithmetic hyperbolic Coxeter polytope of simplest type over Q(sqrt(5)) in dimension 13, count vertices by their entire diagrams. Let V_other count all types except the following five. Then

**39 V_D13 + 31 V_A13 + 9 V_(D5+2D4) + 22 V_(2D5+A3) + 6 V_(3D4+A1) < 276 V_other.**

In particular, these five types cannot account for every vertex of such a polytope. This is a simultaneous frequency restriction, not a dimension bound.

The calculation restricts the complete local lists to the actual five-type special set. The resulting upper payments at special vertices are:

| Whole vertex forest | Possible all-special quad pairs | Half-charge triangle pairs | Total upper payment |
|---|---:|---:|---:|
| D13 | 4 | 3 | 9/4 |
| A13 | 13 | 4 | 17/4 |
| D5+2D4 | 27 | 11 | 39/4 |
| 2D5+A3 | 17 | 9 | 13/2 |
| 3D4+A1 | 18 | 21 | 21/2 |

Some triangle pairs lie outside the quad sets. The source takes the literal maximum of the triangle and quad charge at every root pair, rather than assuming containment or adding the two charges.

No triangle has three special vertices, by the rational normalized Gram and conjugate-positivity argument. Give every nonspecial corner one unit and each of its H4 framing 5-bond pairs another unit. A two-special triangle is then paid by its H4 third corner, if present; otherwise give one-half at each special corner. An all-special quadrilateral receives one-quarter at each corner. Every other triangle and quadrilateral already receives at least two and one respectively from its nonspecial corners. Each nonspecial vertex has 78 root pairs and at most three H4 components, so its total is at most 81.

The strict h-vector bound is 12 V < 2 T + Q. Substituting the five upper payments and 81 for the remaining vertices gives the displayed inequality.

The [composition source](../scripts/fresh_sqrt5_rank13_five_type_incidence_20261017.py) and [saved output](fresh_sqrt5_rank13_five_type_incidence_20261017.json) bind the complete [irreducible](fresh_sqrt5_rank13_irreducible_quad_parameters_20261017.md) and [small-component](fresh_sqrt5_rank13_small_d_quad_parameters_20261017.md) parameter classifications, the [joint ADE triangle census](fresh_sqrt5_rank13_high_ade_triangles_20261017.md), and the [refined third-cost data](fresh_sqrt5_rank13_triangle_third_costs_20261017.md). Every retained triangle occurrence and literal quad map is identified by its source index. The source also saves all-215 upper bounds and the exact destination-type catalogue counts for comparison; these catalogue multiplicities are not vertex frequencies in a polytope.

Generation and a fresh saved replay passed:

```sh
python scripts/fresh_sqrt5_rank13_five_type_incidence_20261017.py
python scripts/fresh_sqrt5_rank13_five_type_incidence_20261017.py --verify
```

Allowing the full 215-type set at the other selected corners instead gives small-component upper payments 83/4, 19 and 24, all exceeding 12. Thus migration to other types is a concrete obstacle to using these simple uniform charges on a much larger set. The remaining costly ADE forests and the [24 identified non-ADE third forests](fresh_sqrt5_rank13_small_d_third_types_20261017.md) need further joint control. No compact ambient realizations are inferred from the retained local data.
