# A joint D13/A13 vertex-frequency restriction

For a compact arithmetic hyperbolic Coxeter polytope of simplest type over Q(sqrt(5)) in dimension 13, let V_D13 and V_A13 count vertices with those entire diagrams, and let V_other count every remaining vertex. Then

**8 V_D13 + 21 V_A13 < 276 V_other.**

This excludes a polytope whose vertices all have these two types. It does not give a dimension bound or exclude the occurrence of either diagram.

The [exact composition source](../scripts/fresh_sqrt5_rank13_irreducible_joint_incidence_20261017.py) and [saved output](fresh_sqrt5_rank13_irreducible_joint_incidence_20261017.json) bind the [complete endpoint lists](fresh_sqrt5_rank13_high_ade_endpoints_20261017.md), the [triangle third-corner census](fresh_sqrt5_rank13_triangle_third_costs_20261017.md), and the [finite-anchor quadrilateral classification](fresh_sqrt5_rank13_irreducible_quad_parameters_20261017.md). Those lists allow a larger set of 215 ADE types at the other selected corners, so they remain valid upper lists when the special set is reduced to D13 and A13.

No triangle can have three special vertices: its normalized full Gram matrix would be rational and have a physical negative direction, contradicting positivity at a nonphysical embedding. Give every nonspecial corner one unit, with an additional unit on the 5-bond pair in each H4 component. A triangle with two special vertices is paid by its third vertex alone when that component is H4; otherwise each special corner receives another half. Give each corner of an all-special quadrilateral one-quarter. All remaining triangles and quadrilaterals already receive at least two and one respectively from nonspecial corners.

At a fixed vertex a root pair determines exactly one two-face. Its triangle and quadrilateral charges therefore contribute their maximum, not their sum. The five half-charge triangle pairs at D13 are among its 35 possible quadrilateral pairs, giving 5/2 + 30/4 = 10. For A13 the corresponding counts are six and 21, giving 3 + 15/4 = 27/4. A nonspecial vertex has 78 pairs and at most three H4 components, so its total is at most 81. The strict h-vector inequality in dimension 13 gives

12(V_D13 + V_A13 + V_other) < 2T + Q
≤ 10 V_D13 + (27/4) V_A13 + 81 V_other,

which rearranges to the stated inequality.

Generation and a fresh saved replay passed:

```sh
python scripts/fresh_sqrt5_rank13_irreducible_joint_incidence_20261017.py
python scripts/fresh_sqrt5_rank13_irreducible_joint_incidence_20261017.py --verify
```

The remaining high-cost forests, including D5+2D4, and the non-ADE third-vertex types still require a joint payment bound. The inequality above deliberately retains their full nonspecial cost.
