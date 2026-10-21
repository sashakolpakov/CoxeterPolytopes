# H4 facet-context compensation still has an exact model obstruction

The enlarged model has the exact lower bound

**508953802651448/14513000101591 ≈ 35.06882099 > 35.**

This certifies an obstruction to the specified payment family. It is a lower bound obtained from a finite collection of vertex constraints, not a fully converged optimum. In particular, its numerical value should not be compared as an optimum with the weaker exact lower bounds from earlier models.

All inherited bonus coefficients are free. H4 compensation distinguishes actual H4 count 1–4, A1 absence or presence, whether the facet root lies inside the paid H4, and, for external marks, whether the marked component is another H4. Compensation at count 5 remains zero. This enlarges the [H4-count model](fresh_sqrt5_marked_h4_count_free_bonus_20261016.md) to 30,356 variables.

The model retains 1,979 inherited baseline-zero quadrilateral rows in 71,036 common contexts. It tests 163 selected triangle families in 139,392 marked contexts, including both possible types of external facet mark. Together these give 47,885 affine face profiles. A complete marked finite-forest dynamic program supplies vertex constraints. Every vertex in the exact dual is also recounted by literal root pairs.

The finite LP lower bound is approximately 35.07583996 after four solves and 5,033 distinct vertex cuts. The exact dual has 1,009 terms: 250 vertex terms, 45 triangle terms, and 714 quadrilateral terms. The distance-two discount is approximately zero and the adjacent-3 discount approximately 0.26895999. Both numerical upper-bound marginals are zero.

The largest vertex masses now occur at 4H4+D5 marked in D5 (about 33.85%) and A1+4H4+D4 marked in A1 (about 24.41%). Large triangle terms use an H4 component with an A-type white tail. This suggests distinguishing H4 compensation according to whether its rank-five non-H4 remainder contains D4 or D5. It is a targeted next model, not a geometric exclusion.

The [source](../scripts/fresh_sqrt5_marked_h4_facet_context_free_bonus_20261016.py) and [certificate](fresh_sqrt5_marked_h4_facet_context_free_bonus_20261016.json) passed generation and a fresh saved replay. The certificate contains the literal dual, variable schema, and numerical history.

**Scope:** restricted payment-model obstruction only. Newly affected triangles and paid-corner quadrilaterals do not yet have complete coverage. No new dimension bound or compact arithmetic example follows.
