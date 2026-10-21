# A feasible partial rank-13 corner-orbit weighting model

The complete saved local constraint set admits an exact nonnegative rule with maximum vertex payment

**21538712201/3333333333 = approximately 6.461613661 < 12.**

This is a partial-model feasibility result, not a dimension-13 exclusion. The rule has already failed two omitted, actual compact local quadrilaterals; their exact rank-13 extensions and the repaired model are being handled in separate namespaces. The failure does not invalidate this certificate's explicitly stated input domain.

The [source](../scripts/fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017.py) and [saved domain and certificates](fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017.json) give an independent nonnegative variable to every unordered marked-root-pair orbit of each entire finite forest appearing in the input faces. All 78 root pairs, with their exact orbit multiplicities, contribute to that vertex's budget. No baseline coefficient or previous payment is fixed.

The domain has 630 vertex forests and 21,492 variables. It includes all 3,558 [old marked triangle constraints](fresh_sqrt5_rank13_triangle_third_costs_20261017.md), all 1,373 new [non-ADE triangle occurrences](fresh_sqrt5_rank13_compensating_endpoints_triangles_20261017.md), and all 510 retained colored quadrilaterals from the [irreducible](fresh_sqrt5_rank13_irreducible_quad_parameters_20261017.md) and [small-component](fresh_sqrt5_rank13_small_d_quad_parameters_20261017.md) parameter certificates. The 5,441 literal rows reduce to 4,324 distinct face inequalities. Each triangle must receive at least two units and each quadrilateral at least one.

The numerical minimum was 6.46161366006. The exact primal above and the exact dual lower bound 2326180916736/360000000019 bracket that numerical value. The primal is checked using rational arithmetic against every saved face and every full vertex budget. The dual is independently valid as a linear combination of the saved constraints: nonnegative face multipliers are charged to whole-vertex budgets using exact maximum orbit loads. The replay reconstructs the domain and checks both certificates without calling the optimizer.

```sh
python scripts/fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017.py
python scripts/fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017.py --verify
```

Generation and fresh saved replay passed. The first solve had a 300-second cap and finished optimally; no unbounded or duplicate cutting loop was used.

The missing coverage is substantial. Only 3,634 variables occur in a saved face constraint, and the recovered minimax rule uses only 839 nonzero variables. Compared with the old baseline, 5,844 orbits are discounted, including 4,310 orbits that appear in no saved constraint. All 630 represented forests have some discounted orbit. A complete finite-type census gives 6,789 rank-13 forests compatible with the Q(sqrt(5)) label restriction; 6,159 are absent here, including 2,174 with baseline strictly above 12. Simply retaining the old baseline elsewhere would therefore fail the target budget.

The exact certificates support only the saved partial model. Remaining triangles, mixed and non-ADE quadrilaterals, and faces at absent vertex types are outstanding. A local admissible Gram configuration is not an assertion of compact ambient realizability, and neither the primal nor its dual proves a dimension theorem.
