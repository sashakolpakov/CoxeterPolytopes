# Freeing every inherited bonus still leaves a restricted-model obstruction

The sparse marked-facet model has an exact rational lower bound

**425633390406733/12138000036414 ≈ 35.06618793 > 35.**

Every one of the 30,169 inherited orbit/context bonus coefficients is a free nonnegative variable. No coefficient from the old dimension-22 primal solution is retained. The model also has 22 discount/exceptional-compensation variables and 33 normalized same-I2(5) context variables, for 30,224 variables in total.

After four numerical LP solves and 4,419 marked-vertex cuts, a finite relaxation has numerical lower bound 35.08886679. Exact rational recovery gives the lower bound above using 999 literal face/marked-forest terms. This is a certified lower bound for the specified model, not an asserted fully converged optimum of its complete vertex polyhedron.

The exact certificate retains all 1,979 inherited baseline-zero quadrilateral rows in 71,036 common spectator contexts. The selected 163 triangle families produce 122,416 marked-context cases; together these yield 41,133 distinct affine face profiles. The vertex separator covers every marked rank 21 finite forest using exact A1/I2(5)/H3/H4 counts and finite-component dynamic programming. Every vertex appearing in the exact dual is additionally recounted by literal individual root pairs, without the component-level marked-hub debit formula.

The two uniform discounts are bounded by 1 so that their baseline-one long-component corner payments stay nonnegative. These are restrictions of this payment family and its proof method, not geometric laws. They are inactive in the numerical obstructing LP: the distance-two discount is 0, the adjacent-3 discount is approximately 0.27221670, and both upper-bound dual marginals are 0. Small positive upper-bound multipliers in the saved exact dual arise from conservative rational recovery, not from a claim that the numerical caps were active.

About 60% of the recovered vertex-dual mass lies on 4H4+D5 with the facet root inside an H4. The largest triangle terms are connected-white extensions core51_tail4 and core46_tail5, and the largest quadrilateral term is prior row 64. This points to compensation contexts that distinguish multiple H4 components and the component containing the facet root.

The [source](../scripts/fresh_sqrt5_marked_free_bonus_proposal_20261016.py) generates and freshly replays the [certificate](fresh_sqrt5_marked_free_bonus_proposal_20261016.json). It records the full variable order, all dual terms and their literal witnesses, rational weak-duality verification, and numerical separation history. The earlier [restricted reallocation diagnostic](fresh_sqrt5_marked_reallocation_proposal_20261016.md) explains why freeing the old coefficients was necessary.

**Scope:** this excludes the specified algebraic payment model with the included obligations. It is neither a geometric nonexistence theorem nor a dimension bound. Only selected connected-white triangle families are imposed; newly affected triangles and paid-corner quadrilaterals still lack complete coverage. The included diagram obligations are used as explicit model constraints, not as a claim that every abstract selected diagram is globally realizable in a compact arithmetic polytope.
