# Exact partial rank-13 rule using simplex packing

Adding 22 proved maximal-simplex incidence charges makes the enlarged flexible-star model exactly feasible below its target of 12. Every preceding face inequality and complete vertex budget is retained. This is a partial weight-rule certificate, not a dimension bound.

The [source](../scripts/fresh_sqrt5_rank13_simplex_packing_repair_20261017.py) combines the frozen [enlarged flexible-star model](fresh_sqrt5_rank13_full_flexible_simplex_lp_20261017.md) with the exact [simplex-port incidence theorem and column export](fresh_simplex_port_incidence_20261017.md). The [result](fresh_sqrt5_rank13_simplex_packing_repair_20261017.json) contains every coefficient and both rational certificates.

| Quantity | Value |
| --- | ---: |
| Whole rank-13 vertex forests | 1,870 |
| Retained simplex aggregate inequalities | 10,185 |
| Retained pair-orbit variables | 67,193 |
| Added packing-charge variables | 22 |
| Nonzero new face coefficients | 45 |
| Exact feasible maximum vertex budget | 239780393396 / 19999999999 |
| Exact dual lower bound | 50353882638360 / 4200000022609 |
| Exact slack below 12 | 219606592 / 19999999999 |

The feasible budget is approximately 11.98901967. The exact primal and dual bracket the numerical optimum; their equality is not asserted. The previous pair-only model had an exact lower bound strictly above 12.

## Why the additional charges are valid

Distinct maximal simplex faces through a vertex have edge-port sets intersecting in at most one port. Otherwise their shared pair determines a triangle with two different maximal-simplex extensions, contradicting the unique-extension theorem.

Choose a port subset `U` and a minimum block size `k` with `2k-|U| >= 2`. At most one maximal-simplex port block contained in `U` can have size at least `k`. Averaging this indicator over all diagram automorphisms keeps its total charge at a vertex at most one. The frozen export calculates its contribution to every old simplex aggregate, including rows outside the old positive dual support.

The successor introduces one nonnegative variable for each of the 22 exported averaged indicators. It does not remove or discount any old inequality. For each new column, it clears that column's rational denominators and multiplies its vertex cost by the same factor. These factors are only 1, 2, or 5. Thus the unscaled physical charge equals that factor times the stored LP variable; every new column still has unscaled vertex cost one.

All 78 original unordered root pairs remain in every vertex budget. The added terms are separate proved incidence costs, not additional root-pair multiplicities. The JSON distinguishes their variable kind and records the original and denominator-cleared face coefficients.

Eighteen of the 22 packing variables have positive weight in the recovered rule. They include the two diagnostic vertices `2 A1 + A2 + 3 H3` and `A1 + A5 + H3 + H4`. Their physical charge coefficients are respectively `20054901654/19999999999` and `9473493516/2857142857`.

## Validation and scope

The optimizer ran once with a 300-second limit and returned an optimum. Rational recovery verifies all 10,185 face aggregates and all 1,870 complete vertex budgets. The dual is also verified by exact rational column-load inequalities. A fresh saved reconstruction replay passed without invoking an optimizer.

The valid packing theorem concerns all maximal simplexes, but this computation tests only the inherited flexible-star face domain. Fixed-discriminant stars, other missing simplex stars, quadrilateral aggregates, and vertex forests absent from these inputs remain outside. No complete face coverage, common compact ambient realization, or dimension bound is claimed. The exact success demonstrates that global incidence information can overcome the previously certified obstruction to pair-only payments.
