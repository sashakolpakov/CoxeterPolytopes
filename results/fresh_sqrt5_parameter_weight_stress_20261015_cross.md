# Partial weights with a unique-A1–E7 cross contribution

Allowing one specific class of disconnected marked-pair weights lowers the exact rank-28 optimum from 79 to **230/3**, which is below the target 77. The same rule has rank-29 maximum **84**. The finite attachment frontier is still incomplete, so these are partial weight results and establish no dimension bound.

The [checker](../scripts/fresh_sqrt5_parameter_weight_stress_20261015_cross.py) and [exact certificate](fresh_sqrt5_parameter_weight_stress_20261015_cross.json) retain all 2,397 literal inputs of the [preceding diagnostic](fresh_sqrt5_parameter_weight_stress_20261015.md): 1,465 distinct parameter-marked graphs, comprising the frozen finite graph collection, subdivision and one-touch necessary patterns, and every singular zero-further-contact obligation in the complete dotted queue through dimension 29. Dotted ratios and feasible external spectator ranks are preserved.

All earlier baseline and connected-pair variables remain available. The only new weights apply when a marked pair joins an entire isolated A1 to a specified root of an entire E7 component, and the **whole vertex diagram contains exactly one A1**. The weight depends on the rooted E7 orbit. All other disconnected marked pairs have weight zero. Four rooted orbits occur among the input corners, adding four nonnegative variables to the previous 918.

The exact forest budget remains additive. If the whole vertex has exactly one A1, every E7 root forms precisely one pair with that root. Hence each E7 component receives the sum of its rooted-orbit bonuses, counted with their node multiplicities. There is no contribution in the zero-A1 or at-least-two-A1 contexts. The same exact finite-type dynamic program, tracking actual A1 count capped at two, therefore covers every forest; no pair-interaction relaxation is used.

The returned rule uses only one nonzero cross bonus, of weight **7/9**. It pays the isolated-A1/E7 corner of dotted state 107 and avoids forcing the no-A1 E8 total to 24. Its rank-28 maxima are all 230/3:

| Whole A1 count | One maximizing forest |
| --- | --- |
| zero | 3E8 + H4 |
| exactly one | A1 + 3D7 + E6 |
| at least two | 2A1 + 3E8 + I2(10) |

At rank 29 the corresponding maxima are 239/3, 84 and 236/3, attained respectively at D5 + 3E8, A1 + 4D7, and 2A1 + 3E8 + H3. All 92 allowed finite component types through rank 29 are checked with their actual bonuses.

There are 2,925 distinct face inequalities and 922 variables. Six cutting-plane iterations produce a rational primal and sparse dual with equal objective 230/3. Fraction arithmetic verifies every inequality, and complete forest separation verifies the bound beyond the finitely many cuts retained by the optimizer. Saved replay rebuilds the graph constraints and checks both certificates and forest budgets without invoking the optimizer.

All new weights are nonnegative, so the established triangle and persistent H4 baseline payments remain valid. The subdivision and one-touch inputs still lack arbitrary-extension sealing claims. Their inclusion is conservative for this diagnostic, not a replacement for complete finite-family coverage. The dotted queue's arithmetic and contact completeness are bound proof dependencies, not independently replayed by this optimizer.

Generation, saved replay, syntax, source-manifest and local-link checks pass. The old 79-obstruction package and every earlier frozen file remain unchanged.
