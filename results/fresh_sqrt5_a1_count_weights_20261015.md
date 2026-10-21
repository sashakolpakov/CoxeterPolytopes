# Exact partial weights distinguished by the number of A1 components

The specified 174-face system has exact rank-28 optimum **153/2**, below the face-count target 77. The same rule has rank-29 maximum **84**, equal to that rank's target. These are partial weight certificates: the remaining finite and dotted quadrilateral families are not covered, so neither number establishes a dimension bound.

The [checker](../scripts/fresh_sqrt5_a1_count_weights_20261015.py) and [certificate](fresh_sqrt5_a1_count_weights_20261015.json) preserve the distance-two/H4 baseline and add nonnegative bonuses to 72 long marked-pair orbits. Each bonus depends on whether the entire spherical vertex diagram has zero, one, or at least two isolated A1 components. This gives 216 variables. Every other bonus, including every disconnected marked pair, is zero.

The input is the same 174 literal face graphs as in the [Boolean-context diagnostic](fresh_sqrt5_rank28_weight_optimization_20261015.md): 112 old sealed cores, 46 sealed singular extensions, eleven sealed attachment states, and five earlier nonsealed states with no further contact. The [reviewed attachment closure](fresh_sqrt5_singular_attachment_closure_20261015_review.md) covers the touching cases of those five starts. The present rule is checked directly on all 174 graphs; it need not dominate the older bonus coefficients.

For a face with orthogonal external spherical forest, a corner's full A1 count is the number already in its core corner plus the number in that external forest, capped at two. The three external-count cases therefore cover all orthogonal spectators. They give 522 literal inequalities, or 274 distinct inequalities. Their extension validity uses the existing sealing and attachment proofs, rather than an assumed monotonicity under arbitrary component enlargement. Nonnegative bonuses leave all previously established baseline triangle, short-path and H4 payments valid.

The rank-28 results for the returned rule are:

| Full A1 count | Maximum vertex weight | One maximizing forest |
| --- | ---: | --- |
| zero | 76 | 3E8 + H4 |
| exactly one | 153/2 | A1 + 3D7 + E6 |
| at least two | 73 | 2A1 + D5 + 3D7 |

The second maximum is also attained at A1 + D6 + 3D7. This resolves the Boolean obstruction: the earlier forcing corner D6 + 2A1 no longer imposes the same bonus in a vertex with exactly one A1. The returned D6 totals are 12, 27/2 and 15 in the three contexts; its entire component and the spectator count, not merely its abstract pair, determine the weight.

The optimization uses all 89 allowed irreducible finite component types through rank 28. For each fixed context, an exact dynamic program tracks total rank and actual A1 count capped at two. It allows arbitrary component multiplicities and retains only the matching terminal count. Thus its maximum ranges over every allowed rank-28 finite forest, not only the diagrams appearing in the face inventory.

Three cutting-plane iterations suffice. Floating-point linear programming only proposes rational coefficients and multipliers. Fraction arithmetic checks every primal inequality, all dual inequalities, equality of primal and dual objective, and complete forest separation. The saved certificate contains the sparse exact dual attaining 153/2; this proves optimality for this particular baseline, orbit set and three-context model. Saved replay rebuilds the face constraints and component tables and verifies the certificate without rerunning the optimizer.

The same literal rule was evaluated, without reoptimization, on every rank-29 forest. The list is enlarged by the baseline-only components A29, B29 and D29, whose weights are respectively 55, 55 and 56. The resulting 92-type dynamic program gives:

| Full A1 count | Rank-29 maximum | One maximizing forest |
| --- | ---: | --- |
| zero | 79 | D5 + 3E8 |
| exactly one | 84 | A1 + 4D7 |
| at least two | 78 | 2A1 + D6 + 3D7 |

The relevant strict face-count thresholds are 77 at rank 28 and 84 at rank 29. Complete triangle and quadrilateral coverage with these vertex budgets would exclude those two dimensions. This package supplies only the stated 174-face payment system and its exact finite-forest budgets; the unresolved fronts remain separate work.

Generation, saved replay, syntax, source-manifest and local-link checks pass. No frozen dependency was changed.
