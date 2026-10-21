# Enlarged partial quadrilateral weights with three A1-count contexts

The enlarged face collection still has exact rank-28 optimum **153/2**, below the strict face-count target 77. The same rule has rank-29 maximum **84**, equal to that rank's target. All added constraints fit without a finer spectator context. This remains a partial payment certificate: complete finite and dotted quadrilateral coverage has not been established, and no dimension bound follows yet.

The [checker](../scripts/fresh_sqrt5_enlarged_a1_weights_20261015.py) and [exact certificate](fresh_sqrt5_enlarged_a1_weights_20261015.json) extend the frozen [174-face rule](fresh_sqrt5_a1_count_weights_20261015.md). The distance-two/H4 baseline stays fixed, with nonnegative bonuses on long marked-pair orbits. Each bonus depends on whether the full spherical vertex diagram has zero, exactly one, or at least two isolated A1 components. Unlisted bonuses and disconnected-pair weights are zero.

The literal inputs are:

* the previous 174 graphs, comprising 169 sealed stars and five closed starts in their zero-further-contact case;
* twenty unpaid sealed stars from the immutable [610-state extraction](fresh_sqrt5_nonsingular_attachment_search_20261015_batch610.json), with its complete rejected contact domains;
* all 31 sealed stars in the immutable [two-chain export](fresh_sqrt5_two_chain_star_closure_20261015_export.json), including its three previously unpaid stars;
* the two two-chain closed starts in their zero-further-contact case.

The resulting 227 literal graphs reduce to **224 distinct marked graphs**. The three previously unpaid two-chain graphs are already among the twenty BFS graphs. Deduplication checks all edge labels, all four marked black nodes and their opposite-pair relation; the literal bijections are saved. No geometric equivalence is inferred merely from matching corner types.

Three possible external A1 counts give 672 face cases and 385 distinct inequalities. There are 94 long-pair orbits and 282 bonus variables. The same sealing and attachment proofs used by the input packages justify arbitrary orthogonal external spherical components; this optimization does not independently replay those arithmetic proofs. At a corner, the full count is its core count plus the common external count, capped at two.

Including the previously paid two-chain states is necessary when changing weights. A preliminary rule that constrained only the new unpaid graphs left children 64 and 65 unpaid in the zero-A1 context. Both are now explicit constraints, together with the other paid states and closed starts. In the final rule each child has face weight exactly one in all three contexts. All thirty supplemental checks pass, so their prior coverage is retained.

The exact maxima of the final rule are:

| Rank | Full A1 count | Maximum weight | One maximizing forest |
| --- | --- | ---: | --- |
| 28 | zero | 76 | 3E8 + H4 |
| 28 | exactly one | 153/2 | A1 + 3D7 + E6 |
| 28 | at least two | 73 | 2A1 + D5 + 3D7 |
| 29 | zero | 79 | D5 + 3E8 |
| 29 | exactly one | 84 | A1 + 4D7 |
| 29 | at least two | 78 | 2A1 + D6 + 3D7 |

All 89 allowed finite component types through rank 28 enter the exact forest dynamic program, with arbitrary multiplicities. Its state records both total rank and actual A1 count capped at two. At rank 29 the same rule, without reoptimization, is checked on 92 types; the added A29, B29 and D29 have only their baseline weights 55, 55 and 56.

The optimizer is used only to propose rational coefficients. Fraction arithmetic verifies every primal and dual inequality and equality of their objectives, while exact complete forest separation checks every rank-28 partition. Three cutting-plane iterations suffice. The sparse dual proves that 153/2 is optimal for this specified baseline, orbit set and three-count model. Saved replay reconstructs the full graph inventory and inequalities, rechecks the exact certificate and both forest budgets, and does not invoke the optimizer.

The unchanged baseline preserves the reviewed triangle and persistent H4 payments. The remaining attachment frontier and dotted cases are outside this certificate. Even though the numerical thresholds at both ranks are met, a dimension conclusion still requires complete face coverage.

Generation, saved replay, syntax, source-manifest and local-link checks pass. All imported frozen packages remain unchanged.
