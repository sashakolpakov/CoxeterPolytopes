# Independent review of the forest-sensitive partial weight rule

**PASS.** The [independent checker](../scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py) and [saved result](fresh_sqrt5_forest_sensitive_weights_20261015_review.json) reproduce the new [partial rule](fresh_sqrt5_forest_sensitive_weights_20261015.md), with no producer or project arithmetic imports. Generation and saved replay pass. The source and proof review found no issue.

The checker reads literal labelled graphs and rational bonus values from the frozen certificate, then uses a separate center-rooted marked-tree invariant. It independently reconstructs the two full-weight and two half-weight D7 orbits, all 112 sealed quadrilateral cores, and both possible external-A1 states. All 224 face sums are at least one. It also reconstructs every pair weight in all 92 allowed irreducible finite component types through rank 29 and the 76 minimal triangular-face markings.

An independent unbounded coin-change calculation tracks both rank and actual A1 occurrence, rather than reserving an A1 in advance. It reproduces every saved dynamic-program value and the sharp rank-29 maxima:

| Full vertex forest | Maximum | Attaining forest |
|---|---:|---|
| No A1 component | \(159/2\) | \(D_6+E_7+E_8+E_8\) |
| Contains A1 | \(84\) | \(A_1+4D_7\) |

The finite-type list is complete under the inherited arithmetic finite-label restriction: A, B, D, E, F, H, and the remaining rank-two types with labels 5, 6 and 10. Rank-two labels 3 and 4 already occur as A2 and B2; a label-2 pair is disconnected. The rank computation therefore includes arbitrary multiplicities and arbitrary finite spherical forests of total rank 29, not only forests seen in the local census.

The geometric step is conditional sealing, not arbitrary extension monotonicity. For each selected core, the independently reviewed one-white census excludes every nonempty unpaid contact. If an external white component touched the core, one of its nodes would have a nonempty direct contact with the core and supply a forbidden principal extension, unless that extension already created a baseline payment. Thus, in the still-unpaid case, the entire external white forest is orthogonal to the core. It is the same white forest at each of the four corners. Its isolated-A1 flag is consequently the only additional datum needed by this new rule. A displayed isolated A1 cannot be silently joined to an external component in that case.

The rule remains at least the original distance-two/H4 baseline at every marked pair. The reviewed triangular-face extension proof therefore still applies: a retained short path stays short in a finite forest, while a connected H4 cannot acquire further nonorthogonal spherical nodes. Existing short-path or H4 quadrilateral payments likewise survive. No monotonicity of the extra D7 bonuses under unrestricted component growth is assumed.

This review **inherits**, rather than repeats, the complete one-white contact enumeration, arithmetic alphabet, and triangle extension theorem from their frozen independent reviews. It newly verifies the forest-sensitive accounting, all literal selected faces, all finite-type weights, and the rank budget. The 130 nonsingular finite unpaid classes and 24 symbolic dotted unpaid classes remain unresolved. The valid upper vertex weight 84 gives no dimension exclusion until those remaining faces and their arbitrary external-white extensions are covered by the same rule.
