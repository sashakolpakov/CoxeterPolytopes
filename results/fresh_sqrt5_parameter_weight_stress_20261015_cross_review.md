# Separate replay of the unique-A1/E7 weight rule

The independent implementation reproduces the exact optimum **230/3** at rank 28 and maximum **84** for the same rule at rank 29. The [reviewer](../scripts/fresh_sqrt5_parameter_weight_stress_20261015_cross_review.py) and [result](fresh_sqrt5_parameter_weight_stress_20261015_cross_review.json) bind the frozen [primary certificate](fresh_sqrt5_parameter_weight_stress_20261015_cross.json).

The reviewer imports only the earlier separate 224-graph review engine, not any producer or optimizer. It reconstructs all 2,397 literal inputs, validates every origin map, and independently confirms 1,465 distinct parameter-marked graphs. Dotted ratios must agree before two graphs can be identified. The 306 connected pair orbits and four rooted E7 orbits are rebuilt from literal labelled graphs using all-root tree canonicalization, without interpreting the producer's canonical strings.

All 4,212 face cases and 2,925 distinct inequalities are reproduced. Feasible spectator cases are independently derived by splitting the residual rank into an exact number of A1 components and a remaining rank that is either zero or at least two, for ambient ranks 28 and 29. This checks the boundary cases rather than blindly imposing all three contexts.

The new cross term is checked separately. When the whole vertex contains exactly one A1, each E7 node makes exactly one disconnected marked pair with that root. The reviewer counts each rooted-node orbit independently and adds its multiplicity only to the context-one E7 total. The zero-A1 and at-least-two-A1 contexts have no cross term. Thus the forest accounting remains additive and is not a relaxation of intercomponent interactions.

Fraction arithmetic checks the complete primal and dual certificates and equality at 230/3. A reverse-type unbounded coin dynamic program, different from the producer's rank-first recurrence, checks every allowed finite forest using all 92 component types through rank 29. All three rank-28 maxima are 230/3; the rank-29 maxima are 239/3, 84 and 236/3. Every saved maximizing forest and each intermediate rank-28 maximum agree. Nonnegative bonuses preserve the previously reviewed triangle baseline.

Arithmetic contact-domain completeness, sealing and the input family classifications remain explicit dependencies. The separate [dotted-parameter review](fresh_sqrt5_dotted_parameters_20261015_review.md) and [438-state queue review](fresh_sqrt5_dotted_rank29_exhaustion_20261015_review.md) are bound in the result. This review does not convert the subdivision or one-touch families into arbitrary-extension closure theorems. The finite attachment frontier remains incomplete, so no dimension bound is claimed. The implementation is separate within the same construction lane; it is not presented as a different mathematical author.

Generation, saved replay, syntax, source-manifest and local-link checks pass. No frozen primary file was changed.
