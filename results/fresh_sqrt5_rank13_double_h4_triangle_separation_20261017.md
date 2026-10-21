# Exact double-H4 triangle separation

Generation and a fresh saved reconstruction both pass, including the exact Gram, ray, center, and payment checks.

The [source](../scripts/fresh_sqrt5_rank13_double_h4_triangle_separation_20261017.py) evaluates all 1,395 occurrences in the [new complete selected-endpoint census](fresh_sqrt5_rank13_double_h4_endpoints_triangles_20261017.md) against the frozen [five-core spectator-repaired partial rule](fresh_sqrt5_rank13_joint_corner_spectator_repair_20261017.md). Of these, 954 have all corner forests represented in the old model, 478 are unpaid, and 441 touch an absent forest. The minimum known payment is zero. Exact [data and the literal geometric witness](fresh_sqrt5_rank13_double_h4_triangle_separation_20261017.json) are saved.

Witness 12 has corner forests `D5+2H4`, `2A1+H3+2H4`, and `A1+D4+2H4`. Its three existing variables 21185, 41383, and 12431 all have payment zero, while a triangle requires total payment two. Thus this is a failure at already represented vertex types, independent of how undefined weights are completed elsewhere.

The witness is checked as actual local geometry. Its 14 norm-two normals have products 0, -1, or -phi. The Gram has physical inertia (13,1), positive-definite conjugate, and 196 integral Cartan entries. The black indices are 1, 2, and 13; deleting any one gives a spherical rank-13 corner. The exact rays are minus the corresponding columns of the inverse Gram. Each is timelike, and all pairwise products are negative, placing them on the same future sheet.

The sum of the three rays is perturbed by `epsilon*G^-1(-1_white)`, with epsilon=1/16, to obtain an exact strict timelike center. On the common-white subspace every feasible vector is uniquely a nonnegative combination of the three displayed rays, because the black inequalities are exactly their coefficients. The resulting hyperbolic face is therefore precisely a compact triangle. Every wall occurs at a finite simple corner and is an actual local facet. A separate agent gave a narrow mathematical read of this argument, without duplicating the arithmetic replay.

This certifies a local admissible integral-reflection chamber and a failure of the saved payment vector's coverage. It does not assert a compact or finite-volume ambient polytope, invalidate the predecessor certificate on its stated rows, or exclude a repaired weighting rule. The entire 1,395-row batch is passed to a fresh successor rather than adding only this witness.
