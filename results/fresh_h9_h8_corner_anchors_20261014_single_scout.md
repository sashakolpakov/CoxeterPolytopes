# Bounded rotation scout for the wall-8 anchored cap

The [scout](../scripts/fresh_h9_h8_corner_anchors_20261014_single_scout.py) and [saved numerical results](fresh_h9_h8_corner_anchors_20261014_single_scout.json) test the exact quadratic cap direction later frozen in the [scalar certificate](fresh_h9_h8_global_cap_gap_20261014.md). Generation and the deterministic saved replay pass.

The projected center of wall 8 is fixed to each of source vertices 0, 32, and 33 in turn. The remaining placement variables are the 28 rotations in its eight-dimensional tangent space. Eight deterministic starts per source corner use at most 400 least-squares function evaluations each. Physical coefficients are evaluated at 90 decimal digits before conversion to double precision.

The residual penalizes old-center containment violations against all eleven source facets, and uncovered source vertices against all twelve old balls. This is only a necessary relaxation: even a zero residual would still require all Coxeter angle relations and coverage of the whole source.

| Anchor source vertex | Best squared residual | Largest center violation at that trial | Largest uncovered-vertex violation |
| --- | ---: | ---: | ---: |
| 0 | 0.909269278 | 0.243442853 | 0.714525255 |
| 32 | 7.611845659 | 0.315754533 | 1.763346014 |
| 33 | 9.026550960 | 0.748874046 | 1.805468042 |

All three displayed best trials reached the evaluation limit. **This numerical no-find is not an exclusion theorem.** The separate [exact radial gate](fresh_h9_h8_corner_anchors_20261014_radial_gate.md) subsequently rules out every rotation for this cap by a simpler necessary geometric condition, independently of this scout.
