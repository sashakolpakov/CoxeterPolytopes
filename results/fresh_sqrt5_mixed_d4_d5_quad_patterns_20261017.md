# Mixed D4/D5 quadrilateral corner patterns

This is a complete finite-corner census for quadrilaterals whose four vertex forests equal one of the four forests below. The opposite wall products remain free. No arithmetic solution or local quadrilateral is asserted merely because its finite corners pass.

| Vertex forest | Dimension | Literal patterns | Colored patterns | Potential corners per vertex |
|---|---:|---:|---:|---:|
| 3D4 | 12 | 66 | 5 | 66 |
| D5+2D4 | 13 | 210 | 38 | 78 |
| 2D5+D4 | 14 | 354 | 68 | 91 |
| 3D5+A2 | 17 | 607 | 60 | 136 |

For each pair of released roots, the checker takes every pair of the complete same-type endpoint attachments in the [endpoint census](fresh_sqrt5_mixed_d4_d5_migration_20261017.md). It tests both possible new–new finite products, zero and minus one, and requires the fourth vertex to have the same forest. Thus copied attachments, opposition changes, and migration between components are all retained. Every deduplication has a literal black-opposite/white colored isomorphism map.

The black order is cyclic, with free products `u=G(0,2)` and `v=G(1,3)`. Label zero in these two positions denotes a free variable, not an assertion that the walls are dotted. Eliminating the common positive-definite white support gives a rational polynomial four-root Schur matrix. Its determinant must vanish. The saved output contains that full matrix and determinant for every colored pattern, checks positive definiteness at all four adjacent corners, and verifies each determinant under the parameter swap induced by its isomorphism.

The fixed label-three graph also records a necessary field distinction. If the endpoints of a free edge are connected by a fixed-edge path, its product itself is a signed cyclic product and lies in the integer ring of Q(sqrt5). Otherwise only its square is automatically forced into that field. Of the literal patterns, 72 for D5+2D4 and108 for 2D5+D4 have one product in the field and only the square of the other forced into it. All remaining patterns have both products in the field. None in this census has two independent square-only products. These statements do not classify admissible values.

The [source](../scripts/fresh_sqrt5_mixed_d4_d5_quad_patterns_20261017.py) and [saved output](fresh_sqrt5_mixed_d4_d5_quad_patterns_20261017.json) have passed generation and a fresh complete saved replay:

```sh
python scripts/fresh_sqrt5_mixed_d4_d5_quad_patterns_20261017.py
python scripts/fresh_sqrt5_mixed_d4_d5_quad_patterns_20261017.py --verify
```

This package inherits the previously frozen finite graph classifier, exact rational polynomial helpers and colored-isomorphism engine. It is not an independent implementation of those tools. Physical Lorentz signature, conjugate positivity in a valid field gauge, the allowed opposite products, and chamber coherence remain to be classified. Consequently this census alone gives neither a mixed incidence inequality nor a dimension bound.
