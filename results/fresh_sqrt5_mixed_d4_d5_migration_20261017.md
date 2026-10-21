# Mixed D4/D5 endpoint migration and triangles

The complete finite endpoint and triangle tests for the three small-component budget obstructions leave a short list of possible two-special triangle corners:

| Full special vertex forest | Dimension | Same-type endpoint attachments, summed over ports | Possible two-special triangle corners per vertex | Colored triangle classes |
|---|---:|---:|---:|---:|
| 3D4, calibration | 12 | 12 | 9 | 1 |
| D5+2D4 | 13 | 24 | 11 | 4 |
| 2D5+D4 | 14 | 30 | 13 | 4 |
| 3D5+A2 | 17 | 38 | 16 | 4 |

Here an endpoint attachment is a literal placement of one new root on the retained roots, not a claim that a framing product or a global compact polytope exists. Triangle acceptance uses all spherical corners and exact necessary physical/conjugate arithmetic. These possibilities form a safe upper list; this package does not reconstruct each accepted local face cone independently.

For each D4 component the three allowed triangle corners are its center–leaf pairs. Number a D5 with center0, short leaves1,2, and long arm0–3–4. Its five allowed pairs are `{0,1}`, `{0,2}`, `{1,3}`, `{2,3}`, `{3,4}`. The A2 edge contributes its single pair. No accepted two-special triangle uses roots from different original components. Thus the counts are exactly the sum of these component contributions, even though endpoint migration itself can cross components.

## Complete finite endpoint enumeration

Every destination component is D4, D5 or A2. Hence the added root has only label-three contacts and degree at most three. The checker enumerates every contact subset of the entire retained root set of sizes one, two and three, then requires precisely the specified finite forest. A disconnected root is impossible because none of the target forests contains A1. This is a complete finite classification for the four displayed forests, with no optimization assumptions.

D4-root replacements are unique copies. D5-root replacements need not be. In particular, deleting the long-arm end of D5 leaves a D4; the new root can extend a different D4 component, moving the D5 component. The 3D5+A2 case also permits exchanging an A2 component with the A2 left after deleting a D5 center. The saved output retains every literal attachment and records the old components recruited by the new component.

## Finite triangles and the field distinction

For every attachment and every other released root, the checker tests all possible finite framing labels2,3,4,5,6,10. This is an exact field restriction, not a numerical label cutoff: the two-cycle places `4cos²(pi/m)` in Q(sqrt5), and the real quadratic cyclotomic possibilities8 and12 belong to Q(sqrt2) or Q(sqrt3) instead.

When the new root migrates between components, its framing bond can be a bridge. Only its square is then forced into the ground field. Assuming the framing product itself belongs to Q(sqrt5) would wrongly omit labels4,6 and10. The checker keeps these labels and applies the frozen exact cycle-field gauge and conjugate tests to each literal graph.

There are no all-special triangles: every root pair of a triangle star occurs in one of its finite corners, so three simply-laced special corners would make the entire norm-two Gram matrix rational. Its required physical negative direction would survive conjugation, contradicting arithmetic conjugate positivity.

The four non-calibration triangle families replace one D5 by `A1+A2+I2(5)`, `A1+H4`, or `A3+I2(5)`, or replace a D4 by `2A1+I2(5)`. In the 3D5+A2 case, the fourth family replaces A2 by I2(5). Orthogonal components persist as shown literally in the output.

The [source](../scripts/fresh_sqrt5_mixed_d4_d5_migration_20261017.py) and [saved output](fresh_sqrt5_mixed_d4_d5_migration_20261017.json) preserve all ports, attachments, accepted triangles and explicit colored isomorphism maps. The arithmetic engine is inherited from the frozen complete triangle catalogue; this is a new exhaustive application, not an independent implementation of that engine. Generation and a fresh saved replay pass:

```sh
python scripts/fresh_sqrt5_mixed_d4_d5_migration_20261017.py
python scripts/fresh_sqrt5_mixed_d4_d5_migration_20261017.py --verify
```

All-parameter quadrilateral closure remains necessary before deriving a mixed-type incidence inequality. No dimension bound, universal corner discount, or compact ambient realization follows from this package.
