# Independent review of the low-polar obstruction

The [separate reviewer](../scripts/fresh_h9_h8_corner_anchors_20261014_polar_review.py) verifies the [uniform low-polar obstruction](fresh_h9_h8_corner_anchors_20261014_polar_gate.md) using exact literal root pairings, separate 128-bit outward-rounded intervals and sparse-polynomial identities. It imports neither the producer, its 192-bit evaluator nor SymPy. Its [certificate](fresh_h9_h8_corner_anchors_20261014_polar_review.json) reconstructs the strict lower bound `r²>2651/918>2`, contradicting the cap normalization.

Coverage forces a small contact with old wall 4. The old roots 0 and 2 are orthogonal to each other and to the Lorentz plane spanned by walls 4 and 11. Their required contacts contribute more positive norm than that plane can offset. The reviewer checks the literal orthogonality, all eight strict interval margins, the source-corner distances and the monotonicity of the Lorentz projection term.

The source contact-gap theorem and independently reviewed two-anchor theorem remain explicit dependencies. This rules out every low wall-11 center for the raw attachment using the twelve old-wall covering balls. Additional nonperpendicular walls, changed retained walls and different source chambers are outside this conclusion.

```sh
python3 scripts/fresh_h9_h8_corner_anchors_20261014_polar_review.py --verify
```
