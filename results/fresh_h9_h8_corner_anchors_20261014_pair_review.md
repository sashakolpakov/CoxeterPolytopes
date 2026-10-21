# Independent review of the two-anchor obstruction

The [separate reviewer](../scripts/fresh_h9_h8_corner_anchors_20261014_pair_review.py) verifies the [two-anchor obstruction](fresh_h9_h8_corner_anchors_20261014.md): any coherent placement of the fixed compact H8 source with the twelve retained walls has **at most one projected norm below `sqrt(5)-1`**. The [certificate](fresh_h9_h8_corner_anchors_20261014_pair_review.json) passes generation and saved replay.

Exact matrix inversion over the rational quartic field reconstructs the three exceptional source rays from their contacts, checks all eleven contacts and derives their Gram matrix. A separate 128-bit outward-rounded interval implementation then proves all 64 ordinary-pair comparisons and all eight special-pair margins. Neither the pair-gate producer nor its 192-bit interval code is imported. The same-ray inequality and monotonicity conditions are checked separately.

The frozen source contact-gap theorem remains a dependency. This is a necessary condition for an attachment, not a construction or a volume certificate.

```sh
python3 scripts/fresh_h9_h8_corner_anchors_20261014_pair_review.py --verify
```
