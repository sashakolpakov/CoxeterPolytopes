# Review and 90-digit replay of the 3,857-cap escape certificate

The [reviewed certificate](fresh_h9_split_multiple_edge_escape_20261013.json) passes a separate run at 90-digit integer interval precision. All 38,570 covector intervals lie inside the original 70-digit enclosures, all 3,857 branch IDs agree, and the common strict center and strict spacelike escape pass for all 72 old walls and all 3,857 reconstructed caps. The [audit result](fresh_h9_split_multiple_edge_escape_independent_20261013.json) records the source hash and counts.

Source review checked outward interval multiplication, reciprocal bounds for both denominator signs, integer square-root bounds, and Gaussian elimination. Replacing pivot entries by exact one and eliminated entries by exact zero is valid because these are algebraic identities of the underlying exact row operation. Certified nonzero pivots give rank nine, so the nine prescribed pairing equations have the archived one-dimensional solution space. The positive discriminant and nonzero leading coefficient define two real norm-two roots; norm two follows from the quadratic identity, not merely from a residual interval containing two.

One branch-identification gap was found and corrected before this replay: an enclosure not wholly inside the numerical proposal box does not automatically prove its exact value lies outside. The current code explicitly proves that the other branch has at least one coordinate interval disjoint from the rational proposal box. The selected branch's full enclosure lies inside it. Thus root identity is certified, even though decimal proposals were used to discover the bank.

The strict spacelike witness and strict timelike center prove infinite volume for the intersection of all these defined walls, hence for any subset retaining the 72 old walls. The other proposed angles need not be certified for this obstruction. No conclusion extends to omitted branches, labels, anchors, new edge generations or arbitrary completions of the seed.

This is a source-level proof review and a separate-precision replay using the existing interval implementation, not a second interval-arithmetic implementation. The [replay script](../scripts/fresh_h9_split_multiple_edge_escape_independent_20261013.py) needs no numerical proposal generation:

```sh
PYTHONPATH=src:scripts:. python scripts/fresh_h9_split_multiple_edge_escape_independent_20261013.py --verify
```
