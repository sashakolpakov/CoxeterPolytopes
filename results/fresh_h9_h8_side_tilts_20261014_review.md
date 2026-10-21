# Independent review of the side-lift rigidity lemma

The [reviewer](../scripts/fresh_h9_h8_side_tilts_20261014_review.py) reconstructs all 121 source Gram entries using separate rational quadratic-field arithmetic, verifies all five strict angle gaps, and checks all 2,048 assignments of right versus nonright side contacts. The [certificate](fresh_h9_h8_side_tilts_20261014_review.json) leaves only the all-right assignment, agreeing with the [side-lift theorem](fresh_h9_h8_side_tilts_20261014.md).

The normalized lift formula follows directly by resolving each side normal into its projection onto the source hyperplane and its normal component. The pi/3 propagation and the forbidden positive pairing for a source-right pair cover the large component. Exact gaps handle labels 3 and 4 at the two remaining nodes; monotonicity of sine handles every label at least 5.

This rules out tilting the eleven actual sidewalls of this fixed source. Additional walls disjoint from its facet are outside the theorem.

```sh
python3 scripts/fresh_h9_h8_side_tilts_20261014_review.py --verify
```
