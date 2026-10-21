# Independent check of the radius-four replacement

The [radius-four replacement](fresh_h20_intrinsic_caps_20261014_radius4_replacement.md) passes a [separate FLINT matrix check](../scripts/fresh_h20_intrinsic_caps_20261014_replacement_review.py), with [saved certificate](fresh_h20_intrinsic_caps_20261014_replacement_review.json). Unchanged geometry is inherited from the independent full [587-wall replay](fresh_h20_intrinsic_caps_20261014_matrix_review.md).

The checker verifies the new norm 16 and all 586 changed pair relations, the surviving strict center and old facet witnesses, and all 587 inequalities at the new endpoint. Its norm is `-8/3`, and its rank-20 incident Gram is positive definite. The other 251 certified finite vertices remain feasible. Thus the one-wall replacement preserves 587 actual facets and a set of 252 certified finite vertices.

The new wall cuts the annular escape strictly. The earlier eighth-vertex escape remains strict, so the replacement still has infinite volume. No complete corner census or bounded-link claim for the new parallel pair is made.

Generation and saved replay pass:

```sh
python3 scripts/fresh_h20_intrinsic_caps_20261014_replacement_review.py --verify
```
