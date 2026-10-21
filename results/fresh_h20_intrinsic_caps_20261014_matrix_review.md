# Independent exact matrix replay of the 587-facet construction

The [independent checker](../scripts/fresh_h20_intrinsic_caps_20261014_matrix_review.py) confirms the [587-facet construction](fresh_h20_intrinsic_caps_20261014_quarter_face_exact.md), including its infinite-volume obstruction. Its [certificate](fresh_h20_intrinsic_caps_20261014_matrix_review.json) binds all three construction archives and this checker by SHA-256.

The checker imports none of the construction arithmetic. It splits every normal into two rational coefficient vectors, `a+b*phi`, and computes the entire Gram matrix using FLINT rational matrix multiplication and `phi²=phi+1`. Exact rational comparisons decide signs in the chosen real embedding.

It independently verifies:

- Signature `(20,1)` from positive rational spatial LDL pivots and the null-plane block.
- All 587 root norms and all **171,991** pair relations, including right angles, labels 3 and 5, parallel and divergent pairs.
- A common strict timelike point.
- All 252 saved finite endpoints, every incident and nonincident inequality, and positive definiteness of each rank-20 incident Gram.
- 565 saved timelike witness points collectively establishing that every one of the 587 inequalities defines an actual facet.
- The final escape vector of square `1/2`, its time orientation, and strict negativity at all 587 walls. Hence the chamber has infinite volume.

The independent replay checks the saved coordinates, not the recipe-selection search or a complete corner census. The quarter-face producer separately reconstructs the 60 selected recipes and verifies their probe coverage. No finite-volume example or new dimension is established.

Generation and saved replay pass:

```sh
python3 scripts/fresh_h20_intrinsic_caps_20261014_matrix_review.py --verify
```
