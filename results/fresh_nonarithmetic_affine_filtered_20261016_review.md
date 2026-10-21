# Independent reconstruction of the filtered 52-wall H9 chamber

The filtered three-cap extension passes independent local reconstruction,
complete-cone generation, and fresh saved replay. It has **52 actual facets,
178 finite vertices, no ideal extreme rays, and 2,279 exterior extreme rays**.
A strict spacelike escape proves infinite volume. This does not establish
that a finite-volume continuation exists.

The [local reviewer](../scripts/fresh_nonarithmetic_affine_filtered_20261016_cluster_review.py)
and [certificate](fresh_nonarithmetic_affine_filtered_20261016_cluster_review.json)
start from the independently reviewed 49-wall geometry. They invert the
eight-root support Gram in reverse order, reconstruct the projected anchor
and its positive orthogonal kernel in the exact nested field, and define
each new cap by its positive square-root formula. These identities prove
the norm and nine prescribed contacts; intervals containing zero are only
sanity checks of those identities. Separate rational interval arithmetic
certifies 120 unspecified old contacts and all three mutual pairs as dotted,
all three cuts of the preceding escape, and all 52 actual facets.

The [complete-cone reviewer](../scripts/fresh_nonarithmetic_affine_filtered_20261016_cone_review.py)
and [certificate](fresh_nonarithmetic_affine_filtered_20261016_cone_review.json)
start with the independently certified complete 49-wall cone. Three exhaustive
clips create 25, 56 and 55 crossing rays. Every crossing has an independently
certified rank-eight common-covector minor and strictly positive convex
coefficients. Every new ray is exterior. Exact inherited and crossing
identities give the active zeros; all nonactive inequalities are strict.
The final pass checks all 127,764 ray/inequality pairs and independently
reconstructs the renewed strict escape. The 178 finite vertices are unchanged,
and none belongs to a new facet.

The affine-circuit exclusion is inherited from the separate
[primary filtered-cluster report](fresh_nonarithmetic_affine_filtered_20261016.md)
and [independent old-49 inventory](fresh_h9_persistent_cusp_screen_20261016.md).
The reviewers here certify local geometry and the complete cone; they do
not rerun the new affine-domain enumeration or the older right-to-33 gate.
The stronger subsequent [all-real cutter review](fresh_nonarithmetic_escape_cutter_domain_20261016_review.md)
uses this complete cone.

```sh
python3 scripts/fresh_nonarithmetic_affine_filtered_20261016_cluster_review.py --verify
python3 scripts/fresh_nonarithmetic_affine_filtered_20261016_cone_review.py --verify
```
