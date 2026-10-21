# Octahedral decomposition: proof and reproduction map

The [nine-page manuscript](../manuscripts/ideal_ra4_octahedral_decomposition.tex)
and [PDF](../manuscripts/build/ideal_ra4_octahedral_decomposition.pdf) prove
a restricted structure theorem. Every facet must be isometric to the
regular ideal octahedron or its triangular-face reflection double. Under
that hypothesis the four-polytope splits into a chain of regular ideal
24-cells, with opposite cuts in each intermediate cell. This is not a
classification of arbitrary ideal right-angled four-polytopes.

The proof was developed directly from intrinsic facet geometry. A
four-dimensional classification, minimum-volume theorem or universal
volume-divisibility assertion is not used to establish the decomposition.
The known volume of the regular24-cell supplies its volume normalization.

## The arguments that must be proved internally

|Manuscript part|Content and scope|
|---|---|
|Cusp rectangles|The doubled facet has three aspect-two cusps; its internal triangle crosses each rectangle at its middle. The permitted ambient boxes are cubic or2:1:1.|
|Closed surfaces|Actual quadrilateral diagonals glue the internal triangles. The four-segment cusp links prove the compactification is a closed two-manifold.|
|Eight triangles|A primal/dual spanning-tree argument proves Euler characteristic<=2. Seven distinct same-component facets arise from three neighbors and three asymptotic opposites; divisibility then forces eight.|
|Six-cusp graph|The graph is simple, four-regular on six vertices, hence octahedral. The eight distinct triangles give the dual cube.|
|Actual hyperplane|Equal unit-normal cusp sums give affine cube normals in a Lorentz4-space. Its six actual cusp rays and convexity identify the entire octahedral section.|
|Valid cuts|A finite section point lies in a chosen facet or ridge interior. An outside facet cannot touch it, so the cut introduces no other sliced facets or finite vertices. Distinct component sections have disjoint closures.|
|Recognizing24-cells|A whole regular octahedral facet and right angles determine the adjacent whole octahedral facet. Propagation through the model24-cell proves equality by both halfspace and convex-hull containments.|
|Opposite cuts|A nonopposite pair of closure-disjoint24-cell facets has a common neighbor. That original facet would contain three disjoint octahedra, contrary to the two allowed intrinsic facet volumes.|
|Uniqueness and counts|The cell tree is a path. Facet symmetries extend to the24-cell, ruling out gluing twists. Face counts follow with no extra cusp identifications.|

The manuscript contains five explanatory figures: the doubled cusp
rectangle, its ambient box belt, the octahedral/cube incidence graphs,
the internal quadrilateral diagonal, and the final path. Captions separate
metric link diagrams from combinatorial schematics.

## Exact checks and their limits

The finite checker is
[check_ideal_ra4_octahedral_path_roots_20261018.py](../scripts/check_ideal_ra4_octahedral_path_roots_20261018.py).
It checks all24 D4 roots,24 model ideal vertices and276 root pairs, including
all96 nonopposite disjoint pairs and1152 reflected octahedral vertices.
It uses exact rational arithmetic. Its successful root execution is
recorded in the [finite-check ledger](../results/ideal_ra4_octahedral_path_roots_check_20261018.json).
This verifies explicit model facts, not the universal surface, cutting or
recognition arguments.

Existing [double-cusp and cube certificates](ideal_ra4_structure_reproduce_and_map.md)
support the same local geometry. They are inherited finite calculations;
the present proof includes its own universal geometric deductions.

There are separate written reviews of the
[decomposition](../results/review_ideal_ra4_octahedral_block_decomposition_mathematical_20261018.md),
[opposition argument](../results/review_ideal_ra4_octahedral_block_path_rigidity_20261018.md)
and [full illustrated translation](../results/review_ideal_ra4_octahedral_decomposition_tex_mathematical_20261018.md).
These are internal mathematical reviews. No external refereeing, novelty
claim or complete geometric Lean verification is made. The existing Lean
integer-volume-gap module does not formalize this structure theorem.

## Reproduction

From the repository root, with Python3 and a TeX installation:

```sh
python3 -B scripts/check_ideal_ra4_octahedral_path_roots_20261018.py
python3 -B scripts/build_illustrated_ra_manuscript_20261018.py manuscripts/ideal_ra4_octahedral_decomposition.tex --check-only
```

Run these serially. Both acquire the shared proof-batch lock and require
40GiB free disk and6GiB available RAM before starting. The build has a
59-second,1GiB process-tree cap, runs one TeX process at a time, and removes
its temporary directory on success or failure. Check-only mode leaves
the preserved PDF and ledger untouched. A fresh PDF need not have the
same binary hash because TeX may embed timestamps.

The [final build ledger](../results/ideal_ra4_octahedral_decomposition_build_20261018.json)
records nine pages, five source figures, two passes, no warnings,85.9MB
peak process-tree RSS and temporary-file removal. The
[freeze manifest](../results/ideal_ra4_octahedral_decomposition_20261018_freeze.json)
binds the manuscript, PDF, direct proof, reviews, checker and build sources.
