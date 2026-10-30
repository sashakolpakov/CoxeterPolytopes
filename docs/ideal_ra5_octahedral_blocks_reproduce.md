# Restricted ideal RA5 octahedral blocks: proof and reproduction

For current reading, see the [current sharp right-angled manuscript](../arxiv/sharp_dimension_bounds_8_and_4/main.tex) and the
[submission guide](arxiv_submission_editions_20261028.md). The supporting
arguments and recorded verification scopes below remain distinct from those
current submission builds. Standalone draft typesetting routes have been retired.

The [new theorem](../results/ideal_ra5_octahedral_blocks_20261018.md)
strengthens the restricted facet counts to36 T and84 D, from30 and66.
It keeps the explicit C/D/T spectrum hypothesis of the
[predecessor](ideal_ra5_cusp_compatibility_reproduce.md). No dimension
ceiling or exclusion of all ideal RA5 polytopes follows.

The supporting argument
contains the universal written argument. The
[certificate](../results/fresh_ideal_ra5_octahedral_blocks_20261018.json)
records the four complete24-vertex ridge stars, their exact Gram matrix,
common null ray and exhaustive six-vertex coloring arithmetic.

```sh
python3 -B scripts/fresh_ideal_ra5_octahedral_blocks_20261018.py --verify
python3 -B scripts/check_ideal_ra5_octahedral_blocks_saved_20261018.py
```

Generation and separate-process replay passed in chunks fb0c22 and44d545,
each taking about0.07 seconds. It uses only the Python standard library
and under one MiB of inherited input data. All predecessors remain unchanged.

|Proof component|Verification scope|
|---|---|
|Actual C,D,T and complete T wall/vertex data|Inherited exact reflection-doubling proof and certificates.|
|Four24-vertex ridges, Gram matrix, null sums and104 polygon faces|Fresh exact reconstruction and saved replay.|
|Hexagon propagation and locally four-cycle graph|Internal written geometry, not inferred solely from local labels.|
|Octahedral components and common-neighbor bound two|Internal graph and Lorentz-space proofs.|
|Five-clique exclusion in an ideal RA5 polytope|Internal strict-interior projection proof.|
|Four-colorability at30 vertices and degree22|Published Andrasfai--Erdos--Sos1974 Theorem1.1, explicit external input.|
|Color-class obstruction|Written proof and all4096 colorings per part-size case.|
|Divisibility6 and y=3z-24|Inherited universal cusp/carrier counting proof for this restricted spectrum.|
|At least36 T and84 D|Composition of the preceding lemmas. No complete Lean verification.|

The [independent graph review](../results/fresh_ra5_t30_abstract_graph_review_20261018.md)
checks the abstract obstruction conditional on the geometric block hypotheses.
The [independent geometric review](../results/ideal_ra5_octahedral_blocks_geometric_review_20261018.md)
checks those hypotheses, the whole composition and the primary coloring source.
The standalone checker reconstructs all788 relevant wall-pair intersections
and104 actual polygons, including rank checks and source bindings; its final
replay passed in87bd0c. The three-page PDF
built in eeadd1 with no final-pass warnings; its unused build files were removed.
The resource-safe Lean runner is available, but Lean follows mathematical
progress and no compilation was needed for this result.
