# Review of the compact H4 prism in the 96-wall H9 chamber

The prism argument is sound. In the common five-facet subspace, the four side normals have the old `3–5–3` tetrahedral Gram matrix. The old base is perpendicular to all four sides; the new base has a negative coefficient in that transverse direction and is strict at each old tetrahedral vertex. Over each ray of the tetrahedral side cone, the two base inequalities bound a finite interval. The resulting cone has exactly four old and four new endpoint rays. All eight are future timelike and satisfy every other chamber inequality, proving an actual compact H4 prism face.

The [review script](../scripts/fresh_h9_face_bridge_20261013_review.py) separately derives the six-facet prism Gram matrix. It verifies its rank-five Lorentzian realization, the dotted base pair, and positive definiteness at all eight vertex supports. A rigorous replay at 160 decimal digits checks that every changed enclosure nests inside its saved 80-digit enclosure and all exact metadata agree. It confirms six new actual finite vertices in total, 19 finite edges and 52 open rays in the ten-vertex reached component, and the surviving strict spacelike escape.

The common five-facet Gram matrix is checked to be `2I`, hence type `5A1`. The construction also exposed a correction to descriptive diagram types in the preceding 94-wall report: its initial rank-five block is `A2+3A1`, its common tetrahedral block is `6A1`, and its vertex types are `H3+6A1` and `A2+7A1`. The archived exact matrices, vertex tests, and volume conclusion were already correct.

This is independent symbolic/source reasoning and a higher-precision replay of the existing interval implementation, not a separate all-root arithmetic implementation. The whole H9 chamber remains infinite-volume.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h9_face_bridge_20261013_review.py --verify
```

See the [review certificate](fresh_h9_face_bridge_20261013_review.json) and [construction report](fresh_h9_face_bridge_20261013.md).
