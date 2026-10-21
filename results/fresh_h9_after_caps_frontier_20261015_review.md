# Independent review of the27-wall frontier

**PASS.** The [independent checker](../scripts/fresh_h9_after_caps_frontier_20261015_review.py) and [certificate](fresh_h9_after_caps_frontier_20261015_review.json) reconstruct the [27 new finite-vertex germs](fresh_h9_after_caps_endpoint_frontier_20261015.json) and the [60 vertex-free lines](fresh_h9_after_caps_vertex_free_edges_20261015.json). Arithmetic uses a separate SymPy algebraic field and192-bit rational outward sign bounds. No producer quartic arithmetic, chamber engine or profile generator is imported.

The full27-wall cone—84 finite vertices,660 exterior rays and no ideal rays—is inherited from the [independent cone review](fresh_h9_vertex_free_extension_20261015_review.md). This review does not repeat that cone clipping.

For the vertex-free lines, the checker independently enumerates every pair of exterior rays sharing at least eight active walls. The2929 candidates split into1518 without a finite-angle eight-clique,1351 with no positive-definite eight-root support, and60 positive cases. Each positive support belongs to exactly its two exterior generators and to no finite vertex. The checker independently reconstructs those generators from the physical wall nullspaces, verifies their strict exterior norms and negative mutual product, and checks the saved future timelike interior point and all wall signs. Their two-dimensional feasible cone contains both future null limits strictly before either exterior boundary; hence its hyperbolic intersection is a complete geodesic one-face with no finite endpoint. All60 supports are distinct. Completeness here depends explicitly on the inherited complete cone census.

For the finite-vertex frontier, all247 finite edges are reconstructed from the84 saved corners. Exactly three vertices use one of the two new facets, and all their27 outgoing directions remain open. Each direction is independently obtained from a column of the full corner inverse Gram. It is orthogonal to the starting vertex and has contact−1 with the released wall and zero with the other eight. Every first possible old-wall contact occurs strictly beyond the opposite null point. The saved strict spacelike witness and all27 physical products are checked for each germ. Consequently its null limit has exactly the eight positive-definite support walls incident; its Euclidean link is unbounded. These are ideal-boundary limits, not new bounded cusps.

The normalized support diagrams and inverse matrices are reconstructed independently. The checker re-enumerates all local positive-definite endpoint profiles:759 finite profiles and36 symbolic isolated-dihedral tails. It checks every contact set, exact Schur cost, positive residual norm, component type, trial count and rejection count. A finite spherical diagram cannot contain a cycle, a four-neighbor new node, a three-neighbor node with a label above3, or a label at least6 in a component larger than two. These standard finite-diagram restrictions justify the enumeration bounds; the isolated `I2(m)` families remain symbolic for every integer `m≥6`. NoI2(10) support component occurs in these27 cases.

Generation and saved replay pass:

```sh
python scripts/fresh_h9_after_caps_frontier_20261015_review.py --verify
```

The scope is the complete60-line census and the27 new-vertex germs of the27-wall chamber. The older235 open germs, the later28-wall chamber, numerical cap scouts, global compatibility of proposed endpoint walls, and affine endpoint extensions are not re-proved here. No finite-volume conclusion follows from these local frontier certificates.
