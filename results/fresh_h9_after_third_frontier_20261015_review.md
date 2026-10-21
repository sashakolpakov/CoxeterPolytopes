# Independent review of the 28-wall frontier

**PASS.** The [independent checker](../scripts/fresh_h9_after_third_frontier_20261015_review.py) and [certificate](fresh_h9_after_third_frontier_20261015_review.json) reconstruct all [24 new finite-vertex germs](fresh_h9_after_third_endpoint_frontier_20261015.json) and the complete [71 vertex-free lines](fresh_h9_after_third_vertex_free_edges_20261015.json). Arithmetic uses a separate SymPy algebraic field and 192-bit rational outward sign bounds. No producer quartic arithmetic, chamber engine or profile generator is imported.

The complete 28-wall cone—87 finite vertices, 727 exterior rays and no ideal rays—is inherited from the [independent cone review](fresh_h9_replacement_escape_extension_20261015_review.md). This review does not repeat that cone clipping.

The line checker enumerates every pair of exterior rays sharing at least eight active walls. The **3,261 candidates** split into 1,666 without a finite-angle eight-clique, 1,524 with no positive-definite eight-root support, and 71 positive cases. Each positive support belongs to exactly its two exterior generators and to no finite vertex. Independent physical nullspaces recover these generators. Their positive norms, negative mutual product, future timelike interior point and all wall signs are checked. Both future null limits lie strictly inside the feasible two-dimensional cone face, so its hyperbolic intersection is a complete geodesic one-face with no finite endpoint. All 71 supports are distinct. Completeness uses the inherited full cone census explicitly.

For the finite frontier, the checker reconstructs all 249 finite edges. Of these, 247 join old vertices, one joins an old and a new vertex, and one joins two new vertices. The three vertices on the added facet consequently give **24 open germs**; 261 older open germs survive and are not re-enumerated as new profiles.

Each new direction is independently obtained from the full corner inverse Gram. It is orthogonal to the starting vertex, has contact minus one with the released wall and zero with the other eight, and points into the chamber. Every first possible old-wall contact occurs strictly beyond the opposite null point. The saved strict spacelike witness and all 28 physical wall products are checked for each germ. Thus its null limit has exactly the eight positive-definite support walls incident and an unbounded Euclidean link. These limits are not bounded cusps.

The checker reconstructs every normalized support diagram and inverse matrix, then independently re-enumerates **247 finite spherical endpoint profiles and 18 symbolic isolated-dihedral tails**. It compares every contact set, exact Schur cost, residual norm, component type, trial count and rejection count. A finite spherical diagram has no cycle, no four-neighbor new node, no three-neighbor new node with a bond above label 3, and no label at least 6 in a component larger than two. These restrictions justify the finite enumeration; isolated `I2(m)` tails remain symbolic for every integer `m≥6`. No inert `I2(10)` support component occurs in these 24 germs.

Generation and saved replay pass:

```sh
python scripts/fresh_h9_after_third_frontier_20261015_review.py --verify
```

This review covers the complete 71-line census and the 24 new germs of the 28-wall chamber. It does not certify a numerical cap scout, global compatibility of an endpoint cap, affine endpoint extensions, or any finite-volume completion.
