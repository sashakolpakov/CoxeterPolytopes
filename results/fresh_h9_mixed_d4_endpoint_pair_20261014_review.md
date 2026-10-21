# Independent review of the mixed D4 endpoint pair

**PASS.** The two proposed normals give a coherent fourteen-wall Coxeter chamber, close both targeted old open edges, and create three finite vertices, including a vertex incident to both new walls. Both cut the previous escape. A new strict rational-coordinate spacelike escape remains, so the chamber still has infinite volume.

The [independent checker](../scripts/fresh_h9_mixed_d4_endpoint_pair_20261014_review.py) and [certificate](fresh_h9_mixed_d4_endpoint_pair_20261014_review.json) bind the frozen [primary source](../scripts/fresh_h9_mixed_d4_endpoint_pair_20261014.py), [primary data](fresh_h9_mixed_d4_endpoint_pair_20261014.json), original chamber and frontier. They reuse the preceding independent `Fraction` implementation of the quartic field, with no primary producer imports. The complete cone is independently reconstructed by Gaussian nullspaces over all **2,002 nine-wall supports**, rather than the producer's incremental clipping algorithm.

The literal cap coordinates lie in `Q(sqrt5)` inside the retained quartic ambient field. Both normals have square 2 and their mutual product is exactly zero. All 25 new pair relations pass independently: 13 right angles, four π/3 angles, two π/5 angles and six dotted pairs. All old pair relations are also rechecked. Their exact pairings with the old strict escape are positive.

All twelve old finite vertices pair strictly negatively with both caps. The specified new center is checked physically, is future timelike, and is strict for all fourteen walls. Projection of this center onto each wall gives an actual facet point: future timelike, incident to exactly that wall, and strict for every other inequality. Thus the wall count is an actual facet count.

The independent complete cone has exactly **15 finite and 91 spacelike rays, with no ideal ray**. The old roots span the ambient ten-dimensional vector space, and their positive sum after extension is timelike and strictly negative on every nonzero feasible cone direction. This gives the bounded affine section needed for the complete ray interpretation. The saved normalized ray vectors, active sets and norms agree with the independent support enumeration.

The three new finite active sets and their spherical types are:

| Active walls | Spherical type |
|---|---|
| 0,2,3,4,5,7,9,10,12 | A8 + A1 |
| 0,1,3,4,5,7,9,10,13 | D7 + I2(5) |
| 0,3,4,5,7,9,10,12,13 | E8 + A1 |

Each active Gram has rank nine and is positive definite. The first two points are the unique new finite endpoints on the saved old open edges 6 and 78. The third is their joint-cap vertex.

The complete finite-corner edge queue is also reconstructed. Every rank-eight support occurs at exactly two cone rays. Their positive sum has precisely those eight equalities and every other inequality strict, proving actual adjacency. For an open direction the checker additionally finds an exact point beyond the light-cone crossing, strict at every nonsupport wall.

This gives 64 directed finite closures, hence **32 distinct finite edges**, and 71 open directions at the 15 finite vertices. The finite-edge graph is connected. In the final ray numbering, the new path is

`14 — 18 — 74 — 67 — 66`.

Its three internal vertices are new and its two endpoints are old. Four finite edges are added, increasing the connected finite-edge graph's cycle rank from 17 to 18. This is a finite-neighbor connection, rather than merely two isolated cap facets.

Finally, the saved rational-coordinate escape is recomputed directly: its square is strictly positive, its pairing with the future center is negative, and all fourteen wall pairings are strictly negative. That certifies the remaining open ideal-boundary region and infinite volume. No arbitrary-cap exhaustion, vertex-free-edge census, compactness or new finite-volume dimension is claimed.

Generation and saved replay pass:

```sh
python scripts/fresh_h9_mixed_d4_endpoint_pair_20261014_review.py --verify
```

No formula, geometric or scope issue was found. The review binds source/data files only, leaving the primary report free to add this link.
