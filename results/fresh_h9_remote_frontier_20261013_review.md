# Review of the 105/106-wall remote frontier construction

**PASS:** the exact cap formulas, compact quadrilateral, added H3+6A1 vertex and surviving strict escape are consistent. I source-reviewed both checkers, replayed their saved 80-digit certificates exactly, and reran the same rigorous interval implementation at 160 digits. All **8,255** changed predecessor/final enclosures nest inside their saved bounds; exact metadata agree. This is a mathematical/source review and precision replay, not a second implementation of all 106 algebraic normals.

The [review checker](../scripts/fresh_h9_remote_frontier_20261013_review.py) and [saved result](fresh_h9_remote_frontier_20261013_review.json) bind the frozen cap recipes, original checkers, interval engine and review source by hash. The construction and its scope are described in the [main report](fresh_h9_remote_frontier_20261013.md).

For each cap, the eight retained edge normals are orthogonal and have norm two. The anchor complement `b = a−(d/2)e` is perpendicular to this span, with `q_b=2−d²/2<0`. The prescribed profile projection has square 1 for the first cap and `(2+phi)/2` for the second. It is perpendicular to b. Subtracting `b/q_b` gives the anchor pairing exactly −1 and particular norm `q_base+1/q_b`.

Projecting the common center off the eight edge normals and then off b leaves a nonzero spacelike vector z. The final radical term therefore restores norm two without changing any prescribed pairing. The review independently checks these norm and anchor identities symbolically. Finite equalities are justified by the orthogonal decomposition, not by an interval merely containing the desired value. Strict inequalities choose the branch and certify all other relations.

The first cap has six right relations, three angles pi/3 and 95 dotted relations. The second has six right relations, two angles pi/3, one angle pi/5 and 96 dotted relations. Their complete 512-subset finite-neighbor censuses have no affine subset. Any subset containing the new wall and a dotted neighbor is already indefinite on that pair, so it cannot hide an omitted positive affine circuit.

In the common rank-seven face, the first new side projects to norm one; the other three side norms are two. Its pairing −1 with the norm-two anchor gives angle pi/4. The other consecutive pairings are zero. Both opposite pairs are disjoint. The checker exhibits all four finite vertices and verifies every other wall inequality strictly. These four points form an actual compact Lambert quadrilateral with labels `(2,2,4,2)`, rather than just a formal four-node Gram pattern.

The second cap has a new rank-nine spherical incident diagram H3+6A1. The review separately checks the positive-definite H3 Gram matrix. All four quadrilateral vertices survive the new inequality strictly, so the compact face remains present at 106 walls. The strict common center and nonpositive wall pairings also ensure that the retained walls remain actual facets by orthogonal projection of the center.

The final escape still has strictly positive norm, negative pairing with the common future center, and all 106 strict halfspace inequalities. It proves infinite volume. Neither the compact face nor the absence of affine diagrams establishes compactness of the full chamber.

```sh
python scripts/fresh_h9_remote_frontier_20261013_review.py --verify
```

The higher-precision replay leaves all source artifacts untouched. Where a predecessor checker compares a regenerated certificate with its saved version, the review requires exact metadata equality and rigorous interval nesting at the new scale. The ordinary 80-digit comparisons are replayed exactly first.
