# Independent review of the Q(√5) triangular-face weights

The proposed rule **passes independent review**: every triangular two-face of a compact arithmetic Coxeter polytope over Q(√5) receives weight at least two, including arbitrary extra white nodes; a full vertex diagram of rank n≥5 receives at most 2n−2. This establishes no quadrilateral lower bound or dimension theorem.

The [review checker](../scripts/fresh_sqrt5_triangle_weights_20261015_review.py) and [certificate](fresh_sqrt5_triangle_weights_20261015_review.json) bind the final [primary source/data](fresh_sqrt5_triangle_weights_20261015.json) and independently reviewed complete Lannér catalogue. The checker uses the standalone Fraction Q(√5) arithmetic in the separate [contact-review helper](../scripts/fresh_sqrt5_quad_contacts_20261015_review.py), whose source is also bound. It imports no primary arithmetic or primary graph code.

## Literal reconstruction

All 76 markings are regenerated directly from the complete catalogue: 12 order-three diagrams, six order-four diagrams with four three-black markings each, and four order-five diagrams with ten markings each. Every corner is checked by a newly reconstructed full Gram and exact congruence inertia, rather than reusing the producer's leaf pivots. All marked distances, bonuses and totals match.

Exactly seven markings have base distance-two total one. Each is repaired by the two H4 corrections alone: one additional unit on its 5-bond and one on its endpoint pair. Thus the smaller H4-only rule also pays every triangle. The full rule's optional bonuses on the 5-bond of I2(5) and H3 remain harmless and may help another face argument.

## Stability under arbitrary white extensions

For a compact triangular face, each white support plus two black nodes is a finite corner. A minimal indefinite subdiagram of the full star must therefore contain all three black nodes. Compactness excludes affine principal subdiagrams, so it belongs to the bound Lannér catalogue. This justifies using the complete 76 markings rather than assuming the full star has at most five nodes.

In a finite spherical diagram the graph is a forest. Any already connected marked pair retains its unique path and hence its distance under spherical extension. Disconnected pairs can only gain a nonnegative base weight. A finite connected diagram containing a 5-bond is I2(5), H3 or H4; its marked 5-bond therefore retains its bonus as the component grows. An H4 cannot grow further, so its endpoint correction is also stable. Arbitrarily many orthogonal components do not affect these conclusions.

The checker independently reconstructs all induced embeddings among I2(5), H3 and H4 and all their relabellings. These finite checks verify the literal bonuses; the finite-spherical classification proves there is no larger untested attachment family. The proof does not confuse an H4-shaped subgraph in a different labelled component with an actual H4.

## Full vertex budget

The independent count uses the tree identity

\[
\#\{\text{pairs at distance 1 or 2}\}
=(r-1)+\sum_v\binom{\deg v}{2}.
\]

This gives 2r−3 on paths and 2r−2 on D/E trees, with zero for A1. The correction gives totals 2,4,7 on I2(5),H3,H4 respectively; F4 remains total five. The only irreducible component with deficit 2r−weight below two is H4, with deficit one and rank four. For total rank at least five, a single component has deficit at least two, while two or more components have at least two combined. Hence the all-rank upper bound is 2n−2. The rank-four H4 exception is retained explicitly.

Every saved component total and the rank-partition maxima through 100 are independently reconstructed. These finite counts check the formulas; the degree identity, finite-spherical classification and deficit argument prove all ranks without a cutoff.

Independent generation and saved replay, recursive source manifests, syntax and local links pass. No issue was found in the triangle, arbitrary-extension or vertex-budget claims. All quadrilateral and global dimension questions remain outside this package's scope.
