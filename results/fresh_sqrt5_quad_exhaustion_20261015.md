# Complete minimal-pair unions over Q(sqrt5)

The [enumerator](../scripts/fresh_sqrt5_quad_exhaustion_20261015.py) and [certificate](fresh_sqrt5_quad_exhaustion_20261015.json) exhaust the union of two opposite minimal Lannér pieces in a compact arithmetic Q(sqrt5) quadrilateral star having **no marked corner pair at graph distance at most two**. There are **5,440 distinct labelled union supports and 1,180 marked isomorphism classes**. This is a complete necessary local list; it does not assert that any listed support extends to a full face star or a compact polytope.

The finite alphabet and minimal pieces come from the [complete catalogue](fresh_sqrt5_lanner_catalogue_20261015.md) and its [independent review](fresh_sqrt5_lanner_catalogue_20261015_review.md). They include all labels 2,3,4,5,6,10, the symbolic dotted order-two family, all 12 triangle types, six order-four types and four order-five types. Dotted contacts are represented by label zero and remain undetermined.

Let the cyclic black labels be 0,1,2,3 and let W be the white support of a quadrilateral. Its four corners have spherical diagrams W plus one black node from each opposite pair. The diagrams W+{0,2} and W+{1,3} are indefinite; choose minimal indefinite pieces in them. Each chosen piece contains both of its opposite black nodes, because deleting either leaves a subdiagram of a spherical corner. Compactness excludes affine alternatives, so the chosen pieces belong to the Lannér catalogue.

For each piece, black names are fixed while white permutations are quotiented. All placements of black nodes in every catalogue graph are retained. This gives 153 representatives: one dotted pair, 49 marked triangles, 40 order-four pieces and 63 order-five pieces. The enumerator tries every injective shared-white identification, including disjoint white sets. Swapping the two opposite pairs justifies using unordered piece indices; no polytope automorphism is assumed.

Every pair not already fixed by a piece is assigned every finite label in the complete alphabet. These free pairs belong to at least one actual corner, so they cannot be dotted or parallel. A direct adjacent-black contact, or a common white neighbor of adjacent black nodes, violates the no-distance-two premise. All four induced corner Grams must be positive definite. The primary checker uses exact forest Schur elimination over Q(phi), including isolated I2(6) and I2(10) factors. Cycles have a nonpositive quadratic witness, while a bond at least six cannot meet another nonright bond in a spherical component. Assigning remaining contacts right is the weakest acute Gram, making early failures safe monotone prunes. If the two chosen indefinite pieces are disjoint and orthogonal, their union has at least two negative directions and is excluded.

The exhaustive traversal has no cutoff. It checks 127,875 shared-white identifications, of which 12,452 pass the initial corner and short-path tests, and visits 207,319 DFS nodes. It rejects 11,781 orthogonal-piece completions and retains 5,502 completions before literal duplicate removal. The 5,440 distinct supports have six to ten nodes. Every isomorphism class includes explicit full label-preserving maps that preserve black/white colors and the two opposite pairs.

The [separate coverage replay](fresh_sqrt5_quad_exhaustion_20261015_incidence_review.md) uses finite Dynkin recognition instead of the primary Schur arithmetic, a different shared-white generator, reversed contact order and an independent canonicalization. Its complete 196,742-node traversal recovers exactly the same 5,440 supports and 1,180 classes. Generation and saved replay pass for both implementations.

Outside white vertices are not enumerated here. Cyclic field membership of the union, conjugate positivity, Lorentz inertia, full-star rank, arbitrary external white components, and the extra triangle-weight corrections require separate arguments. In particular, this list alone gives no new Q(sqrt5) dimension bound. The separately produced arithmetic filters may remove supports, and some surviving local corners already earn an entire-H4 correction; neither fact is silently used to restrict this enumeration.

```sh
python3 scripts/fresh_sqrt5_quad_exhaustion_20261015.py --verify
python3 scripts/fresh_sqrt5_quad_exhaustion_20261015_incidence_review.py --verify
```

The optional running snapshot stores only completed outer-piece ranges. `--resume` checks its source hashes and continues at the next outer piece; an interrupted current piece is recomputed without losing later siblings. The final certificate above is the authoritative completed result.
