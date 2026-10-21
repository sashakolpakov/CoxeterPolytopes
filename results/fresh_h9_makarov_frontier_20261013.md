# Complete frontier of the 12-wall Makarov-face seed

The exact 12-wall seed has **10 finite vertices, 25 finite edges, and 40 open directions** from those vertices. The added cap creates no finite vertex. Every open direction releases one of the four H4 walls. Its only later collision with an existing wall is with the added cap, beyond the light cone, so it has no other finite or ideal endpoint in the current chamber.

The complete projective census gives **50 extreme rays: 10 timelike and 40 spacelike, with no null extreme ray**. All 40 spacelike rays lie on the new cap and are precisely the projective endpoints of those open directions. The cap removes the four pure transverse recession rays of the product seed but replaces them with 40 exterior vertices. This explains both its useful gap cut and its failure to close any edge of the old finite star.

The [exact checker](../scripts/fresh_h9_makarov_frontier_20261013.py) and [certificate](fresh_h9_makarov_frontier_20261013.json) use the separately implemented squarefree-radical arithmetic from the [seed audit](fresh_h9_makarov_face_extension_20261013_independent.md). They import no arithmetic from the original seed producer.

For the finite-vertex census, all 220 nine-wall supports are tested: 192 contain a dotted pair, 18 of the remaining supports are not positive definite, and the other 10 give feasible timelike vertices. Independently of that projection calculation, exact Gaussian nullspaces exhaust the same 220 supports for the full projective cone: 10 are rank deficient, 160 are infeasible, and 50 give distinct extreme rays. A future timelike vector expressed with strictly positive coefficients on all 12 normals proves that the chosen affine slice is bounded, so the latter census is complete.

Every finite vertex supplies nine directed edges. Exactly 50 directed edges pair up into the 25 finite edges; the other 40 have an exact feasible spacelike continuation on their eight-wall support. Along each edge the checker compares all existing collision parameters exactly. The new cap's collision is strictly beyond the null parameter in every open case. The certificate includes all points, norms, pairings, supports, and rational parameters for replay.

## Constructive endpoint patterns

The next useful freedom is to give a new cap several finite nonright contacts with the retained support, especially contacts with two different H4 components after releasing an internal H4 wall. The checker enumerates every possible **spherical** one-node extension of each of the 40 eight-wall supports. There are 1,193 finite profile records and 56 separate dihedral families; these are support/profile counts, not distinct mirrors.

| Released global wall | H4 index | Finite profile records | Unbounded dihedral families |
| --- | ---: | ---: | ---: |
| 7 | 0 | 282 | 9 |
| 8 | 1 | 528 | 19 |
| 9 | 2 | 268 | 19 |
| 10 | 3 | 115 | 9 |

This enumeration has no finite-label cutoff. A spherical Coxeter diagram is a forest, so the new node contacts at most one node in each existing component. Its degree is at most three; degree-three contacts all have label 3. In a component of rank at least three, a degree-two attachment uses labels at most 5. A single attachment to an isolated A1 gives the separate family (I_2(m)) for **every integer (m\ge3)**. The exact Schur complement decides all the finite records. The isolated polar wall is normalized to norm two for this diagram calculation; physical pairings must be rescaled by its actual norm when constructing a cap.

There are 42 three-contact profiles touching two retained H4 nodes. All occur when wall 8 is released. Three concrete families sent to the joint cap search are:

| Edge record | Eight-wall support | Nonright contacts, all label 3 | Residual norm | Proposed endpoint type |
| --- | --- | --- | ---: | --- |
| 6 | 0,1,2,3,5,7,9,10 | 5,7,9 | (1/3) | H4 + D5 |
| 24 | 0,1,2,4,5,7,9,10 | 4,7,9 | (1/6) | H3 + E6 |
| 87 | 1,2,3,4,6,7,9,10 | 1,7,9 | (1/30) | E8 + A1 |

All other support contacts are right. Replacing wall 9 by wall 10 in these contact lists also works. In the last row, replacing wall 1 by wall 4 also works. These patterns are outside the previous single transverse weight-ray template. They are **local possibilities**, not certified globally compatible caps: the new normal must still satisfy its norm equation, the outgoing orientation, every other wall relation, and any proposed mutual cap relations.

Data keys for the next search are `outgoing_edges[edge_index]`, `H4_released_edge_profiles[str(edge_index)]`, and `complete_projective_ray_census.rays`. The last contains the complete 40-point initial gap set. Cutting those points is necessary when retaining this seed, but new intersections must still be checked; it is not a finite-volume certificate for a proposed enlarged chamber.

Generation and saved replay pass. A [separate source/proof review and fresh replay](fresh_h9_makarov_frontier_20261013_review.md) also pass; that review distinguishes the existing-wall collision census from any claim to find a first new mirror.

```sh
python scripts/fresh_h9_makarov_frontier_20261013.py --verify
```

This result completely describes the current cone and the finite endpoint profiles at its existing finite corners. It neither closes the arithmetic or geometric H9 problem nor excludes further caps or replacements. The dihedral tails describe finite endpoints; affine endpoint profiles are not included in the spherical-profile list.
