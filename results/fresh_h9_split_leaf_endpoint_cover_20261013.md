# A solver-free finite-volume exclusion for the fixed split-leaf cap bank

**No selection from the frozen 3,034-mirror bank can complete this seed to a finite-volume Coxeter polytope.** The result retains all sixteen split-leaf seed walls and concerns only this bank. It does not exclude other mirrors or changed seeds.

The bank has 3,040 indexed entries representing 3,034 distinct mirrors. The proof combines the [exact ideal-endpoint census](../data/enumerations/fresh_h9_split_leaf_edge_endpoints_20261013.json), the [certified compatibility supergraph](../data/enumerations/fresh_h9_split_leaf_3034_degree_graph_20261013.json), and one small integer spacelike witness. The graph allows parallel pairs and retains unresolved pair possibilities, so the conclusion covers finite volume as well as compactness.

Each of the two orthogonal sets of eight seed walls has two future ideal endpoints in its rank-two Lorentzian complement. If every selected cap pairs strictly negatively with one of these endpoints, moving slightly past that endpoint and then into the strict interior gives a spacelike escape. Consequently, a finite-volume selection must contain a cap pairing **nonnegatively** with each endpoint. Zero pairings are deliberately included. The endpoint certificate checks all 12,160 signs against the 3,040 exact normals; each of the four necessary cutter sets has fourteen entries, including two zero cases.

Exhausting the `14^4 = 38,416` choices gives exactly 36 distinct compatible minimal covers. They are precisely

```
{i,j},   0 <= i <= 5,   6 <= j <= 11.
```

Every compatible selection satisfying the endpoint conditions contains one of these pairs. Intersecting the two corresponding graph neighborhoods, with the selected entries themselves included, leaves exactly 28 possible cap entries for each pair. The union of all 36 allowed containers has only 332 entries. This is a necessary container calculation; the 28 entries need not be mutually compatible.

Let `u_0,...,u_7,v_0,...,v_7` be the ordered seed normals and let `phi=(1+sqrt(5))/2`. The single vector with integer coefficients

```
w = (8,8,8,1,-3,-3,-3,-3; -3,-3,-3,5,5,5,5,5)
```

in that seed basis has

```
(w,w) = (7242 - 2800*phi)/11 > 0,
(w,S) = -116 - 232*phi < 0,
(S,S) = -64 - 128*phi < 0,
```

where `S` is the sum of all sixteen seed normals. Every pairing of `w` with a seed normal is strictly negative. More decisively, every pairing with any of the **332** allowed cap entries is also strictly negative. Thus this same future spacelike vector survives every possible finite-volume candidate selection left by the endpoint conditions. Its strictness supplies an open escape region, proving infinite volume. No LP result or numerical volume classification enters this conclusion.

The [branch-and-witness certificate](fresh_h9_split_leaf_endpoint_cover_scout_20261013.json) reconstructs the cap covectors from the exact normals and seed Gram. Its standalone checker uses rational pairs for `Q(phi)` and exact signs of `a+b*sqrt(d)`, including dependent or zero radicands. It verifies the allowed normals' norm two, cap-seed acuteness, common strict center, all 332 witness pairings, and the seed Gram signature `(9,1,6)`. It also records an unnecessary but valid strict perturbation `w+S/8`, whose square is `(6912-3460*phi)/11 > 0`.

The endpoint and graph statements are separate certified premises, tied to the same normal archive by source hashes. This branch replay does not import either producer and makes no interval-sign decisions. The name `endpoint_cover_scout` is retained only from its initial assignment; the archived final calculation is exact. The universal witness was simplified from saved escape 120, but verification uses the displayed small integer vector directly.

All sixteen seed inequalities remain actual facets in any selected chamber: orthogonally projecting the common strict timelike point to one seed wall leaves every other inequality strict, because all pairings are nonpositive. Redundant optional caps can be removed. A finite-volume Coxeter chamber would therefore yield a compatible facet subset satisfying the four endpoint clauses, which the preceding argument excludes.

```sh
python scripts/fresh_h9_split_leaf_endpoint_cover_scout_20261013.py --verify
```

The endpoint verifier and compatibility-graph audit should be run alongside this short replay when rebuilding the complete proof. This is a finite-bank obstruction, not a theorem that arbitrary Coxeter caps cannot complete the seed.

The [independent final audit](fresh_h9_split_leaf_endpoint_independent_20261013.md) passes all 12,160 endpoint signs using the reversed radical-extension order, all 38,416 cover choices, the 36 allowed containers, and every pairing of the small witness with the 332-cap union and sixteen seed walls.
