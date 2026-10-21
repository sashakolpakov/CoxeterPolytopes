# Complete edge census at the repaired H20 cusp

The frozen 95-wall chamber has **1,728 edges issuing from the newly completed cusp**. Their opposite ends are now classified exactly:

| Opposite end | Number |
|---|---:|
| Finite vertex | 26 |
| Ideal contact with additional walls, with an unbounded link | 16 |
| Light-cone limit with exactly the 19 positive-definite support walls active | 1,686 |

The 26 finite endpoints are distinct. Six were already among the twelve endpoints recorded in the 95-wall construction; twenty are new. Thus the chamber has **at least 116 verified finite vertices**, together with 26 certified edges from this proper cusp to finite vertices. This is a complete census of the edges from this one cusp, not a complete vertex or edge census of the chamber. No walls are added here.

The original rational strict escape still satisfies every one of the unchanged 95 inequalities, with square `259854917001622621/62500000000000000`. Consequently the chamber still has infinite volume. The new census describes the remaining geometry rather than claiming global closure.

## Exact construction of all cusp edges

The [95-wall certificate](fresh_h20_e7_edge_cusp_repair_20261014_exact.json) records the full incident link as six mutually orthogonal affine components of sizes `9,2,6,2,4,2`, with strictly positive null relations and quotient ranks `8,1,5,1,3,1`. Its bounded horospherical chamber is the product of the six corresponding affine simplices. Every horospherical vertex omits exactly one wall from each component, giving

`9 × 2 × 6 × 2 × 4 × 2 = 1,728`

positive-definite supports of rank 19. These are all the edges from this cusp. The root rescalings used to make component Gram matrices rational are positive and do not change their hyperplanes or inequalities.

Let `ν` be the saved future null cusp point and `H` the saved strict timelike center. For a support `S`, put

`W = H − projection_S(H)`.

Then `W² < 0`, `W·ν = H·ν = h < 0`, and `W` is orthogonal to all 19 support roots. The omitted cusp walls have strictly negative pairing with `W`, as is also checked directly. Their pairings are constant along the line

`X(s) = W + sν`.

All nonincident walls `n` have `ν·n < 0`. Their inequalities are precisely

`s ≥ (W·n)/(−ν·n)`.

Thus the largest of these thresholds is the first contact when walking away from `ν`. The opposite light-cone parameter is `s_null = −W²/(2h)`. Parameters greater than `s_null` give future timelike points.

Only 25 small projections are needed: for each affine component and each omitted node, the checker projects `H` onto the remaining positive-definite roots. Orthogonality of the components gives the combined projection, its square, and all 95 pairings by addition. It does not invert 1,728 separate rank-19 matrices.

The checker uses floating-point evaluation only to propose a maximizing wall. It then compares that threshold with every other threshold exactly. For the winning wall, write `D = −ν·n > 0` and `T = W·n`. The positively scaled first-contact point is

`X = D W + Tν`.

All 95 pairings, its square, and its future orientation are checked in the exact field. This avoids division even when the winning wall is the bridge carrying the additional quadratic radical. The bridge wins seven algebraic thresholds, all lying beyond the opposite light cone; none produces a hyperbolic endpoint.

## Finite endpoints

Every finite endpoint has exactly 20 active walls: the 19 positive-definite support roots and one first-contact wall. The twentieth wall is independent because it pairs nontrivially with `ν`. Its Gram matrix is therefore positive definite in the orthogonal complement of the explicitly timelike endpoint. Every other wall is strict there, so these are actual finite vertices and the intervening segments are actual edges.

Wall 8 supplies 24 finite endpoints; walls 0 and 10 supply one each. Full incident sets establish distinctness and identify the six overlaps with the previously verified 96 vertices. Exact coordinates, norms, future pairings, and active sets are saved for all 26 points. The checker does not infer a complete chamber census from this list.

## Ideal contacts and open opposite ends

The sixteen first contacts occurring exactly on the light cone each have twenty active roots. The separate [ideal-link diagnostic](fresh_h20_e7_repaired_cusp_walk_20261014_ideals.json) reconstructs their complete active Gram matrices. Every quotient has rank 19, but every link also has at least one spherical component. Quotient rank alone therefore does not establish boundedness.

For each of these sixteen points the diagnostic constructs a negative fundamental weight in a spherical component. It is orthogonal to the null point, has strictly positive square, and pairs nonpositively with every incident root, with at least one strict inequality. This is an explicit nonzero horospherical recession direction. None of the sixteen contacts is an ideal vertex with bounded link.

For each of the other 1,686 directions, the checker verifies that the opposite null point

`Z = (−2h)W + W²ν`

has exactly the 19 support roots active and is strict against all remaining walls. It is future, since `Z·H = −h W² < 0`, and null. The active Gram matrix is positive definite. Its nineteen independent Euclidean quotient normals admit a recession direction, so this ideal-boundary end also has an unbounded link. No omitted or nonincident wall becomes active there.

The stored first-contact ties in these open cases concern intersections beyond the light cone. They are not counted as finite or ideal vertices. All 26 finite endpoints and all 16 ideal contacts have a unique additional active wall.

## Certificates and verification

- [Main exact checker](../scripts/fresh_h20_e7_repaired_cusp_walk_20261014.py) and [census JSON](fresh_h20_e7_repaired_cusp_walk_20261014.json): all 1,728 supports and exact first-contact comparisons against all 95 walls.
- [Ideal-link checker](../scripts/fresh_h20_e7_repaired_cusp_walk_20261014_ideals.py) and [ideal-link JSON](fresh_h20_e7_repaired_cusp_walk_20261014_ideals.json): all sixteen full active sets, quotient ranks, affine/spherical decompositions, and explicit positive-norm recessions.
- [Frozen predecessor report](fresh_h20_e7_edge_cusp_repair_20261014.md) and its [independent review](fresh_h20_e7_edge_cusp_repair_20261014_review.md).

Both producers and both saved-result replays passed. Their replay commands are:

```sh
PYTHONPATH=src:. python scripts/fresh_h20_e7_repaired_cusp_walk_20261014.py --verify
PYTHONPATH=src:. python scripts/fresh_h20_e7_repaired_cusp_walk_20261014_ideals.py --verify
```

The source manifests bind the frozen predecessor data and arithmetic helpers. The main JSON stores coordinates only for actual finite/ideal contact points, and uses compact threshold recipes for open directions. This work neither repeats nor claims a complete 95-wall finite-support census.
