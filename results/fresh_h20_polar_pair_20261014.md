# Two rational polar walls and sixteen new finite vertices

The [exact construction](fresh_h20_polar_pair_20261014.json) adds two mutually orthogonal rational walls to the frozen [98-wall chamber](fresh_h20_bridge_edge_closure_20261014.md). All 100 listed walls are actual facets. The three previously bounded cusps survive, and the saved strict spacelike escape survives as well: **this chamber still has infinite volume**.

The separate [finite-support census](fresh_h20_polar_pair_20261014_census.json) finds exactly 16 finite vertices incident to at least one of the new walls. Fourteen lie on each new wall and twelve lie on both. Together with the 118 previously certified vertices, this gives at least 134 known finite vertices. This is a complete census on the two new walls, not a complete census of all old vertices.

## Exact normals and compatibility

In the inherited 21-coordinate rational form, the new normals are

```
a = [7,14,21,18,15,12,9,6,3,1,-1,-3,-5,-6,-7,-8,-9,-10,10,0,0]
b = [14,27,40,35,30,25,20,15,10,5,0,-5,-10,-12,-14,-16,-18,-20,20,0,0]
```

Their squares are `a²=2`, `b²=4`, and `a·b=0`. Their indices in the enlarged chamber are 98 and 99. The 197 new pair checks give:

| Relation | Count |
|---|---:|
| Right | 53 |
| π/3 | 4 |
| π/4 | 3 |
| Parallel | 12 |
| Strictly dotted | 125 |

All finite-angle equalities are coefficient identities. Strict signs and dotted inequalities use outward rational intervals in the inherited radical tower; 192-bit enclosures suffice. The recipes are literal rational vectors, so the numerical polar scout is not a dependency of this certificate.

Each vector satisfies all old 98 halfspaces and has positive square. Its old-wall zero set has rank 20, checked independently, so it is an old spacelike extreme ray. Its own new wall excludes that direction. The certificate saves these zero sets and the complete finite-angle neighbor lists used by the separate vertex census.

## Interior, facets, and retained geometry

Let `T` be the sum of the first 95 physical normals. The hash-bound independent whole-face certificate establishes that `T` is future timelike and strict at those 95 walls. Summing the subsequent certified pairings establishes strictness at walls 95–97; direct rational-vector pairings establish strictness at the two new walls. Thus `T` is a common strict center for all 100 walls.

For every wall normal `N`, the projection

`T − (T·N/N²)N`

is timelike, lies on that wall, and is strict at all other walls by acuteness. The [acute-extension preservation lemma](acute_extension_preservation_20261014.md) therefore proves actuality of all facets and strict preservation of every old finite vertex. The primary certificate additionally saves both new facet points and checks all their 100 signs. It checks the two new pairings at each of the three old bounded cusps directly; all six are strict, so their completed links are unchanged.

The rational census exhausts the complete finite-neighbor graphs on walls 98 and 99. It tests 20,810 distinct twenty-wall supports with exact componentwise positive-definiteness checks and obtains 16 positive supports with distinct primitive rational timelike points. The complete acute Gram and common center turn these supports into actual finite vertices. The old finite vertices are strict at both new walls, so none of these sixteen is an old vertex.

## The remaining global escape

The inherited rational escape `E` has

`E² = 20181174579467302741 / 10^18 > 0`.

Its new pairings are exactly

```
E·a = −361751915426463677 / 500000000
E·b = −1586702382310251073 / 1000000000
```

Both are strict. Together with the frozen old 98 signs and the future-sheet check, this certifies a strict spacelike escape for all 100 walls. The construction supplies new finite intersections and polar truncations, but does not establish global closure.

There is a second exact limitation. Both new walls are strictly dotted to wall 75. For the retained S75 support, the inverse Gram is entrywise nonnegative and its wall-75 diagonal is 5, as checked in the [S75 projection gate](fresh_h20_s75_projection_gate_20261014.md). If either new normal has square `q`, its support pairing at 75 has square greater than `2q`; hence its S75 projection has square strictly less than `q−5(2q)=−9q`. The projected common center fixes the future orientation. Both new walls consequently remain strict throughout the entire S75 hyperbolic plane, which is still unbounded. The sixteen new vertices occur elsewhere and do not provide a boundary for that plane.

## Reproduction and scope

```
PYTHONPATH=src:. python scripts/fresh_h20_polar_pair_20261014.py --verify
PYTHONPATH=src:. python scripts/fresh_h20_polar_pair_20261014_census.py --verify
```

The [primary checker](../scripts/fresh_h20_polar_pair_20261014.py) passed generation and saved replay. It inherits the old mutual-pair and old-vertex certificates, checks every new pair, and proves the actual-facet, cusp, and escape assertions above. The [separate rational census](../scripts/fresh_h20_polar_pair_20261014_census.py) owns the complete new-wall vertex enumeration.

The [independent arithmetic and geometry review](fresh_h20_polar_pair_20261014_review.md) also passed generation and saved replay. It independently checks all 197 pairs and reconstructs the finite-support census using fixed-size clique recursion and exact tree/cycle Schur tests. Direct evaluation gives all 1,600 physical signs at the sixteen new vertices. It also checks the two extreme-ray ranks, the new facet points, the cusp signs and the escape. The old 98-wall census remains an inherited dependency. Frozen predecessor files were not changed.
