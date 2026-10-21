# Boundary effect of the two-cap split bridge

The exact pair separates the former component-B exterior path into two patches. It closes the old edge `[30,44]` and supplies the first finite endpoint of `[30,52]`. It also creates three new open edges. The chamber remains infinite-volume, with the same strict escape.

The [ledger checker](../scripts/fresh_nonarithmetic_s6_split_pair_20261017_edge_ledger.py) reads the [complete 65-wall cone](fresh_nonarithmetic_s6_split_pair_20261017_cone.md) and its frozen 63-wall predecessor. Its [saved JSON](fresh_nonarithmetic_s6_split_pair_20261017_edge_ledger.json) contains every finite-incident edge, every projective S6 edge, exact causal classifications, and the exterior connected components.

| Quantity | 63 walls | 65 walls |
|---|---:|---:|
| Ambient finite vertices | 193 | 196 |
| Ambient finite-incident compact edges | 431 | 434 |
| Ambient finite-incident open edges | 875 | 896 |
| S6 finite vertices | 19 | 22 |
| S6 exterior vertices | 11 | 12 |
| S6 compact geodesic edges | 22 | 25 |
| S6 edges with one open end | 13 | 16 |
| S6 edges with two open ends | 3 | 2 |
| S6 projective edges wholly outside hyperbolic space | 7 | 8 |

Here `S6=[0,1,2,3,5,6]`. In the following supports, the listed side indices are adjoined to S6.

The exact changes to actual geodesic edges are:

- `[30,44]`: one open end becomes a compact edge.
- `[30,52]`: two open ends become one open end.
- `[30,63]` and `[30,64]`: new compact edges.
- `[44,63]`, `[52,64]`, and `[63,64]`: new one-open-end edges.

No old geodesic edge support disappears. The total number of open S6 geodesic ends, counting two for a two-open edge, changes from 19 to 20. The pair supplies a coherent finite bridge but does not establish net closure by that measure.

The former four-node exterior path

```
[8,44,52] -- [30,44,52] -- [30,47,52] -- [34,47,52]
```

is replaced by two components:

```
[8,44,63] -- [8,63,64] -- [8,52,64]

[30,47,52] -- [34,47,52].
```

All three displayed internal links in the new components lie wholly outside hyperbolic space. The separate seven-vertex component A is unchanged. This gives two distinct remaining B patches, centered on old walls 8 and 47, for coordinated continuation.

The side-30 face now has eight finite and three exterior vertices. Each new S6 facet has two finite and two exterior vertices. The previously compact tetrahedral, quadrilateral, and prism faces remain intact.

All finite vertices and all S6 rays are simple. The complete support index therefore identifies each relevant rank-eight edge face and its two extreme rays. For two exterior endpoints, the checker uses the exact two-by-two Gram determinant and mutual product; whenever their positive span meets hyperbolic space, it supplies an explicit future timelike positive combination. This distinguishes a genuine edge with two open ends from an edge of the projective cone entirely outside hyperbolic space. The ambient enumeration is of finite-incident edges; it does not claim a complete edge list for nonsimple exterior ambient vertices.

Generation and fresh saved replay both passed. This is an exact incidence calculation from primary cone certificates, not an independent arithmetic reconstruction of those cones.
