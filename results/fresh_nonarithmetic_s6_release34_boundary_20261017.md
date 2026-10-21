# Boundary comparison for the released-34 construction

The [numerical ledger](../scripts/fresh_nonarithmetic_s6_release34_boundary_20261017.py)
uses one fixed Klein chart on the H3 face supported by
`S6=[0,1,2,3,5,6]`. For every pair of positive projected side normals it
clips their line by all the side halfspaces and the unit ball. This counts
actual hyperbolic edges, including lines with two open ends; it does not
count a wholly exterior projective segment as a hyperbolic edge.

| Quantity | Old 63 walls | Changed 65 walls |
|---|---:|---:|
| Finite S6 vertices | 19 | 23 |
| Compact geodesic edges | 22 | 28 |
| Edges with one open end | 13 | 13 |
| Edges with two open ends | 3 | 2 |
| Open ends, with multiplicity | 19 | 17 |

The three new S6 facets are complete polygons: A and B are quadrilaterals,
and C is a pentagon. The exact construction independently certifies those
three entire two-faces. Unlike the previous split pair, the changed-wall
construction gives a net reduction in this boundary diagnostic.

The reduction includes a cost. Deleting old 34 removes its five old edges
`7/34`, `8/34`, `12/34`, `34/47`, and `34/52`. The formerly compact edges
`7/12`, `7/47`, and `8/12` become open. On the other side, `8/44` and `30/44`
become compact, and `30/52` loses one of its two open ends. The `8/52` edge
gets a new finite endpoint on C but loses its old endpoint on deleted 34,
so it remains one-ended. Thus the count does not ignore the newly exposed
boundary.

Generation and fresh saved replay passed. The
[complete numerical records](fresh_nonarithmetic_s6_release34_boundary_20261017.json)
include all edge endpoints, finite supports, facet incidence and before/after
changes. This is a double-precision diagnostic, not the exact boundary proof
or a complete ambient cone certificate. The
[exact construction](fresh_nonarithmetic_s6_release34_exact_20261017.md)
retains a strict spacelike escape and therefore has infinite volume.
