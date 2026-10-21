# Independent review of the 327- and 328-wall central-ridge constructions

Date: 2026-10-14. **PASS; both chambers have infinite volume.**

The [independent checker](../scripts/fresh_h20_central_ridge_20261014_review.py)
imports neither producer. It uses the basis
`1,sqrt(5),S,sqrt(5)S`, with `S²=406/9−9sqrt(5)` and `S>0`,
and converts the saved `Q(phi,alpha)` coordinates before comparison.
All sign decisions are exact rational comparisons in this field.
The [saved review certificate](fresh_h20_central_ridge_20261014_review.json)
binds both original scripts and certificates by SHA-256.

The review checks the 288 cap norms, all 41,328 cap pairs using literal
vector products for each coordinate-permutation class, every new/core
relation, and all 4,608 finite-vertex/cap pairings. It reconstructs the
additional mirror and its new finite endpoint independently, verifies all
327 new pair relations, and checks that the original sixteen finite
vertices survive. A common strict timelike point and nonpositive wall
pairings establish that the walls are actual facets.

The exact finite-neighbor graphs give sixteen and seventeen finite vertices,
respectively. The counts of elliptic rank-19 faces are 524,608 and 524,627.
The extension closes one finite edge and leaves 338 finite-ended edges
without a second finite vertex or a bounded-link ideal vertex. In addition,
524,288 seed-only vertical edges are outside the finite-vertex components.
Their box-vertex ideal ends have unbounded links. “Open edge” here does not
mean that its hyperbolic ray lacks a point on the ideal boundary.

Both chambers retain the same strict norm-one escaping vector, so neither
is a finite-volume example. The
[noncrystallographic face obstruction](fresh_noncrystallographic_face_gate_20261014.md)
also excludes finite-volume completion within the specified right-angled
ridge scheme. More iterations of that scheme are not a research priority.

Reproduce with:

```sh
python3 scripts/fresh_h20_central_ridge_20261014_review.py --verify
```

See the [construction and edge-step report](fresh_h20_central_ridge_20261014.md)
for the formulas and the precise family retained.
