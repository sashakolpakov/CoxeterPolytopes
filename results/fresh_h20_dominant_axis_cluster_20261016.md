# A compatible pair cutting the two remaining long residual axes

Two coordinated caps extend the
[154-wall dominant-axis chamber](fresh_h20_dominant142_cap_20261016.md)
to **156 actual facets**. Each has pi/3 contact with wall153; the two new
walls are mutually dotted. Exactly four new finite vertices are certified,
so there are now **at least180 finite vertices**. All five proper cusps
and the old compact triangle survive.

The same strict E6 survives, proving infinite volume. The construction
addresses the large141/143 regions rather than that escape. In the fixed
residual chart, a numerical common-direction score estimates about99.11%
less exterior Euclidean volume than the153-wall chamber, and reduces the
maximum squared vertex radius from16.3346 to1.41725. Neither statement
asserts hyperbolic finite volume.

The [source](../scripts/fresh_h20_dominant_axis_cluster_20261016.py) and
[certificate](fresh_h20_dominant_axis_cluster_20261016.json) pass generation
and a fresh saved replay. Independent review is pending. Every earlier
package remains frozen.

## Exact compatible recipes

Let `e_i=n_i/sqrt(q_i)` for the four mutually orthogonal residual roots
141,...,144, and `e0=C/sqrt(-C^2)`. Thus `e_i^2=1`, `e0^2=-1`.
The previous cap is

```
r153 = sqrt2*(-c153*e142+s153*e0),
c153^2-s153^2=1,    s153=.4764935409119465...
```

Set

```
s_star = 1/(2*s153),
c_star = sqrt(1+s_star^2),
r154 = sqrt2*(-c_star*e141+s_star*e0),
r155 = sqrt2*(-c_star*e143+s_star*e0).
```

Both new roots have squared norm2. Their products with r153 equal
`-2*s_star*s153=-1`, giving pi/3. Since `s_star>1`, their mutual product
is `-2*s_star^2<-2`, strictly dotted. Numerically

```
s_star=1.0493321673218596,
c_star=1.4495164701983870,
r154.r155=-2.2021959947527825.
```

Root154 is right to `S16+{142,143,144}`; root155 is right to
`S16+{141,142,144}`. Every other old contact except the stated pi/3
contacts is strictly dotted. All309 new pair checks give **38 right,
2 pi/3 and269 dotted** relations.

The inherited strict future center is checked against all156 walls;
its acute projections certify that all facets are actual. All ten
new-root/old-proper-cusp pairings are strictly negative. The compact
triangle persists because its three old finite vertices remain feasible.

## Complete new finite vertices

Each new root has exactly20 finite old neighbors. For root154 these are
`S16+{142,143,144,153}`; for root155 they are
`S16+{141,142,144,153}`. In either case the only successful19-old-root
subsets omit142 or omit153. The resulting diagrams are

```
A10+A4+A2+2A1+A2,     and     A10+A4+A2+4A1.
```

Each projected point is future timelike and has all156 pairings checked,
including the other new cap. Every other subset contains the dotted
pair142/153 and is not spherical. The new roots are mutually dotted,
so no finite vertex can contain both. Thus exactly four finite vertices
are added; this does not replace the inherited partial ambient census.

## Strict escape and whole-region diagnostic

The unchanged E6 has exact norm `4750061/5000000000>0`, negative
future-chart sign, and strictly negative pairing with all156 walls.
Its two new pairings are approximately `-1.24305183` and `-1.47692473`.
This certifies that the ambient chamber remains infinite-volume.

The new inequalities in the fixed C chart are simply

```
z141 <= s_star/c_star = .7239187611150383...,
z143 <= s_star/c_star = .7239187611150383...,
```

in addition to the prior `z142<=.4301565885384587...`.
The [numerical source](../scripts/fresh_h20_dominant_axis_cluster_20261016_residual_scout.py)
and [saved score](fresh_h20_dominant_axis_cluster_20261016_residual_scout.json)
use every retained inequality, the same fixed chart and the same2^20
Sobol directions at all stages:

| Quantity | 153 walls | 154 walls | 156 walls |
|---|---:|---:|---:|
| Maximum squared vertex radius | 16.3345753 | 6.98324997 | 1.41724998 |
| Estimated exterior angular fraction | .961934 | .492572 | .224974 |
| Estimated exterior Euclidean4-volume | 1.289858 | .305608 | .0114417 |

The last estimate is99.1129% below the153 baseline. These quadrature
estimates have no rigorous error bounds and are not hyperbolic volumes.
The numerical residual vertex census changes from168 to176 vertices:
18 to22 finite,19 near the unit sphere unchanged, and131 to135 exterior.
The pair cuts eight old exterior vertices, creating four finite and
twelve exterior vertices. Counts use a1e-6 sphere band and are not an
exact cone completeness or ideal classification claim.

The largest remaining exterior vertex is approximately

```
(z141,z142,z143,z144)=(.723919,.430157,.723919,.429067),
active projected walls={133,153,154,155}.
```

This is a concrete joint-region target. No closure of it is claimed in
this package, and the strict escape remains an independent obstruction.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_dominant_axis_cluster_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_dominant_axis_cluster_20261016_residual_scout.py --verify
```
