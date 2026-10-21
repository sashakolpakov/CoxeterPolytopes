# A finite corner at the end of an edge from the completed cusp

An exact norm-2 cap extends the [151-wall chamber](fresh_h20_second_cusp_completion_20261016.md)
to an **alternative 152-wall chamber**. It cuts E5 and creates exactly
one new finite vertex, joined by an actual edge to the completed cusp
`nu2`. The five proper cusps persist. A new strict spacelike escape
remains, so volume is still infinite.

This branch does **not** retain the
[optional all-dotted cutoff](fresh_h20_post_cusp_escape_cap_20261016.md).
The [source](../scripts/fresh_h20_second_cusp_pi8_cap_20261016.py) and
[certificate](fresh_h20_second_cusp_pi8_cap_20261016.json) pass generation
and a fresh saved replay. Independent review of this new cap is pending.
No frozen predecessor is changed.

## Exact cap

The new wall 151 is right to S16 and to walls 95, 142 and 145. It meets
wall 150 at pi/8. Every other old contact is strictly dotted. The exact
recipe uses only projections and positive square roots.

Let `U=(n95,n142,n145)`, `a=n95.n145<-2`, and `q=n142^2=309/200`.
The Gram matrix and its inverse are

```
G_U = [[2,0,a],[0,q,0],[a,0,2]],
G_U^-1 = [[2/(4-a^2),0,-a/(4-a^2)],
          [0,1/q,0],
          [-a/(4-a^2),0,2/(4-a^2)]].
```

It is a nonsingular Lorentzian 3-space inside the fixed residual
Lorentz 5-space. Write `proj_U` for its exact orthogonal projection and
define

```
v = n150 - proj_U(n150),
k0 = omega31 - proj_U(omega31),
k = k0 - (k0.v/v^2)*v,
t = -sqrt(2*n150^2)*sqrt(2+sqrt(2))/2,
r = (t/v^2)*v - sqrt((2-t^2/v^2)/k^2)*k.
```

The certificate proves `v^2>0`, `k^2>0`, and `2-t^2/v^2>0`.
By the projection identities, `r^2=2`, r is right to U, and
`r.n150=t`. The negative square-root branch passes all old gates and
cuts E5; the other branch is also tested and rejected. The residual
chart is approximately

```
(x,e,f,omega17,omega31):
(26.1215180208, 0, -25.9143969961, -4.6425694395, -2.0005032057).
```

All 151 new-old checks pass: **19 right, one pi/8, and 131 dotted**.
The least robust dotted contact in the numerical discovery was wall 31;
the primary certificate verifies its strict inequality with rational
intervals, along with every other dotted contact. Exact finite-angle
and right identities follow from the displayed projections, not from
intervals containing the target values.

The retained center is future timelike and strictly inside all 152
halfspaces. Acute projections from that center verify that every facet
is actual. All five previous proper cusp rays pair strictly negatively
with r, so their links are preserved.

## Complete new finite-vertex census

The only 20 finite old neighbors are `S16+{95,142,145,150}`. A finite H20
Coxeter vertex containing r must choose 19 of these old neighbors.
Among the 20 possible supports, exactly one is positive definite:

```
S16 + {142,145,150,r},
diagram A10 + A4 + A2 + 2A1 + I2(8).
```

The other 19 supports contain either the dotted pair 95/145 or the
dotted pair 95/150. This is a complete census of the new finite
vertices, with no search cutoff. The positive determinant of the final
two-node block is `2*n150^2-t^2>0`. Projection of the common center onto
the intersection of the 20 incident hyperplanes gives an exact future
timelike point; its complete 152-wall pairing list confirms the actual
finite corner. All previous finite vertices persist, bringing the
certified total to **at least 172**. The old ambient census is still
partial.

The 19 shared roots `S16+{142,145,150}` are independent and positive
definite. They vanish at both this new finite vertex and `nu2`. Every
other wall is strictly negative on the open segment between the two
endpoints. Hence this segment is an actual finite-to-ideal edge. The
construction completes the finite end of an edge from an existing
proper cusp, rather than merely cutting an isolated open region.

## Escape and whole residual-face diagnostic

The E5 cut is strictly positive, approximately `0.0109865154`. A strict
escape remains:

```
E6pi8 = x - e/100 - (49/50)f - (7/25)omega17 - (1/20)omega31
        - (1/100000)sum(other16 support dual roots).
```

Its exact squared norm is positive, every wall pairing is strictly
negative, and its pairing with the common future center is negative.
The certificate records all these signs; they prove infinite volume.

The separate [numerical residual source](../scripts/fresh_h20_second_cusp_pi8_cap_20261016_residual_scout.py)
and [saved diagnostic](fresh_h20_second_cusp_pi8_cap_20261016_residual_scout.json)
compare all projected inequalities in the four-dimensional residual
affine chart. Qhull reports:

| Residual polytope | Vertices | Clearly inside | Within 1e-6 of sphere | Clearly outside |
|---|---:|---:|---:|---:|
| 151-wall base | 150 | 13 | 19 | 118 |
| pi/8 alternative | 154 | 14 | 19 | 121 |

The cap cuts eight exterior vertices and introduces one finite and
eleven exterior vertices. Some projected incidences are nearly
degenerate, so these floating-point counts are explicitly diagnostic,
not an exact cone or ideal-vertex census. The exact new finite-corner
result above is independent of this computation. The remaining
exterior frontier has not closed; no finite-volume or new-dimension
claim is made.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_second_cusp_pi8_cap_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_second_cusp_pi8_cap_20261016_residual_scout.py --verify
```
