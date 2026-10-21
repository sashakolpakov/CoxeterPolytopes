# Optional rational cutoff after the second cusp completion

The exact wall

```
n151 = 60x - (3/2)e - 59f - 15omega17 - 3omega31,
n151^2 = 1/4
```

extends the [151-wall chamber](fresh_h20_second_cusp_completion_20261016.md)
to **152 actual facets** and cuts its saved escape E5, with exact pairing
`29/200`. This is an optional branch. The 151-wall chamber remains the
base for coordinated finite-contact searches.

The [exact source](../scripts/fresh_h20_post_cusp_escape_cap_20261016.py)
and [certificate](fresh_h20_post_cusp_escape_cap_20261016.json) reconstruct
the physical normals and verify all 151 contacts: 16 are right, at the
support `S16`, and all other 135 are strictly dotted. No finite-angle
rounding or norm-2 assumption is involved. A common strict future center
and acute facet projections certify all 152 facets are actual. The five
proper cusps have strictly negative pairing with the new wall, so their
bounded links persist. All old finite vertices persist, and there are no
new finite vertices because the new wall has only 16 finite old
neighbors. The certified finite-vertex count remains at least 171.

Generation and a fresh saved replay pass. No independent review of this
optional branch is claimed here.

A strict escape remains:

```
E6 = x - e/100 - (49/50)f - (7/25)omega17 - (1/25)omega31
     - (1/100000)sum(other16 support dual roots),
E6^2 = 2500061/5000000000 > 0,
n151.E6 = -1/200.
```

All 152 pairings with E6 are strictly negative, and it lies in the
future chart of the common center. The optional chamber therefore
still has infinite volume.

## Whole residual-face diagnostic

The separate [numerical source](../scripts/fresh_h20_post_cusp_escape_cap_20261016_residual_scout.py)
and [saved diagnostic](fresh_h20_post_cusp_escape_cap_20261016_residual_scout.json)
project all old inequalities to the fixed Lorentz 5-space perpendicular
to S16. In the orthogonal basis consisting of the old corner C and
walls 141 through 144, write

```
X = C - sum(z_i*sqrt((-C^2)/q_i)*n_i).
```

The old cap inequalities give `z_i>=0`; other retained walls give linear
inequalities. The sign of `X^2` is the sign of `sum(z_i^2)-1`.
Floating-point Qhull gives the following diagnostic, with all squared
radii within `1e-6` of 1 left unclassified:

| Residual polytope | Vertices | Clearly inside | Near sphere | Clearly outside |
|---|---:|---:|---:|---:|
| Retained 151 walls | 150 | 13 | 19 | 118 |
| Optional 152 walls | 156 | 13 | 19 | 124 |

The new wall cuts six previously exterior vertices and introduces 12
new exterior vertices, whose squared radii range approximately from
1.00225 to 1.03954. Thus it removes the E5 region without completing a
new residual finite corner or reducing the open exterior frontier.
These are bounded numerical observations, not an exact cone census,
ideal-vertex classification or termination claim. They explain why this
cutoff is not imposed on the next coordinated search.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_post_cusp_escape_cap_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_post_cusp_escape_cap_20261016_residual_scout.py --verify
```
