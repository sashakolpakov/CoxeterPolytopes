# Releasing wall144: an E6 cutter with two finite corners

The frozen [156-wall chamber](fresh_h20_dominant_axis_cluster_20261016.md)
extends to **157 actual facets** by a norm2 wall that cuts E6. The new
wall156 is right to `S16+{95,151}`, has pi/3 contacts with142 and145,
and is strictly dotted to all136 other retained walls. In particular,
its contact with144 is dotted.

Exactly two new finite vertices are certified, bringing the inherited
lower bound to **182 finite vertices**. All five proper cusps survive.
A different exact strict spacelike escape E7 remains, so the chamber
still has infinite volume. The [source](../scripts/fresh_h20_release144_cap_20261016.py)
and [certificate](fresh_h20_release144_cap_20261016.json) passed generation
and a fresh saved replay. Independent review is being performed separately;
the primary source and data are frozen.

## Exact construction

Write `q=q142=309/200`, `a=n95.n145`, and `D=4-a^2<0`. Within
`U=span(n95,n142,n145,n151)`, the only nonzero off-diagonal product
is a; the norms are2,q,2,2. Set

```
p = (a/D)n95 - (2/D)n145 - n142/sqrt(2q),
k = omega31 - proj_U(omega31),
r156 = p - sqrt((2-p^2)/k^2) k.
```

The orthogonal projection is obtained from the displayed two-by-two
95/145 block and the two independent diagonal entries. Thus
`p^2=2/D+1/2`, `k.U=0`, and `r156^2=2` exactly. Certified intervals give
`p^2≈-.3075342124` and `k^2≈.8398940990>0`. Both radical branches are
tested; the negative branch is the unique one passing every retained
pair gate, the common-center condition and the E6 cut.

The prescribed products are

```
r156.n95 = r156.n151 = 0,
r156.n142 = -sqrt(q/2),
r156.n145 = -1.
```

Right relations with S16 follow its orthogonality to the residual
five-dimensional space. The remaining136 products are certified
negative with square greater than the corresponding norm product.
The E6 product is approximately **+.0170151388142**.

The inherited future timelike strict center H remains strict for all157
walls. For each wall n, `H-(H.n/n^2)n` is timelike and strictly inside
every other halfspace, proving that every facet remains actual. Old
finite vertices persist by the acute projection argument: their incident
positive-definite Z-matrix has nonnegative inverse, and the decomposition
of H along that vertex and its incident roots forces each new acute wall
to have strictly negative pairing with the old vertex. Five direct strict
negative cusp pairings preserve the five previously completed links.

## Complete finite vertices on the new facet

The only finite old neighbors are `S16+{95,142,145,151}`, twenty roots.
Of their nineteen-root subsets, only those omitting95 or145 are spherical.
Every other subset contains the dotted pair95/145. The two resulting
twenty-root diagrams are

```
A10+A4+A2+A3+A1,
A10+A4+A2+A2+2A1.
```

The source constructs both points by projection from H. Their norm and
future-chart sign are negative, and all157 wall pairings are checked;
only the twenty incident pairings vanish. The new-root Schur norms are
respectively1 and3/2. This proves exactly two new finite vertices, not a
complete census of every old finite or exterior ray. Their common support
`S16+{142,151,156}` gives an actual finite edge.

## The surviving escape and the reason to release144

Let C and the orthogonal residual roots141,...,144 be the frozen
[145-wall corner basis](fresh_h20_residual_h4_corner_20261016.md), with
`d=-C^2>0` and `q_i=n_i^2`. Define

```
z=(7/10,2/5,7/10,2/5),   epsilon=1/10^10,
E7=C-sum_i z_i sqrt(d/q_i)n_i-epsilon sum_{j in S16}omega_j.
```

Its exact squared norm is
`3d/10+epsilon^2 sum_{i,j in S16}J_ij>0`, approximately
`.0000310409162402`. Every one of its157 wall products is strictly
negative, and `E7.H<0`. The future ray from H toward E7 reaches a strict
open patch of the ideal boundary, certifying infinite volume.

The earlier [all-real residual gate](fresh_h20_residual_cutter_domain_20261016.md)
is relevant without adding any label assumption. The saved escape formulas
give `E6-E4=omega17/50=n17/100`. Acuteness at17 therefore implies
`r.E6<=r.E4`. Any S16-right acute E6 cutter retaining156 consequently
satisfies that gate's E4-cut hypothesis, and must be strictly dotted to144.
The present wall meets this necessary condition. Earlier right144 axis
caps and the separate three-axis bridge alternatives cannot remove E6.

## Bounded region comparison and a whole-face target

The [numerical companion](../scripts/fresh_h20_release144_cap_20261016_residual_scout.py)
and [saved data](fresh_h20_release144_cap_20261016_residual_scout.json)
compare only the three recipes with labels `(2,m,3,2)` on
`(95,142,145,151)`, for m=3,4,5. All three pass the numerical old gates.
With the same fixed C chart and2^20 scrambled Sobol directions, the
exterior Euclidean4-volume estimates are:

| Chamber | Estimated exterior volume |
|---|---:|
| Old156 | .0114417302544 |
| m=3, selected157 | .0114396643157 |
| m=4 alternative | .0114398849612 |
| m=5 alternative | .0114400543842 |

The selected wall has the best volume estimate of these three, but the
gain is only about.018%. The maximum squared residual vertex radius
remains1.417249978. The numerical vertex count changes176→181, with
22→24 inside,19 near the sphere unchanged, and135→138 outside.
These are numerical diagnostics with a1e-6 sphere band and no rigorous
integration error bounds; they are not hyperbolic volume estimates or
an exact ideal/exterior census.

A concrete next face is the actual H2 face with support
`S16+{151,156}`. Its numerical projected polygon has two finite corners
and three exterior vertices, on boundary pairs `(31,95)`, `(31,140)`
and `(140,145)`. Their squared radii are approximately1.001302,
1.002655 and1.002730. A coordinated truncation of this entire boundary
is a more specific target than another unrelated escape-only cut; no
completion of that face is asserted here.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_release144_cap_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_release144_cap_20261016_residual_scout.py --verify
```
