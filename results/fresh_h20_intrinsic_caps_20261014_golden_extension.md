# A coherent half-vertex layer: 527 H20 facets

The [exact extension certificate](fresh_h20_intrinsic_caps_20261014_golden_extension.json) adds 252 mutually compatible facets to the [275-wall golden chamber](fresh_h20_intrinsic_caps_20261014_golden_seed.md). It cuts that chamber's saved escape and preserves all 252 certified finite vertices. **The 527-wall chamber still has infinite volume**, certified by a new strict spacelike escape.

The construction retains the intrinsic affine `A11 × A2 × E6` seed with inradii `1,phi,1`, the radius-two central cap, and the 252 radius-one vertex caps. At every product-alcove vertex `z_v`, add the radius-one cap centered at `z_v/2`:

```text
r_v=(z_v/2,(z_v²/4-1)/2,1),   r_v²=1.
```

For each of its nineteen original incident seed walls, the new pairing is `-b_i/2`. It therefore has seventeen pi/3 contacts and two pi/5 contacts. Its other three seed relations are strictly dotted. Every new cap is dotted to the central cap, all old vertex caps, and every other new cap.

The [exact checker](../scripts/fresh_h20_intrinsic_caps_20261014_golden_extension.py) reuses the frozen seed's Fraction scalar implementation, binds its source and certificate hashes, and verifies all 100,926 new pair relations. The combined census is:

| Relation | New pairs | All 527 walls |
|---|---:|---:|
| Right | 0 | 4,998 |
| pi/3 | 4,284 | 4,324 |
| pi/5 | 504 | 507 |
| Dotted | 96,138 | 128,772 |
| Total | 100,926 | 138,601 |

The exact minimum squared separation of distinct half-vertex centers is `3+3phi>4`; the minimum squared distance from a half-vertex center to an old vertex center is `151/2+phi>4`. Thus there are no tangencies between the new hemispheres or between the two radius-one layers.

The common strict center remains `(0,9/2,1)`, of square−9. Every old finite endpoint is strict on all new caps, and the central facet point remains strict on them. Each new facet has its own point

```text
X_v=(z_v/2,(z_v²/4+1)/2,1),   X_v²=-1,
```

with exactly its own cap inequality zero and every other inequality strict. Consequently all 527 inequalities are actual facets. These new points prove facet actuality; they are not additional finite vertices. The report does not claim a complete corner census.

The first half-vertex cap cuts the prior norm-one escape with exact pairing1. At the quarter-vertex spatial point `z=z_first/4`, however,

```text
E=(z,(z²-1)/2,1),   E²=1,
(E,center)=-215/16-phi/8 < 0.
```

Every one of its 527 wall pairings is strictly negative. The global gap therefore remains despite the new coherent layer. Further progress must pass a global covering test; facet count and the termination of the original 252 edges alone do not supply one.

```sh
python3 scripts/fresh_h20_intrinsic_caps_20261014_golden_seed.py --verify
python3 scripts/fresh_h20_intrinsic_caps_20261014_golden_extension.py --verify
```

Both certificates are separate and frozen. Generation and saved replay pass, with an [independent geometry review](fresh_h20_intrinsic_caps_20261014_golden_review.md) and a [separate matrix arithmetic audit](fresh_h20_intrinsic_caps_20261014_matrix_review.md) covering this and the later 587-wall checkpoint. There is no finite-volume claim.
