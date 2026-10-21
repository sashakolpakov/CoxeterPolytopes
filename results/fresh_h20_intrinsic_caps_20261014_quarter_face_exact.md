# A selected quarter-face layer: 587 H20 facets

The [exact certificate](fresh_h20_intrinsic_caps_20261014_quarter_face_exact.json) adds a coherent selection of 60 caps to the frozen [527-wall golden chamber](fresh_h20_intrinsic_caps_20261014_golden_extension.md). It covers all 252 sampled quarter-vertex positions, retains every existing facet and all 252 certified finite vertices, and cuts the previous strict escape. **A new exact escape proves that the 587-wall chamber still has infinite volume.**

The selected recipes came from the [numerical face-center scout](fresh_h20_intrinsic_caps_20261014_quarter_face_scout.json). The [exact checker](../scripts/fresh_h20_intrinsic_caps_20261014_quarter_face_exact.py) reconstructs only these 60 selected recipes over `Q(phi)`, without using numerical predicates. It does not certify the scout's entire candidate bank or optimization completeness.

For each factor, take the affine-mark-weighted barycenter of the recipe's specified vertices. Sum the three factor barycenters and divide by4 to obtain `c`. The radius is `R=3/2`, so the cap normal is

```text
r=(c,(c²-9/4)/2,1),   r²=9/4.
```

The source wall order, factor marks, and conversion from factor-local recipe indices are explicitly bound and checked. All 33,390 new pair relations pass exactly:

| Relation | New pairs | Full 587-wall chamber |
|---|---:|---:|
| Right | 12 | 5,010 |
| pi/3 | 1,008 | 5,332 |
| pi/5 | 19 | 526 |
| Parallel | 24 | 24 |
| Dotted | 32,327 | 161,099 |
| Total | 33,390 | 171,991 |

Every new cap is strictly dotted to all 505 old spherical caps. Among the 60 new caps there are twelve right pairs and 1,758 dotted pairs. All 24 parallel relations are with old seed walls; no bounded-link or complete cusp-census claim is made for those new parallel pairs.

The common strict center remains `(0,9/2,1)`, of square−9. Every previously certified finite endpoint and strict facet point remains strict on the new caps. A new cap's own point

```text
X=(c,(c²+9/4)/2,1),   X²=-9/4,
```

is strict on every other wall. Thus all 587 inequalities define actual facets. These sixty strict facet points are not claimed to be vertices, and the 252 retained finite vertices are not a complete global census.

For every product-alcove vertex `z_v`, at least one selected ball contains `z_v/4`; the certificate verifies squared distance at most `9/4`, including any equality. This cuts every associated saved norm-one quarter-probe lift strictly. It does not establish a cover of the whole alcove.

For the first vertex, put `z=z_first/8`. The remaining strict escape is

```text
E=(z,(z²-1/2)/2,1),
E²=1/2,
(E,center)=-423/64-phi/32 < 0.
```

All 587 wall pairings are strictly negative. The positive square is deliberately `1/2`: a norm-one lift at this spatial point would fail the central-cap inequality.

## Why another isolated layer cannot close this chamber

The [ball-cover criterion](fresh_h20_intrinsic_caps_20261014_review.md) gives a necessary global condition. Every acute-compatible new spacelike normal retaining the 22 seed walls has a positive `U` coordinate: the three positive affine relations force a nonnegative coordinate, and equality would force its spatial part to vanish, leaving a null vector. It can therefore be scaled to a spherical-cap normal, and finite volume requires the corresponding closed balls to cover the compact product alcove.

The central radius-two ball has a boundary point strictly inside that alcove. Explicitly, on the first-vertex ray take `z*=2*z_first/sqrt(z_first²)`. Since `z_first²=302+4phi>4`, this is a strict convex combination of the alcove's interior origin and its vertex, hence is interior. Every other spherical cap in the 587-wall chamber is strictly dotted to the central cap, so its closed ball is strictly disjoint from the central ball.

Any finite collection of additional caps that remain strictly dotted to the central cap has the same property. Their finite positive separations leave an uncovered open neighborhood immediately outside the central sphere near `z*`. Thus **finite acute-compatible additions all strictly dotted to the retained central cap cannot yield finite volume**. This obstruction is conditional on retaining the seed and central facet with those pair relations; it allows finite-contact bridges, parallel contacts, or replacement. It is not an unrestricted no-completion result.

This is why covering a larger finite probe set is insufficient. The next constructive selection must connect this covering component or change the retained chamber, and must still pass the full global covering test.

```sh
python3 scripts/fresh_h20_intrinsic_caps_20261014_quarter_face_exact.py --verify
```

Generation and saved replay pass. A [separate matrix arithmetic audit](fresh_h20_intrinsic_caps_20261014_matrix_review.md) verifies all 171,991 pair relations, the facet points, the 252 finite endpoint Grams and inequalities, the strict center and the escape. Prior certificates remain unchanged.
