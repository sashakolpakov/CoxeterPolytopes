# H22: complete cusp-direction orbits and bounded inherited links

Adding paired walls over **1,065,708 primitive cusp rays** produces a compatible H22 chamber with **2,252,588 actual facets** and repairs the links at all those inherited cusps. The rays form 15 complete symmetry orbits represented in the saved cusp sample. This is not an exhaustive census of the old chamber's cusps.

A new exact strict spacelike witness proves that the resulting H22 chamber still has **infinite volume**, uniformly for every integer center label `m≥3`. No new finite-volume dimension is established.

## Why the preceding caps could not close the inherited cusps

The [previous construction](fresh_h22_borcherds_extension_20261013.md) retained the 210 old Borcherds H21 walls at added coordinate zero. Its 120,960 orbit caps and two center caps all had future timelike old-space projections. A future timelike vector has strictly negative pairing with every old future null cusp ray `δ`. Consequently none of those caps passes through `(δ,0)`. The old cusp link still has an unbounded new direction, regardless of how many caps from that timelike orbit are added.

The new move instead adds `(δ,+t)` and `(δ,−t)`. Each has square `t²`; the pair is parallel and is orthogonal to every old wall incident to `δ`. The old link has affine quotient rank 20, with positive kernel vectors in every component. The added parallel pair supplies an orthogonal affine A1 of rank one. Their product bounds the full 21-dimensional horospherical link. All other retained walls are strict at this cusp.

## Complete orbits and exact compatibility

The six stored integral isometries preserve the old form, strict center, and old wall set. The [producer](../scripts/fresh_h22_cusp_orbits_20261013.py) completes their orbits on 227 sampled old cusp rays. It obtains these 15 distinct orbits:

| Orbit | Rays | Height t |
|---|---:|---:|
| 1 | 3,780 | 1 |
| 2 | 5,040 | 1 |
| 3 | 6,048 | 1/4 |
| 4 | 60,480 | 1 |
| 5 | 20,160 | 1 |
| 6 | 280 | 1 |
| 7 | 30,240 | 1 |
| 8 | 2,240 | 1 |
| 9 | 30,240 | 1 |
| 10 | 120,960 | 1 |
| 11 | 60,480 | 1 |
| 12 | 120,960 | 1 |
| 13 | 241,920 | 1 |
| 14 | 120,960 | 1 |
| 15 | 241,920 | 1 |

The smaller height in orbit three repairs a forbidden angle with the preceding timelike orbit. Omitting that orbit is unnecessary.

The [standalone verifier](../scripts/fresh_h22_cusp_orbits_20261013_verify.py) proves each saved orbit is both closed under every generator and connected to its primitive representative, using an earlier-predecessor induction. Integer isometries have determinant ±1 and preserve primitivity. Global uniqueness therefore ensures that distinct saved rays are nonproportional.

Two distinct future primitive integral null rays have negative integral product, hence at most −1. For height-one walls this gives an integral same-sign pairing at most zero: right, parallel, or dotted. For every pair involving height one quarter, the nonzero pairings are strictly dotted. Opposite signs over the same ray give the required parallel pair. This proves all cusp/cusp relations without iterating through trillions of pairs.

Relations with the old 210 walls are exact integral checks, transferred from each representative by symmetry. For old timelike caps `(u/2,±sqrt(5/2))`, the representative maxima of `(u,δ)` are at most −6 for height one and at most −2 for height one quarter. These bounds prove strict dotted relations for both signs. The center-cap bounds use their added height at most four and squared norm at most two, so apply to every `m≥3`.

The full chamber has a strict timelike center and nonpositive pairings between distinct walls. Projecting the center onto any wall stays timelike and is strict for all other walls; every inequality is therefore an actual facet. The complete compatibility proof covers **2,537,075,222,578 unordered pairs**, using orbit transfer and the null-product lemma rather than explicit pair enumeration.

## Exact remaining escape and scope

The [certificate](fresh_h22_cusp_orbits_20261013_verify.json) stores an explicit rational vector of denominator `100000000` with

```text
square                         1656895868840709 / 10000000000000000 > 0
center pairing                −374305737 / 100000000
maximum new cusp-cap pairing  −59999 / 100000000 < 0
uniform center-cap bound      −152577893 / 100000000 < 0.
```

It checks every old, orbit, cusp, and center wall strictly. This witnesses an open spacelike patch and proves infinite volume for every center label under consideration. The new walls do cut the previous exact escape: 28 null directions have at least one signed wall that excludes it. That count is of directions, not of individual signed cutting facets.

Numerical constraint-generation linear programs proposed new escapes for labels 3, 4, 5, 6, 8, and 12. Only the separately verified rational witness supplies the uniform volume conclusion. The verifier imports no numerical scout, checks bounds preventing native-integer overflow, and anchors the preceding 121,172-wall compatibility certificate by source hashes.

```sh
python scripts/fresh_h22_cusp_orbits_20261013_verify.py --verify
```

An [independent source review and fresh exact replay](fresh_h22_cusp_orbits_20261013_review.md) passed, including the primitive-ray detail needed for the integral null-product bound. It reuses the same arithmetic implementation.

The [data manifest](../data/enumerations/fresh_h22_cusp_orbits_20261013.json) and [orbit arrays](../data/enumerations/fresh_h22_cusp_orbits_20261013.npz) specify every new wall. The completed 15 orbits are a constructive repair of inherited cusps, not a claim that all cusps have been enumerated, that all escape directions are removed, or that arbitrary further additions are impossible.
