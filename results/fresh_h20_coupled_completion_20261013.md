# Completing one H20 cusp and testing coupled cap orbits

This round reaches an exact **176-facet Coxeter chamber in H20**, but an exact strict spacelike escape still proves infinite volume. It does **not** produce a new finite-volume example. The constructive changes are a complete affine repair of the deficient cusp at the previous 98-wall stage, a cap with a deliberately dotted anchor, a coherent symmetry orbit of that cap, and a second coupled repair at the orbit's remaining gap.

The final cap also gives a useful exact obstruction to the next step: **no finite-volume completion can be obtained by adding only walls dotted to that cap**. It has just 18 intersecting old neighbors. The proof below is independent of any numerical orbit census and leaves intersecting additions, parallel additions, and wall replacements open.

| Exact stage | Change | Volume outcome |
|---|---|---|
| 104 walls | Six walls complete the new cusp to affine D8² × A1³ | The previous exact escape survives |
| 105 walls | A coupled D4 repair frees one anchor to be dotted | Cuts the previous escape; a new exact escape remains |
| 175 walls | Seventy further copies under a 72-element diagram-isometry group | A new exact escape remains after the whole selected orbit |
| 176 walls | A coupled repair with new finite labels 4 and 6 | Cuts the 175-wall escape; another exact escape remains |

All stages retain a strict timelike center and pass the complete new-versus-retained angle gate. Consequently every retained inequality is a true facet: orthogonally projecting the strict center onto one wall leaves a timelike point strict on every other wall, because every off-diagonal pairing is nonpositive. The analogous projection onto a positive-definite two-wall span proves that each finite angle is an actual ridge angle. No complete vertex or cusp census is claimed.

## Full affine completion of the deficient cusp

The input is the [exact 98-wall chamber](../data/enumerations/fresh_h20_gap_next_20261013_exact.json) and its [new-cusp diagnostic](../data/enumerations/fresh_h20_gap_next_20261013_cusp_diagnostic.json), described in the [previous report](fresh_h20_gap_next_20261013.md). The null ray is

```
delta = r18 + r96.
```

Its incident diagram consists of the affine A1 pair `{18,96}` and two spherical D8 components. In the following orders, divide old long roots by 3:

```
D8 component 1: [0,8,7,16,4,19,36,10]
D8 component 2: [2,14,13,26,11,20,5,12]
highest-root coefficients in either order: [2,1,2,2,2,2,1,1].
```

Let the resulting highest roots be theta1 and theta2. The exact construction also gives rational vectors p and q orthogonal to all 18 old incident normals, with p²=5, q²=1 and (p,q)=0. Their coordinates are reconstructed from the earlier recession certificate and an exact rank-20 kernel; no numerical vector is accepted as input. Add

```
-theta1 + 5 delta,   -theta2 + 5 delta,
 p + 5 delta,        -p + 5 delta,
 q + 5 delta,        -q + 5 delta.
```

Their local angles do not depend on the shift along delta. The shift 5 makes every new relation to an old nonincident wall strictly dotted. All 603 new pairs and all 642 signs at the 107 inherited cusp rays are checked exactly.

The repaired cusp now has 24 incident walls in five affine components, with ranks 8,8,1,1,1. Each component has an explicit strictly positive kernel whose normal sum is a positive multiple of delta. These relations, together with quotient rank 19, prove the horospherical link bounded. This is a full affine-link certificate, not a rank-only inference. Previously certified bounded links are preserved by acute-compatible additions.

The [104-wall certificate](../data/enumerations/fresh_h20_coupled_completion_20261013_exact.json) and [replay](../scripts/fresh_h20_coupled_completion_20261013_exact.py) retain the old strict escape of square `698631401/1887436800`. Closing this cusp alone does not close the chamber.

## A genuinely free dotted anchor

The saved escape polar P is normalized to

```
P²=97/16,
(P,r26/3)=(P,r28/3)=(P,r42/3)=-1.
```

Set S=(r26+r28+r42)/3 and C=(sqrt6/3)r49. The three summands of S are orthogonal norm-two roots; C is norm two and has pairing −1 with each. Thus these four old roots form a D4 block. Also

```
(P,C) = -17/4 + 11 sqrt6/12.
```

The new norm-two root is

```
r104 = 3P + S + b(C-P),
1.19 < b < 1.20.
```

The exact quadratic for b follows from r104²=2 and is archived with its positive-real branch. This gives π/3 against all three D4 leaves while allowing both wall23 and the central wall49 relation to be strictly dotted. The earlier finite-anchor scout missed this solution because it constrained the wall23 pairing to a finite label.

All 104 new pairs pass: 13 right, 8 of angle π/3, and 83 dotted. The new root cuts the exact old strict escape and is strict at all 107 inherited cusp rays. Nevertheless the [105-wall certificate](../data/enumerations/fresh_h20_coupled_completion_20261013_dotted_exact.json) contains another exact strict escape. Its [replay](../scripts/fresh_h20_coupled_completion_20261013_dotted_exact.py) reconstructs that escape from 20 exact active equations; the numerical scout supplies indices only.

## Testing a coherent symmetry orbit

Instead of repeatedly choosing a nearby coefficient, the next test transports r104 under a 72-element isometry group of the original H19 wall diagram. The [exact orbit replay](../scripts/fresh_h20_coupled_completion_20261013_orbit_exact.py) checks every listed permutation, group closure, the induced Lorentz isometry, and 72 distinct normal vectors. The cap-field radicand has rational norm `41700835/48`, which is not a rational square; this explicitly proves that the two-coefficient field encoding used for distinctness is injective. No theorem that these are all diagram automorphisms is required.

The original cap and 70 additional copies form a coherent extension of the retained chamber. All 9,765 new pair relations are checked directly, including relations to asymmetric earlier caps and the six cusp-completion walls; no symmetry of the entire 104-wall chamber is assumed. One orbit entry is unused. There is no orbit-maximality claim.

The [175-wall certificate](../data/enumerations/fresh_h20_coupled_completion_20261013_orbit_exact.json) includes a new strict spacelike escape, reconstructed in the cap's quadratic field from 20 independent active equations. Thus this selected symmetry orbit fails to close volume even though its mirrors are mutually compatible.

## A second coupled repair and the remaining gap

At that new gap, normalize the polar P so that

```
(P,r18)=-1,
(P,r30/3)=-3/2,
(P,r46/3)=-1.
```

The three displayed old roots are mutually orthogonal and have square 2. Put

```
u = r18 + (sqrt6/2)(r30/3) + r46/3,   u²=7,
r175 = u + tP,                     r175²=1,
3.53 < t < 3.54.
```

The selected positive branch solves `(P²)t²+2(P,u)t+6=0` exactly. The new root is expressed by adjoining one further square root to the previous cap's field; no minimal-degree claim is needed. Exact compositum arithmetic checks all 175 relations: 12 right, 4 of angle π/4, 2 of angle π/6, and 157 dotted. It cuts the exact 175-wall escape and remains strict at all 107 inherited cusp rays.

The [176-wall certificate](../data/enumerations/fresh_h20_coupled_completion_20261013_second_exact.json) and [replay](../scripts/fresh_h20_coupled_completion_20261013_second_exact.py) then certify a strict future spacelike point of square

```
1591730805825/162076295168 - (299252657/79138816) sqrt6 > 0.
```

This is the final exact checkpoint. It remains infinite volume.

## Why further dotted copies cannot complete this chamber

The [facet diagnostic](../data/enumerations/fresh_h20_coupled_completion_20261013_facet_obstruction.json), recomputed directly by [its replay](../scripts/fresh_h20_coupled_completion_20261013_facet_obstruction.py), finds exactly 18 old hyperplanes intersecting the last cap; the remaining 157 are strictly dotted. It verifies the projected norm signs as well as the pair labels.

Here is a direct obstruction valid for **any finite proposed addition of walls all dotted to r175**, assuming the resulting chamber is full dimensional. Work in W=r175-perpendicular, a 20-dimensional vector space of signature (19,1). Project a strict timelike point of the proposed completion onto W to get a timelike point c strict on every other inequality. This uses that all wall pairings are nonpositive; the point need not be the archived center.

There are at most 19 homogeneous equations requiring a vector w in W to be orthogonal to c and to the 18 intersecting-wall projections. Hence a nonzero such w exists, and w²>0 because c is timelike. The point

```
x = c + sqrt(-c²/w²) w
```

is future null and is strict on all 18 inequalities. A wall strictly dotted to r175 projects to a timelike normal in W. Its negative pairing with c fixes the future orientation, so it pairs strictly negatively with every future null x. This includes every new dotted wall in the proposed finite addition.

Move x slightly beyond the light cone within W, preserving all these strict inequalities, then slightly inward across r175. The resulting point is strict on every wall and spacelike, proving infinite volume. Thus the next successful construction must use at least one new wall intersecting or parallel to r175, replace retained walls, or change the chamber in another way. This theorem does not exclude those alternatives.

## Numerical scope and validation

The bounded numerical scouts are archived separately. The free-dotted-anchor scan produced [62 entries](../data/enumerations/fresh_h20_coupled_completion_20261013_dotted_scout.json); their numerical compatibility graph has no edges, and all leave a [common numerical escape](../data/enumerations/fresh_h20_coupled_completion_20261013_dotted_volume_scout.json). These are numerical observations, not an exhaustive real-parameter or label theorem.

After the second cap, its [72-copy numerical orbit](../data/enumerations/fresh_h20_coupled_completion_20261013_second_orbit_scout.json) appears mutually compatible with the 175-wall stage and still escapes at 247 walls. Those 71 further copies have not been promoted to an exact chamber. Numerically they are all dotted to r175; the exact obstruction above explains why such additions cannot repair its facet, without relying on those numerical pair signs.

The [review note](fresh_h20_coupled_completion_20261013_review.md) records a separate source/formula review and the precise replay scope. The deterministic exact commands are:

```sh
python scripts/fresh_h20_coupled_completion_20261013_exact.py --verify
python scripts/fresh_h20_coupled_completion_20261013_dotted_exact.py --verify
python scripts/fresh_h20_coupled_completion_20261013_orbit_exact.py --verify
python scripts/fresh_h20_coupled_completion_20261013_second_exact.py --verify
python scripts/fresh_h20_coupled_completion_20261013_facet_obstruction.py --verify
```

These proofs depend explicitly on the archived 98-wall certificate and its earlier exact dependencies. This round verifies 10,647 newly introduced pair relations; it does not rerun every old construction from scratch. No global novelty claim is made for the infinite-volume chambers, and no arbitrary H20 finite-volume completion is ruled out.
