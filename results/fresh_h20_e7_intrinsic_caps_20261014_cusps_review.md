# Independent review of the two E7 cusp fillers

**2026-10-14 — PASS.** The two shift-2 caps complete both specified ambient cusp links of the [85-wall replacement](fresh_h20_e7_intrinsic_caps_20261014.md). All 87 facets are actual, all 44 finite vertices of its H13 face survive, and that complete face is unchanged. The saved strict spacelike escape survives, so this 87-wall H20 branch still has infinite volume.

The [separate checker](../scripts/fresh_h20_e7_intrinsic_caps_20261014_cusps_review.py) and [review certificate](fresh_h20_e7_intrinsic_caps_20261014_cusps_review.json) bind the [primary source](../scripts/fresh_h20_e7_intrinsic_caps_20261014_cusps.py), [primary certificate](fresh_h20_e7_intrinsic_caps_20261014_cusps.json), and frozen 85-wall input. They do not import producer or project arithmetic. Sparse coefficients in `Q(sqrt(2),sqrt(79),sqrt(11))`, Fraction matrix elimination, and outward 192-bit rational radical bounds provide the independent arithmetic. The [construction report](fresh_h20_e7_intrinsic_caps_20261014_cusps.md) gives the geometric context.

The review reconstructs the E7 Gram inverse from the physical retained normals. Its unique diagonal entry two is at wall 15. The associated fundamental weight has norm two; the two future null face rays are orthogonal to it. Thus the normals `−ω15 + 2ν1` and `−ω15 + 2ν2` have norm two. Direct checking reproduces all 171 new relations: 39 right, two π/3, and 130 strictly dotted. The common shift two is a verified choice, with no claim that it is the smallest real admissible shift.

The checker independently projects the strict center to each of the 87 wall hyperplanes and verifies all 7,569 resulting inequalities. Each point is future timelike and strict against the other 86 walls. It also checks all 3,828 inequalities at the 44 inherited finite face vertices and reconstructs their positive-definite rank-20 incident Grams. Both new caps are strict at every one of these vertices.

For each ideal ray, all 87 signs are recomputed. The checker derives the connected components from the complete physical incident Gram, without consuming the producer's saved affine partition or kernels. It infers rationalizing scales from the Gram entries, solves for a positive null kernel in each component, and verifies a positive-definite principal minor. The weighted physical normal sums are positive multiples of the specified future null ray. The resulting component sizes are `9+5+8` and `13+8`; they cover the entire incident sets, are mutually orthogonal, and have total quotient rank 19.

These checks prove bounded links, rather than merely full rank. In the Euclidean horospherical quotient, each component's positive relation makes every recession pairing zero; the combined rank-19 span then forces the recession vector to vanish. A strict interior point exists from the global strict center sufficiently far towards the cusp. Every nonincident wall is strict at the null ray and is absent sufficiently deep in that horoball.

The new normals project onto the E7 face as `2ν1` and `2ν2`. These future null inequalities are strict throughout its hyperbolic interior and therefore add no face boundary. Finally, the independent checker reproduces all 87 negative pairings of the saved escape, its future orientation, and its positive square

`1336028658157091827 / 250000000000000000`.

Both the primary saved replay and the separate saved replay pass:

```sh
python3 scripts/fresh_h20_e7_intrinsic_caps_20261014_cusps.py --verify
python3 scripts/fresh_h20_e7_intrinsic_caps_20261014_cusps_review.py --verify
```

The old 85-wall all-pair certificate and ambient Lorentz signature remain explicit predecessor dependencies. This review independently checks all new arithmetic and the listed full incidences and vertices; it does not repeat the old pair census, enumerate every ambient vertex or cusp, or establish a general completion obstruction.
