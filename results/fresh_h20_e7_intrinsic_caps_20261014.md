# A complete H13 face after four-wall replacement

**2026-10-14 — exact generation PASS.** Removing four walls of the frozen [87-wall chamber](fresh_h20_e8_interior_caps_20261014.md) and adding two transferred source facets produces an **85-facet H20 chamber with a complete finite-volume H13 face**. The face has 15 facets, 44 finite vertices and two ideal vertices. Both ideal vertices have unbounded ambient links, and the previous strict spacelike escape survives. Thus the ambient chamber still has infinite volume.

The [exact checker](../scripts/fresh_h20_e7_intrinsic_caps_20261014.py) and [certificate](fresh_h20_e7_intrinsic_caps_20261014.json) preserve all previous artifacts. The construction uses the complete [A11+A2 H13 source](fresh_h13_a11a2_completion_20261014.md) extracted from an E8 face of the Borcherds chamber. It independently repeats the small rational projective census after transferring the source into the present E7 face.

## Transfer and the necessary replacement

Use the [normalized E7 face basis](fresh_h20_e7_face_20261014.md): A11 simple roots at old indices `0,…,10`, A2 simple roots at `18,83`, and a future timelike vector `T` of square −2. Write `λ1,…,λ11` and `μ1,μ2` for the two systems of fundamental weights. The chosen transferred source facets are

```text
R = −λ11 − √2 μ1 + (√2/4) T,
S = −λ3             + (√2/4) T.
```

Both have norm two and are mutually orthogonal. The A2 order is swapped relative to the source certificate. The checker compares all 15 source normals, coordinate by coordinate, under this isometry.

Keeping all 87 old walls prevents this transplant. There are four identifications of the specified source corner with the present corner: A11 reversal and A2 interchange. For all four, the right lifts obey

```text
(R, old26) = √2 − √22/4 > 0,
(S, old26) = −√22/4.
```

The second pairing has normalized squared cosine `11/32`, strictly between `1/4` and `1/2`; hence it is not a Coxeter angle. The only alternative lift permitted by the E7 meeting-profile theorem is `y−2ω33`, of norm eight. Each of these alternative lifts retains a pairing −1 with a norm-two A11 root, giving normalized cosine `1/4`, strictly between right and π/3. Thus all 16 specified individual lift tests fail. This excludes this fixed-corner transplant while retaining all 87 facets, not every possible H13 source embedding or ambient completion.

The coherent replacement removes exactly old indices

```text
[17, 26, 75, 76]
```

and adds the displayed right lifts `R,S`. The certificate saves the old-to-new index map; the new facets have indices 83 and 84. All 167 new pair relations are exact: 38 right, two at π/3, one at π/4 and 126 strictly dotted. The original strict center remains strict. Its orthogonal projection onto each of the 85 walls is future timelike and strict against every other wall, proving that all 85 inequalities are actual facets.

## Complete face geometry

The E7 support is retained. Its face has exactly the thirteen original A11+A2 walls and the two new source walls. Every other retained non-E7 wall projects to a future timelike vector, so it contributes no boundary anywhere in the hyperbolic face. These 63 orientations and norms are checked explicitly.

For the face computation, scale the two A2 basis roots and `T` by √2. The rational form becomes

`diag(C_A11, 2 C_A2, −4)`.

All 15 face normals then have rational coordinates. An explicit positive combination of all 15 normals equals half the timelike basis vector. This gives a bounded projective section; enumerating all `15 choose 13 = 105` active sets therefore exhausts its extreme rays. The independent transferred census finds exactly 44 negative-norm rays and two null rays, with no spacelike extreme ray. All 44 finite points are checked against every one of the 85 ambient inequalities. Their incidence is exactly E7 plus a positive-definite rank-13 face support, so they are actual rank-20 ambient finite vertices.

Each intrinsic ideal link has positive affine kernels covering every incident face wall and quotient rank 12. The checker verifies a positive multiple of the corresponding null ray for each weighted root sum, together with positive-definite principal minors. Thus both intrinsic links are bounded, and the H13 face has finite volume. This is a complete face census; no complete ambient vertex or cusp census is asserted.

## Remaining ambient obstruction

At each of the two face cusps, the complete ambient incident set consists of the intrinsic affine components together with the spherical E7 block. The vector `−ω33` has square `3/2`, is orthogonal to the null ray, and has nonpositive pairings with every incident wall, with at least one strict pairing. It gives a nonzero Euclidean link recession. Therefore the two ambient cusp links remain unbounded even though their H13 links are bounded.

The old global escape also remains strict against all 85 walls. Its squared norm is

`1336028658157091827 / 250000000000000000 > 0`.

The result supplies a complete finite-volume face after a concrete change of retained walls. It does not complete the H20 chamber, and no nonexistence claim for further completion follows from this checkpoint.

```sh
python3 scripts/fresh_h20_e7_intrinsic_caps_20261014.py --verify
```
