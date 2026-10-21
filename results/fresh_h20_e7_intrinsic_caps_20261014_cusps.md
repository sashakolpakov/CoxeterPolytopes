# Completing both ambient E7 cusp links

**2026-10-14 — exact generation PASS.** Two additional facets complete both unbounded ambient cusp links of the [85-wall replacement](fresh_h20_e7_intrinsic_caps_20261014.md). The resulting **87-wall branch** retains its complete finite-volume H13 face and all 44 finite face vertices. A strict global spacelike escape still survives, so the H20 chamber remains infinite-volume.

The [checker](../scripts/fresh_h20_e7_intrinsic_caps_20261014_cusps.py) and [certificate](fresh_h20_e7_intrinsic_caps_20261014_cusps.json) are separate from the frozen 85-wall construction and the earlier, different 87-wall chamber.

Let `ν1,ν2` be the two exact null rays saved in the replacement certificate, normalized to last intrinsic coordinate one. Let `ω15` be the E7 fundamental weight at old wall 15, whose inverse-Gram diagonal is two. The two new normals are

```text
r1 = −ω15 + 2ν1,
r2 = −ω15 + 2ν2.
```

Both have norm two. Each makes angle π/3 with wall 15, is right to the other six E7 walls and all incident intrinsic cusp walls, and is strictly dotted to every nonincident retained wall. They are also mutually strictly dotted. The checker starts with common shift one and doubles once; shift two passes every exact pair and center gate. The 171 new pair relations comprise 39 right pairs, two π/3 pairs and 130 strictly dotted pairs.

The original strict center remains strict. Projection of that center onto each of the 87 walls gives a future timelike point strict against every other inequality, proving all 87 facets actual.

For each cusp, the checker determines the full incident set against all 87 normals. It preserves the intrinsic affine components and replaces the spherical E7 factor by affine E7. The new positive null relation is

`Σ (B_E7⁻¹)3j e_j + r_i = 2ν_i`.

All coefficients are strictly positive. Together with the inherited intrinsic relations, the affine block sizes are `9+5+8` at the first cusp and `13+8` at the second. Every block has a positive-definite principal minor and a positive kernel whose weighted normal sum is a positive multiple of the same null ray. The blocks cover the entire incident set, are mutually orthogonal and have total quotient rank 19. These facts prove both Euclidean links bounded; rank alone is not being used as a boundedness criterion.

Both new walls are strict on all 44 finite face vertices. Their projections onto the E7 face are the future null vectors `2ν_i`, so they add no hyperbolic face boundary. The complete 15-facet finite-volume H13 face is unchanged.

The saved global escape remains strict against all 87 walls and has square

`1336028658157091827 / 250000000000000000 > 0`.

Thus this step repairs the two specified cusps, while leaving a separate infinite-volume region. It makes no complete ambient vertex/cusp census or general completion claim.

```sh
python3 scripts/fresh_h20_e7_intrinsic_caps_20261014_cusps.py --verify
```
