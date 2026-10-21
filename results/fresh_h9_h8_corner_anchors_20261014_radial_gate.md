# The exact wall-8 anchored cap cannot support this H8 attachment

For the isolated cap direction in the [scalar-feasibility certificate](fresh_h9_h8_global_cap_gap_20261014.md), **every coherent perpendicular placement of the fixed compact 11-facet H8 source fails the necessary ball coverage condition for the raw attachment**. Here the chamber consists of the old twelve walls, this cap and the eleven perpendicular source sides. This holds for all rotations about the forced anchor. Additional nonperpendicular walls could cut the uncovered end; they, other cap directions and other sources are outside this theorem.

The [producer](../scripts/fresh_h9_h8_corner_anchors_20261014_radial_gate.py) and [certificate](fresh_h9_h8_corner_anchors_20261014_radial_gate.json) reconstruct the physical geometry using outward rational intervals with 192-bit square-root bounds. Generation and saved replay pass. The exact cap recipe and the source contact-gap classification are frozen dependencies; no numerical placement trial is used in the proof.

## Forced anchor and radial reach

Write `a=(sqrt(5)-1)/4`. The cap has normalized projected norm `kappa_8²=a`, while every other projected old normal has squared norm greater than `4a`. The [source gap theorem](fresh_h9_h8_facet_cover_20261014.md) therefore forces the normalized projected center of old wall 8 to be source vertex 0, 32, or 33, with contact scale one.

Let `c_i` and `rho_i` denote the center and radius of old ball i in the new H8 hyperplane. Any point covered by ball i has distance at most

`d(c_8,c_i)+rho_i`

from the anchor. The maximum of these twelve reaches is attained by ball 7 and is approximately **3.23739078077**.

The source has the following farther vertices:

| Forced anchor | Required source vertex | Exact distance cosh | Distance, approximate |
| --- | --- | --- | --- |
| 0 | 32 | `9+4sqrt(5)` | 3.57964082014 |
| 32 | 33 | `17+8sqrt(5)` | 4.24510024762 |
| 33 | 32 | `17+8sqrt(5)` | 4.24510024762 |

The checker verifies all 36 strict comparisons separately. It reconstructs these source distances from the literal source Gram and the three exceptional contact vectors, rather than taking the displayed decimals as evidence.

## Exact comparison and volume implication

Use the raw cap normal z with square q, old normals n_i with squares d_i, pairings `g_i=(z,n_i)`, and raw projected norm squares `k_i²=g_i²/q-d_i`. Wall 8 already has norm two and `k_8²=a`. For i different from 8,

`D_i = cosh d(c_8,c_i) = (g_8 g_i/q-(n_8,n_i))/sqrt(a k_i²)`.

The anchored reach has cosh

`D_i sqrt((k_i²+d_i)/k_i²) + sqrt(D_i²-1) sqrt(d_i/k_i²)`.

For i=8, the center distance is exactly zero. Every square root and denominator has a checked positive lower bound. The defining quadratic's opposite endpoint signs and nonzero derivative are rechecked on its frozen rational isolating interval.

In each possible placement, the specified source vertex lies strictly outside every closed old ball. Continuity leaves an open patch of the source interior uncovered. The [warped-end coverage lemma](fresh_h9_mixed_d4_global_caps_20261014.md) then gives an infinite-volume end. Thus this is a placement-independent failure of this particular cap/source attachment, rather than merely an optimization failure.

The separate [bounded rotation scout](fresh_h9_h8_corner_anchors_20261014_single_scout.md) also found no feasible relaxation in 24 starts. It is retained as numerical provenance and is not a certificate dependency. The scalar-feasibility report remains correct: its weaker unanchored diameter and projected-norm conditions do pass.
