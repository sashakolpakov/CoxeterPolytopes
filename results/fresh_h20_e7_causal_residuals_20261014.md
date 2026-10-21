# Mixed E7 residuals and an exact two-anchor cap

**2026-10-14 — exact construction PASS; H20 volume remains infinite.** A new norm-two cap retains all walls of the [completed 87-wall E7-face branch](fresh_h20_e7_intrinsic_caps_20261014_cusps.md), makes actual π/3 ridges with old walls 61 and 75, and cuts both previously saved global escapes. Those two old facets are dotted to one another. The resulting 88-wall chamber preserves its complete finite-volume H13 face, all 44 finite face vertices and both completed ambient cusp links. A new rational strict spacelike escape survives. The new facet has no finite vertices.

The final [checker](../scripts/fresh_h20_e7_causal_residuals_20261014_bridge.py) and [certificate](fresh_h20_e7_causal_residuals_20261014_bridge.json) are separate from every frozen predecessor. The [independent bridge review](fresh_h20_e7_causal_residuals_20261014_bridge_review.md) records its own verification scope. This package concerns new finite contacts and changed core choices; it does not provide a finite-volume H20 example.

## Causal reduction from the complete H13 face

All indices below refer to the completed 87-wall input. Its E7 roots are `[12,13,14,15,18,31,65]`; write their fundamental weights as `ω0,…,ω6`. The thirteen corner walls of the H13 source are `[0,1,2,3,4,5,6,7,8,9,10,17,79]`, with fundamental weights `λ0,…,λ12`. The future timelike vector `T` is orthogonal to this entire rank-20 corner and has square `−2`.

Every additional acute normal has a decomposition `r=−w+y`, where `w` is a nonnegative combination of the E7 fundamental weights and `y` lies in the H13 face space. Its pairings with the fifteen retained face normals are nonpositive. The [complete source certificate](fresh_h13_a11a2_completion_20261014.md) and [physical transfer](fresh_h20_e7_intrinsic_caps_20261014.md) prove that this linear cone has precisely 44 timelike and two null extreme rays. A positive combination of all fifteen source normals gives a future timelike vector, so its nonzero vectors admit the same bounded projective normalization used by the source census. Consequently `y` is future causal, or zero: the entire cone lies in the convex future light cone. This conclusion uses the full cone certificate, not just the signs at a sample of vertices.

For a norm-two cap, therefore, `y²=2−w²≤0`. This permits a search in intrinsic causal directions while allowing nonzero, mixed E7 contacts. It does not restrict those directions or their coefficients to a rational lattice.

## Final 88-wall bridge

The cap is

```text
r = −√3(ω4+ω6) − a λ0 − b λ12 + t T,
t = (22√3+1)/(2√11).
```

The positive coefficients `a≈2.0606667215` and `b≈2.1623151140` are specified exactly by

```text
(11/12)a² + (2/3)b² − 2t² = −125/2,
(11/12)a + ((99−√869)/24)b
  + ((√79−33√11)/8)t + 75√3/2 = −1.
```

The source eliminates `a` linearly, records the resulting quadratic and its positive discriminant, and chooses its minus branch. It checks `a>2`, `b>2`, all norm and angle identities, and both escape cuts exactly. Coefficients are represented in `K(√3,√Δ)`, where `K=Q(√2,√79,√11)` and the explicit positive `Δ` is saved. No assertion about the minimal degree of this representation is needed.

The E7 contact magnitudes are `√3` at positions 4 and 6, giving π/6 with current walls 18 and 65, and zero at the other five E7 walls. The two anchor equations give pairing `−1` with each of walls 61 and 75. The 87 new relations are:

| Relation | Count |
| --- | ---: |
| Right | 16 |
| π/3 | 2 |
| π/6 | 2 |
| Strictly dotted | 67 |

The old strict center remains strict. For every wall, projection of that center onto its hyperplane is future timelike. All other inequalities remain strict by acuteness: in `p=H−(H·r/r²)r`, their pairings can only decrease. This proves all 88 facets actual. The checker also gives explicit timelike relative-interior points of the two new π/3 ridges and checks them against all 88 inequalities.

The new H13 projection is future timelike with square `−125/2`, so the additional inequality is strict throughout the hyperbolic face. The checker verifies strictness at all 44 finite face vertices and at both null cusp rays. Thus the existing positive affine link components and quotient rank 19 at those two cusps are unchanged. This is preservation of the specified complete H13 face and two bounded ambient links, not a census of every ambient cusp.

There are exactly twenty old finite neighbors of the new cap. Every potential finite vertex on it contains the cap and nineteen of those neighbors. All twenty such supports are checked: an old nonfinite pair, a non-positive-definite old Gram, or a negative new Schur complement excludes each. Hence the new cap creates **no finite vertex**. Its two actual finite ridges connect retained facets, but do not establish edge closure or volume closure.

The cap cuts the original input escape and the later [common escape for 42 proposals](fresh_h20_e7_replacement_escape_20261014.md). A [fresh numerical proposal](fresh_h20_e7_two_anchor_escape_20261014_scout.json), rounded to rational ambient coordinates and then checked exactly in the bridge certificate, remains strict against all 88 walls. Its square is

`259854917001622621 / 62500000000000000 > 0`.

The chamber is therefore still infinite-volume. The four actual old facets in the [conditional non-quasi-arithmetic witness](fresh_h20_e7_intrinsic_caps_20261014_cusps_nonquasi.md) all remain; any finite-volume completion retaining those facets would be non-quasi-arithmetic. No arithmetic finite-volume H20 search is being proposed.

## Two smaller changed-core checkpoints

The [86-wall checker](../scripts/fresh_h20_e7_causal_residuals_20261014_exact.py) and [data](fresh_h20_e7_causal_residuals_20261014_exact.json) delete only old walls 61 and 75, then add

`−ω2−ω4−ω6−√2 λ0−2λ12+5T`.

It cuts the original escape and preserves the face and both filled cusps. All 85 new pairs and all actual facets are exact. The later common escape survives.

The [one-deletion refinement](../scripts/fresh_h20_e7_causal_residuals_20261014_contact.py) and [data](fresh_h20_e7_causal_residuals_20261014_contact.json) instead retain wall 75 and use

`−ω2−ω4−ω6−√2 λ0−√(86/11)λ12+(17/√11)T`.

It deletes only wall 61, gives an actual π/3 ridge with retained wall 75, and again preserves the face and specified bounded links. Both constructions have no finite vertex on their new facet and have the same exact surviving escape of square `4374521734010421633/10^18`. Their [independent review](fresh_h20_e7_causal_residuals_20261014_review.md) passes.

The separate [continuous certificate](fresh_h20_e7_replacement_escape_20261014_continuous.json) proves that escape survives every real parameter in

`−ω2−ω4−ω6−√2 λ0−xλ12+sqrt((71+x²)/3)T`.

The final bridge changes both the E7 contact profile and the free source coefficients, so it is outside this closed curve.

## Bounded numerical search scope

These numerical scouts propose formulas; none supplies an exhaustive all-label exclusion:

| Scout | Scope | Result |
| --- | --- | --- |
| [Mixed edge events](fresh_h20_e7_causal_residuals_20261014_scout.json) | 315 source edges; 56 two-/three-node π/3 E7 profiles; one anchor label 2,3,4,6 | 195 acute cutters; no accepted numerical Coxeter proposal |
| [Wall-65 edge events](fresh_h20_e7_causal_residuals_20261014_single_scout.json) | All 315 edges; wall65 labels 3,4,5,6; same anchor pool | 340 acute cutters; no accepted proposal |
| [Interior contact profiles](fresh_h20_e7_causal_residuals_20261014_interior_scout.json) | 64,662 sparse source profiles and 66 E7 profiles | 1,128 acute cutters; 42 pass after deleting61/75; none mutually compatible numerically |
| [Two-anchor events](fresh_h20_e7_causal_residuals_20261014_two_anchor_scout.json) | 245 E7 profiles; residual in `span(λ0,λ12,T)`; labels2,3,4,6 at61/75 | 29 numerical caps retaining87; no mutually compatible pair; selected14 exactified above |

The numerical finite-angle screening used labels through 64. The exact construction checks only its actual identities and strict dotted inequalities, without relying on a cosine fit. Root's common-42 theorem is exact; the common-29 scout remains numerical except for the selected cap's independently checked rational escape. The next constructive requirement is a new finite endpoint or a different residual support, rather than additional caps from these mutually incompatible batches.

Generation and saved replay of the exact constructions:

```sh
python3 scripts/fresh_h20_e7_causal_residuals_20261014_exact.py --verify
python3 scripts/fresh_h20_e7_causal_residuals_20261014_contact.py --verify
python3 scripts/fresh_h20_e7_causal_residuals_20261014_bridge.py --verify
```
