# Finite attachments beyond the right-S16 H20 wall family

The finite graph census is complete for a norm-two new wall meeting the actual H3 face with support
`T = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,61,87,151]`.
This support is `A10 + A4 + A2 + A1`. There are **57 fixed labelled contact profiles**, together with the single unbounded family that meets only root151 at angle `pi/m`, for any integer `m >= 3`. The subsequent wall search is numerical and bounded; it is not a completeness claim about all walls or chambers.

Generation and fresh saved replay pass for both the [profile/cache source](../scripts/fresh_h20_spherical_attachment_profiles_20261016.py) and [bounded scout](../scripts/fresh_h20_spherical_attachment_profiles_20261016_scout.py). Their outputs are [the profile/cache JSON](fresh_h20_spherical_attachment_profiles_20261016.json) and [the scout JSON](fresh_h20_spherical_attachment_profiles_20261016_scout.json). Source/data/report are final for this package. No predecessor was edited.

## Why the finite profiles are complete

The connected finite component containing the new root contains whole old irreducible components: deleting the new vertex cannot separate vertices that were already joined by old edges. Its deletion is therefore exactly a subset of the four labelled old paths. Its rank is at most18. The source enumerates every marked deletion and every edge-labelled isomorphism in the finite irreducible types A, B, D, E, F, H, through rank18, including the isolated new A1. A rank-two component can use only the isolated old A1, root151; the entire I2(m) family is recorded symbolically instead of being truncated.

All roots in this graph calculation have square2. Thus a label m contributes magnitude `2 cos(pi/m)`, including the B/F and H cases; no root-length gauge is inferred merely from graph names. The exact rational support inverse gives the projection cost `d^T J d`. Every saved fixed profile has cost strictly below2, certified symbolically. For the parametric family the residual square is `2 sin(pi/m)^2 > 0`.

There are four E8 profiles, joining either end of A4, either end of A2, and root151, all at label3. Their cost is `4/5 + 2/3 + 1/2 = 59/30`, so the projected normal in the H3 Lorentz space has square `1/30`. These profiles were absent from every earlier bank constrained to be right to S16.

## Native ambient gates

The scout does not apply only the residual Lorentz Gram. For support projection `p = -T J d` and residual root y, its pairing with every old root is `y.N_projected + p.N`. The cache reconstructs all original107 physical covectors in the saved radical tower and the later physical coordinate enclosures through wall157. It verifies the full158-by-158 Gram decomposition at75 decimal digits before conversion to scout doubles. The maximum relative reconstruction error is approximately `2.42e-62`.

Native positive root squares are saved separately. Direct double-precision subtraction of large nearly null coordinates is unsafe for some old roots and is not used to infer their squares. The cache remains a numerical aid, not an independent exact geometry certificate.

## Bounded search outcome

The three-additional-old-neighbor scout uses the actual old H3 facets `[95,142,145,150,152,156,157]` and prescribed labels in `[2,3,4,5,6,8,10]`. Its 17,075 contact patterns yield551 numerical roots satisfying all158 native gates, of which226 are genuinely nonright to S16. No single tested wall closes the H3 face. The E8 profiles have no successful three-old-neighbor root in this bounded domain, but36 E8 rays with every other old contact dotted are found using nine old-exterior directions and all four placements.

The old face has16 rays in the scout's fixed future-center chart, with maximum squared radius `3.06383257946`. Some E8 rays lower this maximum to `2.36107631150`. These values must not be compared to radii in earlier, different affine charts. Sphere classifications use a numerical band of `1e-6`.

The saved comparison initially exposed integer dictionary keys becoming strings in JSON. This serialization issue was corrected before final generation/replay; the numerical census was unchanged. The fixed label lists, chosen ray directions, tolerance tests, and finite sampling are explicit bounds on this experiment, not all-parameter exclusions. A coordinated E8 triangle experiment is recorded under a separate prefix and is not part of this frozen package.
