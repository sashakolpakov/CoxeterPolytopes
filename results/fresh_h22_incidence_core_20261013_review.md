# Independent review of the H22 incidence core and binary orbit

**PASS within the stated construction scope.** The four incidence caps and the subsequent 560-wall orbit have exact compatible pair relations and a common strict timelike center. The first replacement retains a known strict spacelike escape. The 560-wall addition cuts that point; its pair certificate alone makes no statement about whole-chamber volume.

Reviewed sources are the [incidence core](../scripts/fresh_h22_incidence_core_20261013.py), [binary orbit](../scripts/fresh_h22_incidence_binary_20261013.py), and their [core](fresh_h22_incidence_core_20261013_full.json) and [orbit](fresh_h22_incidence_binary_20261013.json) certificates. Both full saved-result replays pass.

I also wrote a [separate standard-library integer/Fraction checker](../scripts/fresh_h22_incidence_core_20261013_review.py), which imports neither producer nor NumPy and has its own pair classifier. Its [saved audit](fresh_h22_incidence_core_20261013_review.json) checks the following directly from the authoritative old wall coordinates and frozen new normals:

- The 42 short walls split into two orthogonal 21-element parts, with incidence matrix satisfying `A Aᵀ=4I+J` and row weight five. The two prescribed integral projections have norms `-2`, mutual pairing `-5`, and sum equal to the old center.
- The four height-two caps have norm two and all 846 changed core pair relations: 84 right, 86 π/3, and 676 dotted. The center has strictly negative pairing with every wall.
- All 280 additional projections have norm `-2` and center pairing `-8`; both height signs therefore give 560 norm-two walls. All 117,600 pairs with the old 210 walls and all 2,240 pairs with the incidence caps are compatible.
- All 156,520 pairs within the 560 new walls are checked: 35,280 right, 40,320 π/3, 2,520 parallel, and 78,400 dotted.
- Every stated generator preserves the old form, center, and wall set; it permutes the 280 projections exactly. An explicit graph traversal reaches the entire set from one representative. Eight signed new caps cut the preceding escape.

## Transfer to the retained banks

The representative gate covers the whole 560-wall orbit because the retained banks are invariant under **the same** generator matrices. The earlier U and R certificates explicitly check generator closure of their projection arrays. The K certificate checks each cusp orbit separately and makes its assigned height constant on that orbit. Both added-coordinate signs are retained throughout.

For the repair banks, the archived tangent frames are equivariant. Every old finite component and its negative highest root are transported by the old-wall permutation. The lowered shift depends only on component size and norm, both invariant. Consequently the full shifted wall set is invariant even when a generator permutes component types within the archive. These facts justify transporting each representative new/old pair gate; no transitivity of a bank larger than its certified orbit is assumed.

The new gate checks integer products with overflow bounds, uses rational upper bounds for the signed radical terms, and handles the K-family right and π/4 contacts explicitly. Those contacts must not be described as all dotted. The lower certificate's source hash and its source manifests, including the exact frame and highest-root arrays consumed by the helper, are now explicitly checked. This hardening was added during review.

The resulting count is 5,397,550 after replacing the old two center walls with four incidence caps, and 5,398,110 after the binary orbit. Every retained wall is an actual facet: a common strict timelike center and all nonpositive distinct pairings let one project that center onto any wall while keeping every other inequality strict. Nonpositive pair checks also exclude positive-proportional duplicate normals.

This is an exact construction and pair-compatibility review, not a claim of a new finite-volume polytope, literature novelty, or completeness of arbitrary binary or continuous angle profiles. The full-bank representative calculations were replayed using their existing implementation; only the finite incidence/orbit algebra and pair calculations listed above were independently reimplemented.

```sh
python scripts/fresh_h22_incidence_core_20261013.py --verify --full-bank
python scripts/fresh_h22_incidence_binary_20261013.py --verify
python scripts/fresh_h22_incidence_core_20261013_review.py --verify
```
