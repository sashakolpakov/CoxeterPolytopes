# Exact frontier of the 27-wall joint chamber

The frozen [27-wall H9 chamber](fresh_h9_makarov_joint_20261013.md) has **exactly 65 finite vertices**. Its finite-vertex graph is connected and has **144 finite edges**. Another **297 edges reach the ideal boundary without a second existing finite or ideal vertex**. Each of these directions has an exact feasible spacelike continuation, so the present chamber remains infinite-volume.

This replaces the numerical frontier proposals with an exact census of all finite corners and every edge incident to one. It does not enumerate exterior projective vertices or edges having no finite endpoint. The numerical scout's 776 exterior incidence ambiguities are therefore outside the claimed census.

## Complete finite-corner and edge census

Every finite vertex of a Coxeter chamber is simple. Its nine incident normals have positive-definite Gram and every pair meets at a finite angle. The exact finite-pair graph of the 27 walls has only 351 nine-cliques. Exhausting these avoids the much larger unrestricted `27 choose 9` search:

| Exact result | Count |
| --- | ---: |
| Finite-pair nine-cliques | 351 |
| Not positive definite | 286 |
| Positive definite, feasible finite corners | 65 |
| Previously certified corners | 50 |
| Additional finite corners | 15 |
| Distinct edges between finite corners | 144 |
| Open edges with one finite end | 297 |
| Feasible ideal extensions of these edge supports | 0 |

For each positive-definite support, projecting the strict center away from its nine normal directions constructs its exact timelike vertex. All 27 wall signs are then checked, including exactly nine zero pairings. The 15 additional corners involve the cap pairs `(13,14)`, `(17,18)`, `(21,22)`, `(23,24)`, and `(25,26)`, contributing respectively 3, 4, 3, 3, and 2 corners. They are not merely intersections proposed by numerical tolerances.

Releasing each of the nine walls at each corner gives 585 directed incidences and 441 distinct rank-eight supports. For every support, all remaining existing walls are tested. A second finite endpoint requires a positive Schur complement; an ideal endpoint requires a zero complement and feasible null vector. Dotted pairs cannot occur in either case. There are no affine extension candidates in this census. The resulting finite-endpoint incidences exactly agree with the complete corner list.

For each of the 297 open supports, the archive additionally gives a future-oriented spacelike point satisfying all inequalities, with precisely those eight equalities. Its positive norm and all signs are exact. The segment from its finite corner to this point crosses the light cone while remaining feasible, certifying an open end independently of an approximate collision ordering. Every such support contains current caps: 214 contain one cap and 83 contain two. These are new directions after closure of the original forty seed edges.

## Profiles for constructive continuation

The [certificate](fresh_h9_joint_frontier_20261013.json) records **10,057 finite spherical endpoint profiles**, together with **418 isolated dihedral families** admitting every integer label `m >= 3`, for the 297 open supports. These are per-edge profile records, not distinct new mirrors.

All walls are normalized to norm two, including polar wall 6. A positive-definite Coxeter diagram has no cycles and no node of valence four. A node of valence three has three label-three edges. A label at least six can occur only in an isolated two-node component: adding any further incident edge already gives a nonpositive three-node Schur complement. Thus the new node can attach to at most one node of each old component, using labels 3, 4, 5 except for the isolated dihedral families. The checker exhausts these possibilities and keeps precisely the positive Schur complements. All-right profiles are included.

This is a complete **local spherical** profile census. It supplies necessary contacts for a finite endpoint; it does not prove that the corresponding cap exists, is globally compatible, or closes volume. New affine endpoint profiles are not included.

A low-field target sent to the joint cap search is **edge 235**. It starts at additional corner 34, releases polar wall 6, and retains

```text
[0, 1, 3, 4, 7, 10, 13, 14].
```

Both caps 13 and 14 already lie in `Q(sqrt(5))`, as does this entire support. A candidate with label-three contacts to walls 3 and 13, and right contacts to the other six support walls, has residual norm `2/3`. Exchanging cap 13 for 14 gives another allowed profile. This is an exact new open-edge target rather than another cut of the old forty directions.

Other examples with two current caps are edge 182, with support `[1,3,4,5,8,9,17,18]`, label-three contacts to 1 and 17, and residual norm `1/2`; and edge 354, with support `[0,2,3,4,6,8,25,26]`, label-three contacts to 0 and 25, and residual norm `1`.

## Exact arithmetic and independent replay

The [producer](../scripts/fresh_h9_joint_frontier_20261013.py) uses the prior separate `Fraction` arithmetic for `Q(sqrt(5))`, with sparse products of positive radicals for actual coordinates. Positive rescaling of the polar normal puts every finite Gram entry in the base field. Radical generators need not be algebraically independent: exact recursive signs remain valid with dependencies. Floating-point values only propose the rational step used for a spacelike continuation; every accepted identity and inequality is checked exactly.

The [separate checker](../scripts/fresh_h9_joint_frontier_20261013_independent.py) imports no producer arithmetic. It reconstructs every normalized mirror from the frozen physical normal, verifies coefficient identities, and proves every strict sign with outward rational square-root intervals at 160-bit precision. It independently exhausts the 351 nine-cliques, rechecks all 65 points and 441 edges, re-enumerates all spherical profiles, and verifies the source chamber's strict spacelike escape. Its [audit certificate](fresh_h9_joint_frontier_20261013_independent.json) binds the sources by SHA-256.

```sh
python3 scripts/fresh_h9_joint_frontier_20261013.py --verify
python3 scripts/fresh_h9_joint_frontier_20261013_independent.py --verify
```

Both exact implementations pass. A further [source/proof review and fresh replay](fresh_h9_joint_frontier_20261013_review.md) found no issue in the corner, edge, or profile completeness arguments. The frontier and profile data have been handed to the ongoing constructive cap search; no new finite-volume or compact polytope is asserted.
