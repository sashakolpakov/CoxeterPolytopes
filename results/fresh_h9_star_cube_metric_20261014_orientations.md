# Sixteen metric branches of the star incidence

All sixteen choices of which leaf has own contacts `(3,3)` in each fixed B2 block are excluded for compact H9. The opposite leaf may have any admissible own contacts. This is a family of metric obstructions for case 0 of the tree catalogue; it is not an exclusion of that incidence or of the other 32 incidence representatives.

The proof retains the four internal B2 labels 4 and the literal positive support `[3,4,6,7,9,10,12,13]`. A free leaf and its `(3,3)` partner must be orthogonal: a nonright mutual contact would give a nonpositive all-ones quadratic on an actual three-wall finite face. The resulting actual paths exclude every own label at least 6, and the Lannér condition leaves eleven ordered profiles. Both retained-own-root choices are checked for every leaf, so a `(3,3)` leaf has seven possible projections and every other profile has only its own two contacts. The 24 cross-block right-contact proofs use the B2 labels alone, with all leaf contacts relaxed to right.

No leaf swap is assumed to be a polytope automorphism. For each of the sixteen literal choices, the checker recomputes the two finite-contact upper bounds from cut 14 to leaves 11 and 17, relaxing every unused free leaf to right. There are **2,112 exact actual-vertex upper-bound witnesses**. All restrict the two contacts to labels 2 or 3.

The eight orientations in which cut 14 is free have 8,092 residual cases each. The eight in which it is fixed have 3,332 each, for **91,392 occurrences** and **12,852 distinct cached parameter cases**. In each orientation precisely four residual determinants vanish; all four are excluded by exact nonpositive pivots inside actual finite vertices. Thus **91,328 occurrences fail rank at most ten**, and the remaining **64 fail finite-vertex positivity**. These are overlapping family counts, not disjoint incidence types.

The [certificate](fresh_h9_star_cube_metric_20261014_orientations.json) and [checker](../scripts/fresh_h9_star_cube_metric_20261014_orientations.py) use exact Q(√2,√5) arithmetic and literal incidence witnesses. Generation and saved-result replay pass. A [separate Fraction arithmetic review](fresh_h9_star_cube_metric_20261014_orientations_review.md) reconstructs all changed domain witnesses, all 12,852 cached residual cases and all 64 finite-vertex failures, without importing producer arithmetic. The preceding [complete cut-profile branch](fresh_h9_star_cube_metric_20261014_all_cut_profiles.md) is frozen separately and has an independent Fraction arithmetic replay.

The [bounded continuation](fresh_h9_star_cube_metric_20261014_orientation_search_case_00.json) reaches another mandatory-label assignment in 15 nodes. It changes the final positive leaf to `(3,4)` and its partner to `(5,2)`, so neither is `(3,3)`. A separate [exact metric certificate](fresh_h9_star_cube_metric_20261014_next_metric.json) excludes this particular state: all 30 domain witnesses are rebuilt, both anchor projections have no foreign contact, and all 28 third-leaf/contact cases have nonzero residual determinant. This additional result is a single-assignment exclusion. No further star search is claimed here.

A better next metric target is the saved T-shaped case 8. Its mandatory assignment uses three B2 blocks and one genuine I2(m) label category, m≥7. The potential positive support `[1,2,4,5,7,8,12,13]` could reduce its metric problem to two residual dimensions with the parameter cos(π/m), once its actual support and forced cross contacts are checked. It should be treated parametrically; category 7 is not the single angle π/7. This offers a different metric problem from further serial changes to the star leaves.

Replay:

```sh
python scripts/fresh_h9_star_cube_metric_20261014_orientations.py --verify
python scripts/fresh_h9_star_cube_metric_20261014_orientation_search.py --case 0 --limit 100000 --verify
python scripts/fresh_h9_star_cube_metric_20261014_next_metric.py --verify
```
