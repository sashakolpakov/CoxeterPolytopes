# Review of the orthogonal H9 pair

The [90-wall construction](fresh_h9_crossplane_repair_20261013.md) passes a separate source and algebra review and a 160-digit replay. It cuts both specified old gaps, but a new strict spacelike escape proves infinite volume.

The construction projects into the span of eight prescribed old normals and the exact vector `D=seed6−seed7`, with `D²=4`. Writing the first cap as `r=u+t v`, where `v` is perpendicular to that span, the chosen `t²=(2−u²)/v²` gives `r²=2`. The constraint `r.D=2` then gives the second cap `s=r−D`, with `s²=2` and `r.s=0`. The verifier anchors the physical seed6 and seed7 coordinates and checks the unchanged prescribed walls are orthogonal to D; the reflection exchanges the two remaining seed constraints.

Each new cap has exactly eight finite old neighbors. Any new positive-semidefinite diagram must lie in the finite-neighbor set of each new cap it contains. Enumerating the 256 subsets for each single cap and the 128 subsets for both therefore covers all 640 possibilities. The exact Q(phi) checks find no affine diagram and no new rank-nine positive-definite diagram. Thus these caps introduce neither an affine circuit nor a finite vertex.

Both target vectors are checked for positive square, future orientation and all 88 old strict inequalities, before the four new strict cuts are checked. All 90 inequalities hold strictly at the final exact escape. Strict center projection and nonpositive wall pairings establish that the walls are actual facets.

The [review script](../scripts/fresh_h9_crossplane_repair_20261013_pair_review.py) checks the norm/reflection identities symbolically and runs the existing interval engine at 160 digits. All 465 changed enclosures are nested in their archived 80-digit enclosures; exact metadata agree. This is a source review and precision replay, **not a second implementation of interval arithmetic**. The [machine-readable review](fresh_h9_crossplane_repair_20261013_pair_review.json) records the source hashes and scope.

```sh
python3 scripts/fresh_h9_crossplane_repair_20261013_pair_review.py --verify
```
