# Symmetric 862-facet application of the one-side reduction

The [exact checker](../scripts/fresh_h22_one_side_20261013.py) reflects the 25 cap types of the preceding 839-facet π/7 construction across the old H21 hyperplane. Its 50 signed caps together with the 812 rational parents give **862 actual facets**, containing all 839 preceding facets. The [certificate](fresh_h22_one_side_20261013.json) and [source/proof/replay review](fresh_h22_one_side_20261013_review.md) pass.

All 1,225 signed-cap pairs are verified: 1,217 are strictly dotted, six are right and two are parallel. Height reflection exactly permutes the parent normals and preserves the strict center. The two previously bounded cusp links and the actual compact quadrilateral survive; reflection supplies a second compact quadrilateral.

The [general one-side theorem](fresh_h22_extension_feasibility_20261013.md) shows that this is a legitimate search reduction for extensions retaining all old vertical facets: any coherent upper cap family can be reflected without introducing forbidden cross angles. A finite-volume extension exists within that scope only if one side can cover the entire old chamber. The reduction includes the exceptional pure-vertical cap case.

The integer spacelike witness of squared norm 1,586 remains strict on all 862 walls, and a specified inherited cusp remains untouched. **This chamber has infinite volume.** The symmetry improvement does not establish a new finite-volume example. The [separate π/6 replacement](fresh_h22_multiquadratic_20261013.md) changes the angle field and is a different construction, not an extension retaining these π/7 normals.

Replay:

```sh
python3 scripts/fresh_h22_one_side_20261013.py --verify
```
