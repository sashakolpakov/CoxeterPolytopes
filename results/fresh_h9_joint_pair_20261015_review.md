# Independent review of the 31-wall polar extension

The exact polar cap passes all 30 old-wall gates and creates one new finite vertex. The complete resulting chamber has **31 actual facets, 90 finite vertices, no ideal vertices, 858 hyperideal extreme rays, 251 finite edges and 308 open finite-vertex directions**. A new exact strict spacelike escape proves that its volume is still infinite.

The [independent checker](../scripts/fresh_h9_joint_pair_20261015_review.py) and [saved certificate](fresh_h9_joint_pair_20261015_review.json) bind the [polar certificate](fresh_h9_joint_pair_20261015_polar.json), the [31-wall primary certificate](fresh_h9_joint_pair_20261015_extension.json), and the [independently reviewed 30-wall starting cone](fresh_h9_nine_contact_extension_20261015_review.json). The arithmetic uses the existing independent SymPy algebraic-number-field implementation and 192-bit outward rational sign bounds. It imports no producer arithmetic or producer chamber engine.

The new normal is the literal old hyperideal ray with active set

    [0,1,2,3,5,6,7,10,29].

Its positive norm is retained as given; it is not rescaled to norm two. Its new pair relations are nine right angles and 21 strictly dotted pairs. The review also reconstructs the complete 31-wall Gram matrix and all 465 pair relations. Projection of the common strict timelike center onto each wall gives 31 explicit timelike facet points, with all other wall inequalities strict. Thus every displayed wall is an actual facet.

For the cone census, the independently certified completeness of the old 30-wall cone is the starting dependency. The checker reconstructs all 894 starting rays directly from their active-wall nullspaces and rechecks every old inequality and norm sign. The new wall retains 876 rays, removes 18, and contains none of them. A separate clipping routine finds 72 crossing edges by testing their complete common-face ray sets, rather than using the producer's adjacency test. It produces exactly 948 rays. Every resulting ray has its full active set, rank and all 31 inequalities rechecked. This verifies the complete new census by one exact cone cut; it is not a second construction of the old cone from scratch.

All 89 old finite vertices remain. The unique new finite vertex has incident walls

    [0,1,2,3,5,6,7,10,30].

It joins the old vertex with last incident wall 11 along the common rank-eight support `[0,1,2,3,5,6,7,10]`. Hence the cap closes an old open edge. The review checks all 90 positive-definite incident Grams, all 2,790 vertex inequalities, every finite edge, and the complete open-direction list. The incidence identity is `9×90 = 2×251 + 308`.

The old strict escape is cut only by the new wall. The replacement escape is reconstructed as a positive-center perturbation of the hyperideal ray with active set `[0,1,2,3,5,6,7,12,30]`. Its square is strictly positive, its pairing with the future center is negative, and all 31 wall inequalities are strict. These checks establish infinite volume without relying only on the presence of hyperideal rays.

Generation and saved replay pass. This review establishes the stated polar extension and its complete cone; it does not exclude other caps or a later finite-volume completion. An earlier bounded numerical endpoint scout missing this cap is a sampling limitation, not an all-parameter obstruction.

```sh
python scripts/fresh_h9_joint_pair_20261015_review.py --verify
```
