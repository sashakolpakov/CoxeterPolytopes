# Independent review of the two-cap 32-wall chamber

The two exact caps have square two, meet at angle π/3, and pass every old-wall Coxeter gate. Both cut the old strict escape. Their common finite vertex is genuine and requires both caps. The complete alternative chamber has **32 actual facets, 90 finite vertices, no ideal vertices, 882 hyperideal extreme rays, 250 finite edges and 310 open finite-vertex directions**. A new exact strict spacelike escape proves infinite volume.

The [independent checker](../scripts/fresh_h9_joint_pair_20261015_two_review.py) and [certificate](fresh_h9_joint_pair_20261015_two_review.json) bind the [two-cap formulas](fresh_h9_joint_pair_20261015_twocaps.json), the [primary 32-wall geometry](fresh_h9_joint_pair_20261015_two_extension.json), and the [independently certified 30-wall starting cone](fresh_h9_nine_contact_extension_20261015_review.json). The checker uses the separate SymPy quartic-field implementation and 192-bit outward rational sign bounds. It imports neither producer arithmetic nor the producer's cone engine.

Both formulas use contact set `[0,1,2,3,5,6,7,10,29]`. The first cap has label 5 at wall 7 and label 3 at wall 29; the second has label 5 at walls 3 and 7 and label 3 at wall 5. The other prescribed contacts are right. The independent calculation reconstructs the common kernel, checks the particular solutions and the two literal positive square roots, verifies the square identities, and checks both solution branches. In each case exactly one branch satisfies the old Coxeter gates, the inherited strict-center inequality and the escape-cutting inequality. All coefficients stay in the existing quartic field.

All 496 pair relations are checked directly. The 61 new relations comprise 13 right, six finite (three π/3 and three π/5), and 42 strictly dotted pairs. Projection of the common strict timelike center onto each wall gives 32 actual facet witnesses, with all other inequalities strict. Both new normals have strictly positive pairing with the frozen polar normal, so this is an alternative to the [31-wall polar branch](fresh_h9_joint_pair_20261015_review.md), not an extension of it.

The complete cone review starts from the previously independently certified 894 old rays. Every starting ray is reconstructed from its exact active-wall nullspace and all old signs are rechecked. The first cut retains 880 rays and removes 14, producing 64 crossing rays and 944 rays in total. The second cut has 900 negative, 20 positive and 24 zero old rays; 48 new crossings give 972 final rays. Crossing adjacency is determined by the complete ray set of each common face, rather than the producer's adjacency criterion. Every final ray has its full active set, rank, norm and all 32 wall inequalities rechecked. Thus the new census is complete conditional on the already certified starting cone; the old cone is not regenerated from scratch.

All 89 old finite vertices remain. The unique new finite vertex has incident walls

    [0,1,2,5,6,10,29,30,31].

Its positive-definite diagram is H3+D6. The H3 edges are `(0,1):5` and `(1,2):3`; the D6 edges `(5,6)`, `(5,10)`, `(5,31)`, `(29,30)` and `(30,31)` all have label 3. The vertex has no finite incident edge and all nine of its directions are open. The review checks all finite-vertex Grams and 2,880 vertex inequalities, the complete finite-edge list and the open-direction list. The incidence identity is `9×90 = 2×250 + 310`.

The new escape is reconstructed as a positive-center perturbation of the hyperideal ray with active set `[0,1,2,3,5,6,7,10,30]`. It has strictly positive square, negative pairing with the future center, and strictly negative pairing with every wall. This proves infinite volume for the displayed chamber. The local construction does not claim a finite-volume completion or exclude further caps.

Generation and saved replay pass; source bindings, syntax and report links are checked.

```sh
python scripts/fresh_h9_joint_pair_20261015_two_review.py --verify
```
