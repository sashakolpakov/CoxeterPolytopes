# Three escape-cutting caps without an affine circuit

This alternative retains the independently checked49-wall chamber and adds three exact compatible caps. All three cut its previous strict escape. Their complete52-wall cone has **178 finite vertices, no ideal extreme rays, and2,279 exterior rays**, with a new exact strict spacelike escape. It therefore still has infinite volume. Unlike the separately archived [first three-cap attempt](fresh_nonarithmetic_escape49_20261016.md), this52-wall chamber contains no affine circuit at all.

The [exact cap data](fresh_nonarithmetic_affine_filtered_20261016_cluster.json) and [checker](../scripts/fresh_nonarithmetic_affine_filtered_20261016_cluster.py) prescribe the following finite contacts. Every omitted label on `S+{33}`, where `S=[0,1,2,3,5,6,7,30]`, is2.

| New wall | Nonright contacts |
|---|---|
|49|old30:5, old7:3, old33:10|
|50|old3:3, old7:5, old33:10|
|51|old5:6, old6:4, old33:10|

Every other old contact is strictly dotted. All three new mutual pairs are strictly dotted. Exact square-root recipes give norm2 and the nine prescribed contacts; outward rational intervals certify the40 remaining old contacts for each cap, the three mutual products, strict center signs and positive pairings with the previous escape. The common center and acuteness also prove that all52 walls are actual facets. The cap formulas use `2cos(pi/10)=sqrt((5+sqrt(5))/2)` and the positive norm-equation root, with no approximate equation solving in the acceptance certificate.

The [complete cone certificate](fresh_nonarithmetic_affine_filtered_20261016_cone.json) and [checker](../scripts/fresh_nonarithmetic_affine_filtered_20261016_cone.py) reproduce the frozen36→48→49 cone before the three new clips. The clips create25,56 and55 crossing rays. All178 finite vertices of the49-wall chamber survive, and none is on a new cap. The caps cut5,21 and10 old extreme rays respectively, with22 distinct old directions removed. No new finite vertex is created. The renewed escape is the surviving ray with active set `[0,1,2,3,5,6,8,34,44]`, perturbed by `H/262144`; its norm and all52 strict signs are certified. The complete clipping is a producer implementation with saved replay, not an independently reconstructed52-wall census.

## The affine filter is complete for this chamber

The [exact graph gate](fresh_nonarithmetic_affine_filtered_20261016_affine_gate.json), produced by [this checker](../scripts/fresh_nonarithmetic_affine_filtered_20261016_affine_gate.py), binds the separate [complete old49 affine inventory](fresh_h9_persistent_cusp_affine_20261016.json). That inventory contains no irreducible affine circuit anywhere in the retained chamber, not merely none orthogonal to an old forbidden-angle edge.

A new irreducible affine circuit contains at most one new cap, since every pair of new caps is dotted. All its old nodes must belong to that cap's nine finite neighbors. Apart from the rank-two parallel case, irreducible affine Coxeter diagrams have only labels2,3,4,6. There are no parallel pairs here. The gate enumerates every connected subset containing a new cap and using those labels; all14 matrices have a nonzero exact symbolic determinant. This rules out every new affine circuit. It also covers a new forbidden-angle ridge orthogonal to an old affine circuit, because the old inventory is empty.

This check was necessary: testing only affine triangles had missed an affineF4 path in an earlier numerical candidate. Absence of affine circuits does not establish finite volume or guarantee that a completion exists. It only avoids the permanent forbidden-angle cusp obstruction that invalidated the first52-wall attempt.

## A finite-volume continuation must meet the new caps

Each new cap has an actual pi/10 ridge with old33, but its present facet contains no finite vertex. In a finite-volume Coxeter polytope, an ideal vertex cannot lie on a pi/10 ridge: the bounded Euclidean cusp link would contain two reflections whose product has irrational linear trace. Hence that ridge must have a finite vertex. This is the previously reviewed [forbidden ridge-label argument](fresh_h20_weighted_binary_20261014_dotted_completion_review.md).

Any additional wall disjoint from a new cap has a future nonspacelike projection onto that cap's supporting hyperplane, so its inequality is strict on the entire finite facet. This also holds for a parallel wall. Therefore additions that all miss a particular new cap cannot supply the finite vertices required on its pi/10 ridge. **A finite-volume continuation retaining this chamber must contain an added wall meeting each of the three new caps.** One added wall could meet more than one; no lower bound of three added walls is asserted.

The [bounded bridge diagnostic](fresh_nonarithmetic_affine_filtered_20261016_bank_bridge.json) finds no numerical finite-angle partner for any selected cap among the earlier1,846 proposals, and no additional proposal numerically compatible with all three. The respective strict-dotted partner counts are38,1,295 and317. This is a numerical diagnostic only, with near-parallel cases kept unresolved; it does not exclude new directions, contacts or fields. It explains why extending this chamber requires a new contact family rather than another selection from the same saved bank.

All four scripts (`cluster`, `cone`, `affine_gate`, `bank_bridge`) passed generation and saved `--verify` replay, source hashes, syntax and report links. The exact local construction, full cone and affine gate have distinct scopes. No finite-volume example or new dimension is claimed.
