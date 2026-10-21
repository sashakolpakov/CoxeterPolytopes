# Replacing the label9 wall59 bond: the old prism is the blocker

Releasing old walls59,60,61 admits89 numerical replacements for59 in the tested family, including47 with a new58/67 label below9. Every accepted replacement rebuilds the complete wall55 quadrilateral on the S6 target face. Their target cones have25 finite rays and9 open ends. This is **not a certified improvement over the original63-wall,199-finite-vertex,10-end control**: the deleted old prism has not been rebuilt, and these normals have not been promoted to exact ambient geometry.

The [producer](../scripts/fresh_wall59_bond_replacement_20261018.py), [saved records](fresh_wall59_bond_replacement_20261018.json) and [independent saved-record checker](../scripts/check_wall59_replacement_saved_20261018.py) preserve the bounded experiment. The [companion successor](fresh_wall59_low_bond_companions_20261018.md) tests the changed five-edge cavity. No frozen predecessor is edited.

## Exact debit before interpreting the search

The authoritative control is the [exact replacement chamber](fresh_nonarithmetic_release49_62_replacement_20261017.md), with its [complete S6 boundary](fresh_nonarithmetic_release49_62_boundary_20261017.md). Its199 ambient finite supports give the following old-certificate debit:

|Released walls|Old finite supports lost|Compact2-face certificates touched|Compact3-face certificates touched|Lost supports containing8|
|---|---:|---:|---:|---:|
|59|7|6|1|0|
|59,60|7|6|1|0|
|59,60,61|7|6|1|0|

These are exactly the same seven finite supports and seven compact-face certificates, not just equal counts. Every old finite support involving60 or61 already contains59. The one compact S6 side touched is55. Thus releasing the two prism walls60/61 alongside59 costs no further old certificate in this ledger. The compact wall8 facet's old finite supports are untouched. This is a support-deletion statement; it does not by itself certify a replacement chamber or preservation against new numerical walls.

Without a replacement, deleting59 alone gives23 finite S6 rays and11 open ends. Deleting59/60, or59/60/61, gives23 finite rays and9 open ends. The latter decrease is already caused by deleting walls and must not be counted as closure progress.

## A different finite-angle family

The [original wall59](fresh_nonarithmetic_pi8_quad_completion_20261017.md) is right to S6+55, has native30 label5 and native58 label9. Its later prism walls60/61 are right to it. The new family replaces that bond rather than adding another cap constrained by it.

Work in the intrinsic Lorentz4 complement of S6, using the inherited seven-profile/native-offset implementation. The new wall67 has profile0, hence square2 and zero S6 contacts. With

    F=(55,30,58),  b=(0,-2cos(pi/m30),-2cos(pi/m58)),

solve F·v=b and v²=2. The frame has Lorentzian signature(2,1); its orthogonal line is positive. The linear projection and the two signs of its remaining norm therefore give both numerical norm branches. The prescribed finite corners are55/30/67 and55/58/67. No artificial30/58 adjacency is imposed.

The bounded native labels are

    m30 in {2,3,4,5,6,8,12},
    m58 in {2,3,4,5,6,7,8,9,10,11,12,16,24}.

All91 equations and182 branches are retained. Norm/contact residuals must be below2e-5. Every surviving branch is screened against every retained old native wall, with the inherited spacing-relative finite-angle gate; near-parallel ambiguity is rejected. The inherited strict common center chooses the outward branch. The root solve and these tolerance checks are numerical, not an interval proof or an all-label classification.

## What the retained walls force

Exactly91 branches have the wrong center sign. Among the91 forward branches:

-88 fail native admissibility only at60 and61;
-2 fail at60,61 and66;
-1 is native to every old wall except the intentionally removed59.

The last branch has labels(55,30,58)=(2,5,9) and reproduces the original59 normal within2e-6. It is the only accepted branch when retaining either both60/61 or61 alone. Thus the initial restricted experiment is preserved inside the final records, not discarded after releasing additional walls.

For illustration, keeping m30=5 and reducing m58 to8 gives native products approximately+0.05664 and+0.11026 with60 and61. At m58=10 the products are approximately-0.04068 and-0.07918, which lie in a nonnative angle gap. Only the sampled value9 restores their right contacts. This is evidence within the stated finite sample, not an all-parameter rigidity theorem.

After releasing all59/60/61,89 branches pass every remaining native wall and the strict center. Each complete S6 cone agrees with the independent positive-definite support incidence computation, and each has exactly the four finite wall55 rays

    (7,30,55), (7,55,58), (30,55,67), (55,58,67).

There are no other finite, ideal or exterior rays on55. This certifies the whole quadrilateral only at the stated numerical target-face level.

All89 replacements have the same five exterior S6 supports:

    (7,30,47), (30,47,52),
    (30,33,66), (30,66,67), (58,66,67).

Their nine one-open edges split into the persistent four-edge B region and the changed five-cycle

    30–67–58–66–33–30.

This actual boundary is the input for the companion experiment, not an assumed patch boundary. No full ambient finite-vertex census, exact normal certificate or finite-volume claim is made for a replacement.

## Verification scope

Generation96304 and separate fresh replay17641 both exited0/PASS. Both runs enumerate91 equations/182 branches, accept1/1/89 replacements in the three release scopes, and record zero numerical unresolved cases. Zero unresolved cases refers only to these sampled equations and the implemented numerical gates.

The saved checker independently reconstructs the exact debit from the199 old supports, verifies saved native values and profile0 norm residuals, all branch indices, old55 whole-face ray sets, boundary ledgers and saved cone/PD agreement. It does not reconstruct the physical coordinate engine, rerun a cone or prove the inherited finite-diagram classifier. The combined saved checker exited0/PASS (terminal c43eb7):23 source bindings,4 exact deletion scopes,182 branches,91 accepted replacement records,10 disk topologies,272 companion caps and7,200 pair indices across16 saved ledgers. The companion full generation9380 and fresh replay99757 also exited0/PASS. Root completed a full mathematical/source read of both producers, both reports and the saved checker: PASS within the stated bounded numerical scope. The read checked the identical7/6/1 deletion debit, retained native filtering, restored55 quadrilateral versus unreconstructed prism, actual five-cycle and separate B-region ledger, both disk duties including the deferred Z stage, complete7,200-pair accounting, and the helper's independent old199-support and saved-index checks. It did not reconstruct physical normals, rerun norm roots or cones, or prove the inherited spherical classifier. The independent sibling review checked only the two literal disk topologies.
