# Independent complete cone of the 48-wall batch

The [independent checker](../scripts/fresh_nonarithmetic_joint_regions_20261016_cone_review.py) and [certificate](fresh_nonarithmetic_joint_regions_20261016_cone_review.json) pass generation and fresh saved replay. The complete cone has **172 finite rays, no ideal rays, 2,109 hyperideal rays, and no opposite-sheet causal rays**. Its finite graph has **395 edges and 758 open finite-corner germs**. The separately checked strict escape proves that this is still an infinite-volume chamber.

The starting 1,168-ray list is inherited from the [independent complete 36-wall review](fresh_nonarithmetic_escape35_20261016_cone_review.md). The new normals and coherent 48-facet chamber are checked by the [independent batch review](fresh_nonarithmetic_joint_regions_20261016_batch_review.md). No new producer arithmetic is imported.

All twelve clipping stages are independently reproduced:

| Added wall | Finite rays | Hyperideal rays | Old rays exactly on the added wall | New crossing rays |
|---:|---:|---:|---:|---:|
|36|101|1,219|29|279|
|37|110|1,381|22|304|
|38|115|1,514|28|246|
|39|132|1,641|48|341|
|40|136|1,724|0|145|
|41|148|1,772|4|100|
|42|151|1,841|0|111|
|43|153|1,913|0|111|
|44|155|1,941|0|60|
|45|158|2,016|1|113|
|46|170|2,065|5|91|
|47|172|2,109|0|64|

## Independent adjacency and zero handling

Full dimension and pointedness follow from the inherited rank-ten old covectors and the independently checked common strict center. For two distinct current extreme rays, their sum lies in the relative interior of their common-zero face. Its dimension is therefore ten minus the rank of their common active covectors. They are adjacent exactly when that rank is eight.

Every accepted negative–positive crossing pair has an independently certified nonzero eight-by-eight minor of its common active covectors. Fewer than eight common walls cannot suffice. Five remaining candidates at wall 39 fail adjacency: for each, the certificate supplies three linearly independent rays annihilated by all common covectors. This proves common rank at most seven, without using the producer's exact-two-rays-in-the-common-face acceptance test. All other candidates with at least eight common walls pass the rank-eight test.

For incoming contact values `a<0<b`, the review constructs

```
w = (b*u-a*v)/(b-a).
```

Its coefficients are strictly positive and sum to one. Its old zero set is exactly the intersection of its parents' zero sets; the new-wall contact vanishes identically. The new equation is independent of the old rank-eight support, so its active rank is nine. The code asserts mask uniqueness; it never silently discards an additional ray.

The 137 exact zero incidences are checked with the frozen independent nested-radical engine. Such rays are retained with the incoming wall added to their active mask. Strict signs use outward rational intervals, starting from the independently certified old normalized-coordinate enclosures. Exact fallback expressions may use different positive scales from the interval convex expressions: the crossing construction is projectively invariant under independent positive rescaling of both parents. This avoids any need to assume the two arithmetic representations have matching magnitudes.

The ordinary exhaustive halfspace-clipping theorem now propagates completeness at each stage. All 109,488 final ray–wall inequalities are checked, including every inactive strict inequality, and every final norm is strictly signed. The final active-set census agrees literally with the producer.

## Finite graph and scope

Every finite vertex has nine incident walls. The review independently checks a rank-eight common active span for each finite edge and counts 395 edges. Simplicity ensures that every finite edge appears this way. Hence the remaining number of finite-corner germs is exactly

```
9*172 - 2*395 = 758.
```

The proof is complete for the specified 48 walls, with the independently complete 36-wall cone as an explicit starting dependency. It proves no finite-volume completion and supplies no new dimension bound. Completeness of the old cone's sixteen compatible polar proposals is established separately by the [all-label polar review](fresh_nonarithmetic_joint_regions_20261016_polar_gate_review.md).
