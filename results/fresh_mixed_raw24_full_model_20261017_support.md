# Support and exact restrictions of the new mixed-raw24 dual

The [new finite-model dual](fresh_mixed_raw24_full_model_20261017.md) remains above12 even after every ambiguously merged Q row is removed: the exact restricted ratio is **12.395014913347106**. By contrast, retaining only Q rows already covered by the exported mixed raw23/24/25/39 cohort lowers the ratio to **11.59248703405827**. These are exact restrictions of the current dual, not reoptimized values. They locate a concrete missing geometric domain without invoking another optimizer.

The [source](../scripts/fresh_mixed_raw24_full_model_20261017_support.py) and [certificate](fresh_mixed_raw24_full_model_20261017_support.json) rebuild the full435,793-column matrix, verify both original certificates and recount every one of the2,319 appended input exports. The diagnostic binds the new model hash `65e1790b7a9cae086d2829dc6084faf860cd2aed4924b5cb23a5d319df72d6fb`.

## Changed support and explicit coverage

The optimized dual has4,647 positive face masses and807 positive vertex allowances. Its support comprises2,747 simplex/edge rows and1,900 Q rows. Relative to the prior fixed-mass witness,1,109 faces are newly positive,1,521 have disappeared and all3,538 shared positive masses have changed.

Of the new simplex rows,2,122 appeared in the preceding complete packing-feature domain and625 did not. Even for the2,122 shared rows, old packing mixtures cannot be assumed to cover the changed marginals. Of the1,900 Q rows,138 are in the exported mixed actual-corner cohort and1,762 are outside it. The138 split as raw23:65, raw24:15, raw25:49 and raw39:9. Only378 Q rows are in the old joint-feature export;1,522 are outside it.

The complete source-context lists distinguish416 explicit Q rows,1,088 rows with exactly one saved context, and396 ambiguously merged rows. No arbitrary representative is selected for the latter. Exact active/tight counts for the appended input exports are:

| Exported family | Columns | Positive loads | Exactly tight positive loads |
|---|---:|---:|---:|
|Local simplex charges|78|76|9|
|Selected joint-port charges|68|41|0|
|Within-raw24 capacities|380|134|0|
|Mixed raw24 capacities|1793|1490|37|

Every load is recounted from the new alpha. All exported columns pass. This does not verify the full implicit packing or joint families at the changed masses.

## Two exact restricted witnesses

Both tests retain all2,747 current simplex/edge masses, remove the specified Q masses and recompute every vertex allowance from all435,793 columns. The source directly verifies that no retained face mass changes and no new simplex row is restored. Each resulting vertex allowance is at most its original allowance. All96,176 zero-budget signed transfers and52 signed columns with positive vertex cost are included.

| Retained Q domain | Kept Q | Dropped Q | Exact ratio |
|---|---:|---:|---:|
|Already exported mixed actual-corner cohort|138|1762|32134239971667960/2771988433306747 =11.59248703405827|
|Every explicit or uniquely saved Q row|1504|396|11452991257108520/923999796464609 =12.395014913347106|

The second exact dual proves that ambiguous merged provenance is not necessary for the surviving finite-model obstruction. The first test does not prove that all obstructions supported inside the mixed cohort are removed; it evaluates this particular restricted mass vector. Neither test establishes ambient realizability, complete implicit-family coverage or a dimension theorem.

## Concrete next geometric cohort

Ranking the uncovered explicit/unique Q classes by their exact contribution to the dual objective gives a manageable first target:

| Raw core | Canonical core | Classes (explicit+unique) | Core order | Actual corners | Demand-weighted mass |
|---|---:|---:|---:|---:|---:|
|178|158|18 (17+1)|9|144|0.0922718646|
|296|267|16 (14+2)|10|128|0.0918986118|
|297|268|13 (10+3)|10|104|0.0817261374|

Together these47 classes contribute376 actual corners and objective mass0.2658966138. Every finite label is2,3 or5. The dotted squares are6+2*sqrt(5) for raw178/297 and(7+3*sqrt(5))/2 for raw296, with positive square roots1+sqrt(5) and(3+sqrt(5))/2 respectively, both in K=Q(sqrt5). Thus their common-frame contact projections can use K arithmetic, without introducing another independent contact radical. Their timelike-height square roots still require exact sign and fixed-K two-cycle checks.

The top raw178 face60822 has provenance1355, context86657 and mass6431479/2000000000. The top raw296 face68072 has provenance1623, context102878 and the same mass. For raw297, face68119 has provenance1637, context102925 and mass28907257/10000000000. These literal identifications and every top cohort's complete context status are saved.

A useful next step is to reconstruct these whole finite corners and their complete forest embeddings, then test actual-normal compatibility within and between the47 classes. It must generalize the previous seven-root core code to orders9 and10, preserve whole-face corner multiplicities and handle the corresponding rank6/rank5 spectators. No such new compatibility conclusion or census is claimed by this diagnostic. The recommended domain is selected by actual uncovered mass and exact arithmetic simplicity, rather than an unchanged optimization retry.

## Verification scope

Generation55387 and fresh saved replay66078 both terminate exit0/PASS. The replay reconstructs the full matrix, both exact restrictions, all appended loads, source-context classification and literal cohort metadata without an optimizer. The primary finite-model generation, single solve and fresh saved verification are separately complete. No old-alpha full packing cover or joint scan is reused at the new masses.
