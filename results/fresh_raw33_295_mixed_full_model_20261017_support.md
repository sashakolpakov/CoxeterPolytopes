# Support of the enlarged raw33/295 finite-model dual

The new exact dual has 4,713 positive face masses and gives **12.382248317539402 > 12**. Restricting it to all positive simplex/edge rows and the 184 positive quadrilateral rows in the latest 243-class geometry export still gives **12.076081050570941 > 12**. Thus the remaining finite-model obstruction is not solely caused by quadrilateral rows outside that export.

The [source](../scripts/fresh_raw33_295_mixed_full_model_20261017_support.py) and [saved diagnostic](fresh_raw33_295_mixed_full_model_20261017_support.json) rebuild the [enlarged finite model](fresh_raw33_295_mixed_full_model_20261017.md), verify both exact certificates and recount all 795 appended capacity loads. No optimizer is invoked.

## Changed support and actual coverage

| Statistic | Count |
|---|---:|
|Positive maximal-simplex/edge rows|2720|
|Positive quadrilateral aggregate rows|1993|
|Positive vertex allowances|986|
|Newly positive face rows|1328|
|Previously positive rows now zero|1046|
|Changed common positive face masses|3385|
|Positive simplex rows outside the old complete feature domain|836|
|Quadrilateral rows with explicit actual provenance|377|
|Quadrilateral rows with one unique saved context|1214|
|Ambiguous merged quadrilateral rows left unexpanded|402|

The 2,720 simplex/edge rows comprise 1,534 maximal-simplex aggregates and 1,186 maximal-edge incidences. The latest named geometry export consists of exactly the 243 class identifiers in the raw33/295 package; its inherited 220-class export is a subset. There are 184 positive Q rows in the 243-class domain and 1,809 outside it. Of these 184, 175 belong to the inherited 220-class domain. This is a comparison with those precise exports, not a claim to enumerate every historical geometry class.

Only 292 current positive Q rows lie in the old joint-port export, while 1,701 lie outside it. The old complete simplex feature domain contains 1,884 current positive S rows. Neither membership nor an unchanged class identifier makes an old-mass packing mixture valid at the new masses. No old complete local cover or implicit joint scan is reused. Every inherited *exported* constraint remains present in the full finite matrix.

The 443 current220 columns have 439 positive loads, including 15 exactly tight positive loads. The 352 raw33/295 columns have 266 positive loads, including 14 exactly tight positive loads. All 795 inequalities pass exactly. This is a coefficient-load recount, not a new physical projection, orientation or pair-compatibility census.

## Exact restrictions

Each restriction retains all 2,720 positive S masses, sets excluded Q masses to zero, and recalculates vertex allowances using the entire saved matrix. All 96,176 zero-budget signed constraints and all 52 positive-cost signed constraints remain in the verification. No outside simplex is restored, no retained face mass changes, and each resulting vertex allowance is at most the original allowance.

| Retained Q domain | Kept Q | Dropped Q | Positive faces | Vertex allowances | Exact bound |
|---|---:|---:|---:|---:|---:|
|Latest 243-class geometry export|184|1809|2904|758|12.076081050570941|
|All explicit or unique-context Q rows|1591|402|4311|969|12.333600087212057|

The corresponding exact ratios are

    652755403518090600 / 54053579202106229,
    2666717340473877600 / 216215648441433661.

Both are exact obstructions for the full saved finite charge matrix. Neither establishes simultaneous global realizability or complete coverage by the implicit local and joint families. In particular, dropping all ambiguous merged Q rows does not remove the finite-model obstruction.

## Concrete next geometry

The largest explicit or unique-context cohorts outside the latest export are:

| Raw family | Classes | Actual corner occurrences | Demand-weighted mass | Native dotted square |
|---|---:|---:|---:|---|
|4|77|462|0.0271199631|(7+sqrt(5))/2|
|110|5|30|0.0224740650|15+6sqrt(5)|
|25, new class identifiers|17|136|0.0217542564|6+2sqrt(5)|
|125|10|60|0.0196324446|13+5sqrt(5)|
|76|19|76|0.0168557380|(7+sqrt(5))/2|
|57|33|132|0.0138104259|6+2sqrt(5)|
|205|9|54|0.0137951982|6+2sqrt(5)|

These are complete saved-metadata counts for the indicated current explicit/unique cohorts; they are not new capacity theorems. Raw110 is a small, concentrated next target: all five rows have unique saved contexts, its core has order eight, and the cohort gives only 30 actual corners. Its leading row is face 48317, source context 67927, mass 40746569/10000000000. Its dotted square is three times the raw88 parameter; this observation alone does not identify its radical with the raw88 radical over Q(sqrt(5)). A new exact square-class embedding and complete common-frame geometry would still be required.

Raw4 has the largest uncovered demand mass and shares the existing raw23 dotted parameter, but its 77 classes give a larger reconstruction domain. Its leading row is face 14792, unique source context 7841, mass 15063047/5000000000. Raw125 has ten unique-context classes and 60 corners; its leading row is face 51173, source context 70783, mass 2626973/500000000.

The diagnostic also identifies a different structural domain: 264 uncovered explicit/unique classes across 58 raw families have finite labels only 2, 3 and 5 and all dotted contact magnitudes in K=Q(sqrt(5)). They contribute 1,332 actual corner occurrences and demand-weighted mass 505833469/10000000000. Their common-field contact arithmetic could support one joint geometry treatment rather than separate arbitrary representatives. No projection or capacity reconstruction for that prospective union is asserted here.

There are 114 positive classes of previously treated raw families outside the latest exported class list, including 25 raw39, 24 raw24 and 17 raw25 classes. The exact represented-Q restriction above also leaves an obstruction. Consequently, either changed-mass complete local/joint separation on that restriction or new actual-corner compatibility is a substantive next step. Repeating the unchanged finite optimization is not justified by this diagnostic.

## Verification

FINAL/FROZEN: generation 80313 and fresh saved replay 61482 both terminate exit 0/PASS. Both rebuild the full matrix, verify both parent certificates, recount every appended load and verify both exact restrictions. A separate checker importing no project producer passes all 795 saved loads, both restriction objectives and normalizers, retained-mass identities and vertex-allowance monotonicity. All 15 direct source bindings, three local links and source syntax pass. The parent independently read the complete diagnostic delta and its restriction semantics; this is a source/composition review, not a second matrix or geometry implementation. No new optimizer, arithmetic field coverage or dimension claim is made.
