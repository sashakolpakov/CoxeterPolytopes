# Support of the enlarged current-raw23 finite-model dual

The new exact dual has 4,431 positive face masses and still gives **12.391408496497405 > 12**. Restricting it to all simplex rows and 220 quadrilateral rows in the three named actual-corner domains gives another exact obstruction, **12.025638499759506 > 12**. Thus the remaining obstruction cannot be attributed solely to quadrilateral rows outside those domains.

The [source](../scripts/fresh_raw23_mixed_full_model_20261017_support.py) and [saved diagnostic](fresh_raw23_mixed_full_model_20261017_support.json) rebuild the complete [enlarged finite model](fresh_raw23_mixed_full_model_20261017.md), verify both exact certificates and recount all 2,004 appended column loads. No optimizer is called and no old-mass implicit cover is reused.

## Changed masses and coverage

| Support statistic | Count |
|---|---:|
|Positive maximal-simplex/edge rows|2600|
|Positive quadrilateral aggregate rows|1831|
|Positive vertex allowances|761|
|Newly positive face rows|1225|
|Previously positive rows now zero|1204|
|Changed common positive face rows|3206|
|Positive simplex rows outside the old complete feature domain|740|
|Quadrilateral rows with an explicit actual provenance|384|
|Quadrilateral rows with one unique saved context|1046|
|Ambiguous merged quadrilateral rows left unexpanded|401|

The three named geometry domains are the older mixed-raw24 domain, the order-nine/ten raw178/296/297 domain, and the current raw23/87/88/89/178/296/297 union. Their union contains 220 positive Q rows; 1,611 positive Q rows are outside it. This is an explicit named-domain comparison, not a claim that the union exhausts every historical corner class. Of the current Q rows, 149 belong to the latest raw23 union and 316 belong to the old joint-port feature export.

Complete local packing mixtures and full implicit joint scans at older masses are not current certificates. All common positive face masses changed; 740 positive simplex classes are outside the old complete feature domain as well. Their features and marginals would need reconstruction for a new full-family claim. The finite matrix itself retains every inherited exported column and all new columns.

The 1,002 raw87 exports have 719 positive loads, the 263 mixed-large exports have 243, and the 739 current raw23 exports have 707. All pass exactly at the new masses. Respectively 3, 2 and 12 positive loads are exactly tight in the rational repaired dual. This is a load recount, not a new pair-geometry or orientation census.

## Exact restrictions

For each restriction, all 2,600 simplex masses are retained and the excluded Q masses are set to zero. The complete matrix verifier recomputes every vertex allowance, including all signed vertex-cost loads, while retaining every zero-budget transfer. No simplex outside the retained set is restored, and no retained face mass changes.

| Retained Q domain | Kept Q | Dropped Q | Exact bound |
|---|---:|---:|---:|
|Three named actual-corner domains|220|1611|12.025638499759506|
|All explicit or unique-context rows|1430|401|12.33809262465147|

The corresponding exact ratios are

    2600114203061924640 / 216214232875362329,
    889230048650259680 / 72071921949554903.

The first restriction has 2,820 positive face rows and 733 vertex allowances; the second has 4,030 face rows and 753 allowances. Both certify obstructions only for the full saved finite charge matrix. Neither supplies complete implicit local coverage or simultaneous global geometric realizability.

## Concrete next geometry

The highest-demand explicit/unique Q cohorts outside the three named domains are:

| Raw family | Classes | Actual corners | Demand-weighted mass | Native dotted square |
|---|---:|---:|---:|---|
|33|10|60|0.0407114379|20+8√5|
|295|13|104|0.0335062710|3+√5|
|4|68|408|0.0278838315|(7+√5)/2|
|25, newly represented classes|11|88|0.0206446620|6+2√5|
|125|8|48|0.0182372703|13+5√5|

These are exact metadata counts and masses, not newly proved capacities. Raw33 is especially economical: its core has order seven, its ten represented classes give only sixty actual corners, and

    20+8√5 = 4(5+2√5) = 4 tau88.

Its dotted contact therefore lies in the existing raw88 square class with positive coefficient two. The current field helper requires exact native-parameter equality, so this dependent-square-class embedding must be implemented and checked explicitly rather than treating two equivalent radicals as independent. Raw295 has order ten and finite labels only 3 and 5; its thirteen classes are another small concrete domain. Existing contact, height and whole-component orientation methods are relevant, but literal provenance and multiplicities must be reconstructed before exporting a new constraint.

The leading raw33 row is face 27388, unique source context 31659, mass 24492627/5000000000. The leading raw295 row is face 68011, explicit provenance 1604 and source context 102817, mass 1713297/1000000000. Complete ranked rows, dotted parameters and prospective common-K contact domains are saved in the JSON.

There are also 67 newly positive classes of already treated raw families outside the named export union, including eleven raw25 classes and seventeen each from raw24 and raw39. Completing mixed compatibility on the restricted 220-Q witness is another direct structural target. These observations justify new geometry or changed-mass full-family separation; they do not justify repeating the unchanged finite optimization.

## Verification

FINAL/FROZEN: generation 25877 and fresh saved replay 73441 both terminate exit 0/PASS. Both rebuild the full model, verify its exact primal and dual, check all appended loads and recalculate the two restrictions. A separate small check confirms that restricted masses are unchanged, vertex allowances do not increase, both exact objective/normalizer ratios agree and both target flags are correct. All 18 direct source bindings, three local links and source syntax pass. This verification does not reconstruct old implicit packing cones, new pair incompatibilities or a broader arithmetic field. No second optimizer is invoked.
