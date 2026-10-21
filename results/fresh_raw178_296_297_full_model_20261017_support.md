# Support of the enlarged large-core finite-model dual

The new optimized dual remains above12 after deleting every ambiguous merged quadrilateral row: its exact restricted ratio is **12.356546235052296**. Thus unresolved provenance choices are not the sole obstruction. In contrast, retaining only the170 positive quadrilateral classes in the two latest exported actual-corner domains gives **11.880391289435755**. Both restrictions retain all2,578 positive simplex/edge rows and pass every exported model constraint.

The [source](../scripts/fresh_raw178_296_297_full_model_20261017_support.py) rebuilds the [new finite model](fresh_raw178_296_297_full_model_20261017.md), verifies its exact certificates, and saves the support, every appended-column load and both restricted certificates in the [data](fresh_raw178_296_297_full_model_20261017_support.json). It performs no optimization, new projection census, or implicit packing/joint-family separation.

## Changed masses and coverage

The dual has4,410 positive face masses:2,578 simplex/edge rows and1,832 quadrilateral rows, with984 positive vertex allowances. Relative to its predecessor,1,128 faces become positive,1,365 disappear and3,281 common positive masses change. The preceding complete local-packing certificates cannot be reused at these changed marginals.

Of the2,578 positive simplex/edge rows,666 are outside the old complete local-feature domain. Only341 positive quadrilateral rows belong to the old selected joint-port export. This comparison identifies exported domains; membership alone does not certify a whole implicit family at the new masses.

The quadrilateral support contains397 explicitly split provenance rows,1,066 rows with exactly one saved context, and369 merged rows with several possible contexts. The diagnostic never selects an arbitrary representative of a merged row. It records each concrete context and canonical/raw core identifier when those are unambiguous.

The two named latest geometric domains are the181-class [mixed raw24/raw23/raw25/raw39 cohort](fresh_mixed_raw24_actual_corners_20261017.md) and47-class [raw178/raw296/raw297 cohort](fresh_raw178_296_297_actual_corners_20261017.md). Their union contains170 currently positive quadrilateral rows:128 from the first and42 from the second. This is explicitly a comparison with these two named domains, not a claim that all historically exported actual-corner classes lie in their union.

All469 newly added capacities pass exactly at the optimized masses;426 have positive loads and34 positive-load constraints are exactly tight. Every inherited exported model constraint, including the52 positive-cost signed incidence columns, is also checked.

## Exact restrictions

For each restriction, all positive simplex/edge masses and the retained quadrilateral masses are unchanged. Dropping other quadrilateral masses relaxes every zero-budget transfer because their positive entries occur on quadrilateral rows. Exact dual reconstruction introduces no additional simplex row and only decreases vertex allowances. The complete finite-model verifier checks the result; this is not merely a restricted-objective estimate.

| Retained quadrilateral domain | Q rows | Exact ratio |
|---|---:|---:|
| Two named latest actual-corner cohorts |170|197593316016201840 /16631886206635733 ≈11.8803912894|
| Every explicit or unique-context Q row |1,463|381669018512349360 /30888001489417327 ≈12.3565462351|

The full optimized ratio is335527348961076090 /27027016235751353 ≈12.4145168684. These restricted witnesses concern the same finite exported column family. Neither restriction establishes feasibility of target12, complete implicit coverage, or simultaneous geometric realizability.

## Next actual-corner domain

Ranking unambiguous positive classes outside the two named latest cohorts by their current demand times face mass identifies the following leading families. The ranking measures current objective contribution; it does not predict a future capacity improvement.

| Raw core | Core order | Classes | Explicit / unique | Actual corners | Demand mass |
|---|---:|---:|---:|---:|---:|
|87|8|36|12 /24|288|0.0534350580|
|89|8|26|19 /7|208|0.0477189162|
|88|8|27|20 /7|216|0.0459108018|
|295|10|13|11 /2|104|0.0330557640|

Raw87/88/89 together give89 classes,712 actual corners and demand mass0.147064776. All their finite labels are2,3,5. Their dotted squared products are respectively

    (9+3√5)/2,   5+2√5,   6+2√5.

The third has positive square root1+√5 in the fixed fieldK=Q(√5). The first has norm9, but its rational square-root test would require15/4 or3/4 to be a rational square; neither is. The second has nonsquare norm5, and the quotient has nonsquare norm9/5. Thus the first two generate independent square classes overK, so the existing pairwise biquadratic method is a suitable candidate for their mixed products. A useful next geometric task is one common-frame actual-corner compatibility domain across these three families, retaining all current explicit/unique classes and actual corner multiplicities. It should test mixed pairs directly and then inspect common finite forests with already covered families. No such orientation or compatibility census is claimed here.

Concrete high-mass rows are raw87 face45517, provenance954/context63363; raw89 face45881, provenance1019/context63727; and raw88 face45685, provenance988/context63531. The complete saved rows include every other current class. Raw295 is a separate13-class/104-corner candidate, including face68011, provenance1604/context102817.

The diagnostic also records a broader prospective domain of263 unambiguous classes/1,316 corners across53 raw families whose native finite contacts and positive dotted roots already lie inK. Its total demand mass is581186279/5000000000. This is a metadata selection only; it contains no new projection, orbit, pair-compatibility or capacity certificate. There are35 newly positive unexported classes of previously treated raw families, likewise explicitly listed rather than silently assumed covered.

## Verification

Generation20038 and fresh saved replay30076 both terminate exit0/PASS, each rebuilding the complete inherited model and checking its certificates, the469 loads and the two exact restrictions. Direct source hashes, Python syntax and local links pass. No optimizer is invoked by this package. The next geometric recommendation changes the represented actual-star domain; it does not repeat optimization of unchanged constraints.
