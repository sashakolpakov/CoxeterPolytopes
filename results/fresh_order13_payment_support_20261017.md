# Exact order-13 payment support and the remaining B/D obstruction

The current rank-12 envelope remains **6**, corresponding to all78 pairs in a factor of order13. The complete selected-pattern diagnostic plus the immediate affine alternatives of the existing D14 transfer require all78 pairs for both the affine B and affine D targets. The analogous unions contain only65 pairs for affine A and73 for affine C. Lowering A or C alone therefore does not lower the worst rank-12 cost.

This package is a diagnostic, not a new payment theorem or dimension bound. In particular, an absent preferred-map demand is not by itself a proof that a pair can be removed from the complete face-payment rule.

The [source](../scripts/fresh_order13_payment_support_20261017.py) and [exact saved data](fresh_order13_payment_support_20261017.json) reconstruct the complete new target-order domain from the frozen [2,614-context classification](fresh_distance8_high_overlap_20261017.md). They reuse the same exhaustive induced-map algorithm as the [established small-factor census](fresh_ideal5_small_factor_payments_20261017.md), checking both labelled edges and nonedges. No temporary exploratory JSON is an input.

## Complete selected-pattern domain

A selected component can embed in an affine factor of order13 only if its proved full-component bound is at least13 and its displayed marked component has at most13 nodes. Exactly801 of the current2,614 contexts satisfy these conditions:

|Proved full-component bound|Eligible patterns|
|---:|---:|
|13|556|
|14|108|
|15|87|
|16|50|

These include606 displayed components of order13,160 of order12,33 of order11 and two of order10. The larger-bound entries cannot be discarded merely because the target is smaller. All injective induced labelled maps of every eligible pattern into each of the four order13 classical affine types are enumerated. The standard affine rank is12 in each case.

|Target|Basic ideal5/R5 pairs|Selected marked-pair union|Basic plus selected plus immediate union|Contexts with a map|Emitted maps|
|---|---:|---:|---:|---:|---:|
|A|65|39|65|67|1742|
|B|73|9|78|107|240|
|C|72|5|73|13|30|
|D|74|8|78|305|2344|

The basic rule pays distance at most5, or distance6 through9 with a mark within5 of a branch or nonsimple-bond feature. The selected union is not assumed disjoint from that basic set. The saved certificate contains all maps, marked images and whole marked-core keys.

For affine A, the13 pairs absent from this union are exactly its distance-six pairs. Their removal would still require the universal cycle13 replacement in the short8-first stage. The existing cycle theorem starts at larger order, so this note does not silently extend it. For affine C the five absent pairs are `(0,10),(0,11),(0,12),(1,12),(2,12)` in the saved canonical numbering. No new reduced C rule is promoted here either.

## The extra B and D orbits

For the B target, the fork center is0 with short leaves11,12; the spine is0–1–…–10 and the bond9–10 is doubled. Its five extra pairs form three symmetry orbits:

|Representative|Distance|Orbit size|Selected supporting contexts|
|---|---:|---:|---:|
|0,10|10|1|2|
|9,11|10|2|23|
|10,11|11|2|4|

The singleton0,10 is required by contexts160 and1206. They are the two orientations of L10 with marks{3,4} joined to L25 with marks{0,3}, sharing L10 node5 with L25 node9. Their preferred13-node component is already affine B12. Unlike the completed order15 reduction, it is not an abstract affine completion of an untouched finite component. Removing this singleton needs an actual new local transfer.

The distance11 orbit has contexts166,170,1744,2137, also orientations and marked choices of this same radical L10/L25 core. The remaining B orbit has23 contexts across several core families. The source saves every literal key rather than classifying these cases by counts alone.

For the D target, the fork centers are0 and8, with leaves9,10 at0 and11,12 at8. The only extra orbit consists of the four cross-fork leaf pairs, each at distance10. It is required by127 selected contexts. Thus the B obstruction is small in one orbit, while the D obstruction already has a substantial complete support domain.

## Immediate alternatives in the existing D14 transfer

The [universal affine D14 long-pair transfer](fresh_affine_d14_long_pair_transfer_20261017.md) has an original factor of order15. Deleting a fork center leaves finite D12 plus two isolated leaves. Its replacement profiles `inner4` and `inner_next3` create connected affine factors of order13, already paid in that proof by the whole-small-factor rule.

This package reconstructs both literal13-node labelled diagrams, matches their frozen positive-kernel records, and enumerates all marked images for each of the three possible opposite marks. They give:

|Profile|New affine target|Marked-pair roles|
|---|---|---|
|inner4|B12, order13|doubled terminal against the far fork center or either short leaf|
|inner_next3|D12, order13|one new short leaf against the far fork center or either short leaf|

The B alternatives therefore include both the singleton0,10 and the two distance11 pairs. The D alternatives include all four cross-fork leaf pairs, as well as already-basic distance-nine pairs. All these immediate images are already in the unions in the preceding table, so the numerical counts do not rise. Their provenance is nevertheless different: proving a new replacement only for the selected L10/L25 or other core contexts would not automatically discharge these universal D14 alternatives.

The lower D13 transfer has immediate factors of order12, so it does not introduce another order13 demand. The previously sealed order16/17 transfers end in finite payments or short paths. The diagnostic does not reprove those frozen conclusions.

## A precise next local obligation

For a possible universal affine D12 long-pair theorem, the old profile formulas suggest a useful extreme case. The center/center, center/leaf and leaf/leaf resistances would be8,10,12. Their extreme costs sum to12,14,16, respectively. Conditional on forcing both heights to1/2, the zero-determinant equation requires unscaled spectator energy1. The frozen energy-one classification is therefore directly relevant to this extreme branch.

It does **not** close the entire theorem. Merely reusing the earlier nonextreme bounds, center at most3 and leaf at most4, gives normalized gaps `-1/8`, `-1/10` and `0` in these three end types. Thus there are additional nonextreme profiles and equality cases to analyze; neither the old positive-gap proof nor a height1/2 reduction may be imported without that work. For reducing the four extra D pairs alone, only the leaf/leaf branch is required, so its zero-gap cases are the narrower relevant obligation. The negative center cases matter only to a stronger all-nine-pair theorem. The certificate records these exact conditional comparisons to prevent treating the extreme spectral match as complete coverage.

The formal dimension74 boundary profile `(1,12,12,12,12,12,12)` explains why this is a meaningful obstruction: it uses six rank-12 factors. Nevertheless this report changes no factor cost and draws no new scalar or dimension conclusion.

## Verification scope

Generation and a separate fresh saved replay both completed with terminal PASS for all801 eligible patterns,4356 emitted maps, six literal immediate-profile/mark records and the exact orbit counts. Final checks passed for Python syntax, all17 direct source bindings and all5 report links. The checker validates completed input statuses and every recorded direct source digest before using them.

The root agent completed a separate full source and written-mathematics review: the eligible domain includes the larger component bounds, the basic/selected/immediate unions remain distinct, the six literal D14 interfaces are accounted for, and the B/D orbit and conditional D12 gap statements have the stated scope. That review did not repeat the4356-map enumeration. Neither this package nor that review repeats the82,197 overlap construction or its all-label geometric proof. The map enumeration is complete within the explicitly defined selected-component domain; it is not a classification of all actual chambers or all possible face stars. No new local payment rule or dimension bound is asserted. The source, data and report are final and frozen; no frozen predecessor was edited.
