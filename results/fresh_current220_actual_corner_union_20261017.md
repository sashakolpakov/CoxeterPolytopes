# Joint actual-corner capacities on the current 220 represented Q classes

Testing the three previously separate geometry domains together gives **443 new capacity columns**, of which 61 violate the current exact witness. Keeping all 4,431 face masses fixed and increasing 46 vertex allowances repairs its ratio from 12.391408496497405 to

    893073974017141280 / 72281373619198223
        = 12.35552023017915... > 12.

A separate calculation starts from the frozen restriction retaining all 2,600 simplex rows and only these 220 Q rows. The same columns repair that restricted ratio from 12.025638499759506 to

    2600114203061924640 / 216842261971938119
        = 11.990809261150437... < 12.

Thus the new mixed geometry removes this particular restricted obstruction, while the full witness remains above 12. Neither repair computes an enlarged optimum or proves target-12 feasibility or a dimension bound. No optimizer is invoked.

The [source](../scripts/fresh_current220_actual_corner_union_20261017.py) and [certificate](fresh_current220_actual_corner_union_20261017.json) bind the [current finite model](fresh_raw23_mixed_full_model_20261017.md) and its [exact support restrictions](fresh_raw23_mixed_full_model_20261017_support.md). All frozen predecessors remain unchanged.

## Represented domain and actual multiplicities

The selected domain is exactly the 220 positive Q rows identified by `in_geometry_union` in the frozen support diagnostic: the union of the mixed-raw24, large-core raw178/296/297, and current raw23 mixed domains. It is not every positive class of the ten raw families. In particular, newly positive classes outside those named exports and ambiguous merged rows are not silently assigned representatives.

There are 149 explicit actual-provenance rows and 71 rows with exactly one saved source context. Their 1,860 actual corners give 33,526 distinct literal stars on 510 complete rank-13 finite vertex forests.

| Raw family | Classes | Actual corners | Literal orientations | Forests |
|---:|---:|---:|---:|---:|
|23|62|620|7268|340|
|24|16|128|1459|78|
|25|45|360|4404|262|
|39|6|24|130|18|
|87|16|128|5378|119|
|88|20|160|8186|128|
|89|21|168|848|143|
|178|15|120|645|104|
|296|12|96|5004|90|
|297|7|56|204|49|

Core orders are seven for raw23/24/25/39, eight for raw87/88/89, nine for raw178 and ten for raw296/297. Spectator ranks are respectively eight, seven, six and five. Raw23 has ten actual corners, raw39 four, and every other family eight. Each retained corner core is embedded componentwise into the complete finite forest, including every matching component copy and every internal weighted automorphism. Passive components have zero external contacts.

Full unordered external contact pairs identify literal stars. Coincident keys must agree in occupied quadrilateral pairs, exact squared heights and the native product of the two omitted roots. Complete orientation orbits are disjoint and each is attached to one Q class. Column coefficients count actual corners of the whole prism; they do not count orientation presentations.

## Exact fields and inherited pair ownership

The ground field remains K=Q(√5). The nonsquare native dotted parameters are

    tau23=(7+√5)/2, tau24=(23+9√5)/2,
    tau39=(41+17√5)/2, tau87=(9+3√5)/2,
    tau88=5+2√5.

Each parameter and each ratio of two distinct parameters is checked to be nonsquare in K. Other native dotted roots already lie in K. Any wall pair uses at most two of these radicals. Its immutable biquadratic field is selected from native edge encodings; all projections and signs are calculated in that field. No field-specific coefficient vector is transported between different fields.

At the common finite vertex, with incident positive-definite Gram V and complete external covectors p,q, the outward side fixes the actual product

    g = p^T V^-1 q - sqrt(h_p h_q),
    h_p = p^T V^-1 p - 2 > 0.

Native selfchecks recover the original omitted-root product. Identical covectors denote one wall and are skipped. The exact gate tests signs before squaring, membership of g² in K, algebraic integrality and the finite Coxeter alphabet where applicable. Distinct full stars occupying the same quadrilateral port pair are incompatible. Every pair not excluded by these necessary conditions remains allowed.

The three predecessor domains overlap. Reuse therefore has an explicit ownership rule: **the first predecessor map containing both stars owns that pair, including an unknown result**. A later map cannot replace an inherited unknown by an exclusion. The rule is applied recursively at every ancestor. Transport requires equality of the complete forest, unordered native covectors, occupied pairs, both heights and native omitted product; retained face and provenance records also agree. Inclusion-exclusion of map intersections checks the number of pairs owned by each map.

The saved `incompatibility_certificates` list contains only newly computed exclusions. A consumer must expand `inherited_pair_domains` in the saved order with this ownership rule before checking the full relation. Counting only the new list, or summing overlapping predecessor domains without ownership, is incorrect.

## Pair census and capacity certificates

All 4,785,958 literal pairs are covered: 4,321,030 inherited pairs and 464,928 new tests. The new tests give 448,198 exclusions:

| First applicable reason | New exclusions |
|---|---:|
|Positive product of distinct actual walls|161216|
|Distinct full stars sharing a quadrilateral pair|59568|
|Actual product in the open interval (−1,0)|77406|
|Squared native product outside K|60022|
|Nonintegral native two-cycle|89986|

Together, inherited and new certificates exclude 4,686,927 pairs. Of 2,044,552 cross-family pairs, 1,983,596 are excluded and 60,956 remain unknown. The specific cross-cohort comparison between raw23/87/88/89/178/296/297 and raw24/25/39 has 576,546 pairs on 265 common forests, with 556,036 excluded and 20,510 unknown.

Mutually incompatible literal stars give capacity-one orbit groups. Universally incompatible groups give clique columns. Proper colorings of the conservative literal compatibility graph give additional integer capacities: simultaneous actual stars form a clique, so at most one can occur per color. Every color class is checked pairwise incompatible. These certificates require neither an optimal coloring nor sufficiency of the retained compatibility tests.

After removing 1,956 exact predecessor tuples, the new export has 194 clique columns and 249 coloring columns. Of the 443 columns, 366 join the two cohorts; 59 of those violate the baseline. Mixed columns cause 44 additional vertex increases beyond the two within-cohort increases.

The largest increase is at forest 3002,

    A1 + A3 + A4 + H3 + I2(5).

For example, column 394 has cost one and coefficients

    Q26410 + 2 Q26474 + 2 Q26709 + Q28824
        + Q45649 + Q45845 + Q45848.

It joins raw23/24/39/88/89 and raises the vertex allowance from 138239475833/9900000000000 to 143530907/10000000000. The coefficients two are distinct actual corner occurrences. The complete literal clique witness is saved.

## Repair, preservation and limits

For each vertex the new allowance is the maximum of its baseline value and every new load divided by its capacity. Face masses and the objective numerator remain fixed. All 438,266 inherited finite-model variables remain valid, including 96,176 zero-budget signed columns and 52 positive-cost signed columns. Every primal charge variable is nonnegative; signed refers to its face coefficients. Its dual condition is a load inequality, not an equality. Increasing a vertex allowance relaxes positive-cost inequalities and leaves zero-budget inequalities unchanged.

The restricted calculation uses its own frozen baseline allowances and objective. Every new paid Q row is retained there with the same mass, so the new loads are identical, but the two normalizers and ratios are kept separate. Both exact repairs and their 46 updates are saved. Omitting the new cross-cohort columns leaves the full ratio at 893073974017141280/72072059559311473, still above 12.

Deduplication against historical exports does not assume that every such export was appended to the current matrix. All 4,266 distinct recursively bound historical capacity tuples are also checked directly at the restricted initial masses and allowances; all pass. Consequently omitting a historical duplicate loses no capacity in either repair.

This is a finite actual-corner capacity export at fixed masses. Older complete local-packing mixtures and implicit joint scans at different masses are not asserted to survive. The selected geometry does not cover all 1,831 positive Q rows, all possible arithmetic fields, or globally realizable polytopes. The full repaired ratio remains above 12, so this package does not justify an unchanged repeat of the global optimization.

## Verification

FINAL/FROZEN: final generation 27796 and fresh saved replay 19077 both terminate exit 0/PASS, in 108.13 and 79.04 seconds. They include the explicit 4,266-column restricted-baseline assertion and reconstruct the complete current geometry and all new pair tests. The separate saved-capacity checker 10670 also terminates exit 0/PASS. It imports no producer, expands the inherited relations, counts ownership by explicit pair enumeration rather than the producer's inclusion-exclusion, and checks all 1,860 occurrence records, 443 capacities, coefficients, loads, both repairs and ratios. This is independent saved-certificate verification conditional on the pair certificates, not an independent reconstruction of their physical projections or orientation census. All 41 direct source hashes, four local links and source syntax pass.

Root's complete narrow source read and construction's full source/math/provenance read pass the stated domain, immutable pair fields, ordered first-owner reuse, occurrence coefficients, conservative capacities and separate repairs. Construction identified the historical-dedup accounting point above, which is resolved by the explicit saved-load check. Neither review claims an independent orientation or 4.8-million-pair geometry census.

    python3 scripts/fresh_current220_actual_corner_union_20261017.py
    python3 scripts/fresh_current220_actual_corner_union_20261017.py --verify
