# Current raw23 mixed actual-corner capacities

The current raw23 classes, tested jointly with raw87/88/89 and raw178/296/297, give 739 new capacity columns after removing 1,678 exact duplicates of frozen exports. Of these, 272 violate the saved witness. Keeping all 4,410 face masses fixed and raising 148 vertex allowances changes its exact ratio from 12.200030329472785 to

    1342109395844304360 / 116030537828718707
        = 11.56686352540649... < 12.

The improvement comes from actual mixed geometry. If the 664 columns joining raw23 to the other six families are omitted, the ratio is instead 12.20003032807404, still above 12. This calculation repairs one fixed partial-model witness. It does not determine the enlarged optimum, establish target-12 feasibility, or prove a dimension bound. No optimizer is invoked.

The [source](../scripts/fresh_raw23_mixed_current_corners_20261017.py) and [certificate](fresh_raw23_mixed_current_corners_20261017.json) bind the [previous mixed witness](fresh_mixed_raw87_large_core_corners_20261017.md) and the [frozen finite model](fresh_raw178_296_297_full_model_20261017.md). All predecessor files remain unchanged.

## Current represented domain

Every positive class in the seven named raw families having an explicit provenance or exactly one saved source context is included. No ambiguous merged row is assigned a representative. The current support contains 198 such classes: 149 explicit and 49 unique-context classes. The represented domain agrees exactly with the corresponding entries of the frozen support diagnostic.

| Raw core | Core order | Spectator rank | Classes | Actual corners | Literal orientations | Forests |
|---:|---:|---:|---:|---:|---:|---:|
|23|7|8|66|660|7226|369|
|87|8|7|36|288|9530|239|
|88|8|7|27|216|9248|176|
|89|8|7|26|208|1042|173|
|178|9|6|18|144|764|124|
|296|10|5|14|112|5233|105|
|297|10|5|11|88|324|77|

The union has 1,716 actual prism corners, 33,367 distinct literal orientations and 568 complete rank-13 finite vertex forests. Raw23 shares 266 forests with the other cohort, giving 996,792 literal cross-cohort pairs. The six-family baseline is unchanged. Relative to the older raw23 domain in the [mixed raw24 package](fresh_mixed_raw24_actual_corners_20261017.md), the current raw23 domain adds faces 26262, 26275, 26308 and 26412 and drops nine zero-mass faces; the exact lists are saved.

Raw23 consists of five-dimensional canonical dotted prisms with ten actual corners; the other six families have dimension four and eight corners. For each corner, deleting its two omitted roots and adjoining its entire spectator forest gives thirteen incident walls. Every retained whole component is embedded into every matching whole component of the target forest, including all weighted automorphisms and identical-component choices. Full unordered external covector pairs identify literal stars. Coincident keys must also agree in occupied quadrilateral pairs, native omitted-root product and exact squared heights. The complete orientation orbits are disjoint and each belongs to one face class. Column coefficients count actual prism corners, never the number of orientations.

## Exact pairwise fields and actual walls

Let K=Q(√5). The only nonsquare dotted parameters needed are

    tau23=(7+√5)/2,
    tau87=(9+3√5)/2,
    tau88=5+2√5.

Each is nonsquare in K, and every ratio of two is nonsquare. These statements are checked exactly. Their norms are respectively 11, 9 and 5; the three ratio norms are nonsquare rationals. The other dotted roots are already in K:

    sqrt(tau89)=sqrt(tau178)=sqrt(tau297)=1+√5,
    sqrt(tau296)=(3+√5)/2.

Each external contact covector uses at most one of the three nonsquare radicals. A pair of covectors therefore lies in one of the three exact biquadratic fields K(√tau_i,√tau_j). The immutable field is selected from the native contact encodings for that specific pair. Harmless independent padding is used if fewer than two radicals occur. Every projection, native-product selfcheck, sign comparison and arithmetic test uses that selected field. A coefficient vector from one field is never interpreted in another. Saved radical expressions name their two field generators explicitly. Pairwise incompatibility suffices to rule out a simultaneous collection; this argument does not require constructing a global three-radical coefficient representation.

In the common complete finite vertex frame with positive-definite incident Gram V, an external norm-2 wall with full covector p has squared timelike height

    h_p = p^T V^-1 p - 2 > 0.

The outward side fixes its future sign. Two distinct external walls thus have forced product

    g = p^T V^-1 q - sqrt(h_p h_q).

Selfchecks reproduce each native omitted-root product. Identical full covectors denote the same actual wall and are skipped. All supported finite components occur in the projection; passive components contribute zero.

The frozen exact field routines determine the sign of g before squaring. They then impose the necessary fixed-ground-field condition g²∈K, algebraic integrality, and the finite Coxeter alphabet where applicable. If sqrt(h_p h_q) lies in the selected biquadratic field, the unique positive physical square-class representative is used. If it lies outside, a nonzero projection makes g² leave that field and hence K. Two distinct complete stars sharing an occupied quadrilateral port pair are also incompatible because that actual two-face has a unique complete wall star. All unexcluded pairs remain allowed.

## Recursive reuse and new exclusions

Every current face and orientation orbit is reconstructed. A retained face's source context, provenance and old corner fields must agree with its frozen record after canonical JSON serialization. A pair relation is transported only after matching the full forest, full unordered covectors, occupied pairs, both exact squared heights and native omitted product. This comparison uses native encodings, not field-dependent projections or old masses.

The immediately preceding mixed certificate itself references older pair domains. The new producer recursively expands all such domains, validates each injective star map and its exact reused-pair/exclusion totals, and checks every bound source hash. An inherited unknown pair remains unknown. There are 3,783,066 reused pair tests, including 3,738,906 exclusions, and 1,009,380 genuinely new pair tests.

The new tests exclude 978,090 pairs:

| First applicable reason | New exclusions |
|---|---:|
|Distinct full stars sharing a quadrilateral pair|109833|
|Positive product of distinct walls|296078|
|Product in the open interval (−1,0)|317103|
|Nonintegral normalized two-cycle|218720|
|Squared native product outside K|36356|

Across the full union, 4,716,996 of 4,792,446 pairs are excluded and 75,450 remain unknown. Among raw23 cross-cohort pairs, 965,502 are excluded and 31,290 remain unknown.

**Saved-certificate convention:** `incompatibility_certificates` contains only newly computed exclusions. Consumers must recursively expand every `inherited_pair_domains` entry, with its saved native-key transport, before checking a capacity or reusing a pair relation. The new exclusion list alone is not the complete relation. This convention avoids duplicating millions of frozen pair records while retaining their exact provenance.

## Capacity columns and repair

An orbit whose distinct literal stars are mutually incompatible has capacity one. Universally incompatible such orbits form capacity-one cliques. More generally, a proper coloring of the conservative literal compatibility graph bounds simultaneous stars by its number of colors. Every color class is checked pairwise incompatible. Unknown pairs are retained as compatibility edges; no optimal-coloring or geometric-realizability claim is required.

The 739 new columns comprise 517 clique columns and 222 proper-coloring columns. Their costs are 528 ones, 157 twos, 25 threes, 18 fours, seven fives and four sixes. Exact tuples of vertex, cost and all face coefficients are deduplicated against every recursively bound predecessor export. All omitted old tuples are directly checked against the current baseline allowance.

There are 664 mixed columns, of which 269 violate the baseline. They cause 145 additional vertex raises beyond the three tiny within-cohort raises. The largest raise occurs at vertex 904, with forest

    3 A1 + A2 + 2 H3 + I2(5).

For example, column 321 gives a capacity-one clique on six actual occurrence classes: raw23 faces 26304, 26362 and 26392 together with faces 45517, 45625 and 45713 from raw87/88. The exact load is 93609129/10000000000, exceeding the old allowance 65028699/10000000000. Its literal orbit and all pair exclusions are saved. Thus the improvement is an actual common-vertex constraint, not a new abstract port-family bound.

The repaired allowance at each vertex is the maximum of its old allowance and each new load divided by its positive cost. The objective numerator remains 122903790828233/9900000000000. The new normalizer is 116030537828718707/108108000000000000, giving the ratio displayed above.

All 436,262 exported columns of the frozen finite model, the previous 1,002 raw87 columns, the 263 previous mixed columns, and the other bound predecessor exports remain valid. The 96,176 zero-budget signed columns depend only on unchanged face masses. The 52 signed positive-cost columns satisfy load≤cost·mu because their primal variables are nonnegative; increasing mu preserves them. Here “signed” describes face coefficients, not free variables. No complete implicit local or joint-port cover from different face masses is reused. Current full implicit coverage, other actual-corner families and broader arithmetic-field/domain coverage remain unfinished.

## Verification and review

FINAL/FROZEN: generation 57196 and separate fresh saved replay 50446 both terminate exit 0/PASS, in 143.44 and 182.12 seconds respectively. The no-producer saved-capacity consistency walk 57887 also terminates exit 0/PASS. It expands all inherited pair relations and verifies 1,716 occurrence records, all 739 clique/coloring capacities, exact loads, 148 allowance increases and the repaired ratio conditional on the saved pair certificates. It does not independently reconstruct their physical projections or orientation census. All 38 direct source hashes, five local links and source syntax pass.

Construction's narrow source/math delta review passes immutable per-pair field selection from native contacts, K-valued heights, recursive ancestor transport and completeness counts, current provenance and variable-dimensional corner multiplicities. This review does not duplicate the arithmetic census or inherited field implementation. Any subsequent review is recorded in the successor aggregate rather than changing this frozen report. No predecessor is edited and no optimizer is run.

    python3 scripts/fresh_raw23_mixed_current_corners_20261017.py
    python3 scripts/fresh_raw23_mixed_current_corners_20261017.py --verify
