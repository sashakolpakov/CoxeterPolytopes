# Raw33 and raw295 mixed actual-corner capacities

Adding the 23 positive explicit/unique-context classes from raw33 and raw295 to the [current 220-class union](fresh_current220_actual_corner_union_20261017.md) produces **352 new exact capacity columns**. At the unchanged 4,431 face masses, 119 columns violate the baseline. Raising 72 vertex allowances repairs the full fixed witness from 12.35552023017915 to

    893073974017141280 / 75815153215391173
        = 11.779623678658464... < 12.

This improvement requires mixed geometry. Omitting the 208 columns joining the new families to the old ten-family cohort leaves

    893073974017141280 / 72504009047677973 > 12.

This is an exact repair of one fixed partial-model witness, not an enlarged optimum or target-12 feasibility certificate. No optimizer is invoked and no dimension bound follows. The [source](../scripts/fresh_raw33_295_mixed_actual_corners_20261017.py) and [certificate](fresh_raw33_295_mixed_actual_corners_20261017.json) preserve every frozen predecessor.

## Added actual domain

The baseline contributes exactly its 220 represented Q classes. The new cohort includes every positive raw33 or raw295 row in the frozen [support diagnostic](fresh_raw23_mixed_full_model_20261017_support.md) having an explicit actual provenance or exactly one saved source context. Ambiguous merged rows are not represented by arbitrary choices.

| New raw family | Core order | Prism dimension | Spectator rank | Classes | Actual corners | Literal orientations | Forests |
|---:|---:|---:|---:|---:|---:|---:|---:|
|33|7|3|8|10|60|492|43|
|295|10|4|5|13|104|1281|83|

The full union has 243 classes, comprising 160 explicit and 83 unique-context rows, 2,024 actual corners, 35,299 distinct literal stars and 527 complete finite rank-13 vertex forests. There are 96 forests shared by the new and old cohorts, with 172,364 literal cross-cohort pairs.

Every retained old face is reconstructed and compared with its frozen provenance and corner records. The two new cores undergo the same complete componentwise embedding: every whole retained component maps into every matching whole component of the finite vertex forest, including internal weighted automorphisms and all identical spectator copies. Core order and spectator rank are checked separately; no order-seven assumption is applied to raw295. Each dimension-d prism has exactly 2d distinct actual corners.

Full unordered external contact pairs identify literal stars. Occupied quadrilateral pairs, native omitted-root products and both exact positive squared heights must agree whenever a key coincides. Complete orbit domains are disjoint, and each orbit belongs to one Q class. Global coefficients count actual corners, not the number of orientation presentations.

## Dependent native radicals are embedded explicitly

The ambient arithmetic field is K=Q(√5). Raw33 has

    tau33 = 20+8√5 = 4(5+2√5) = 4 tau88.

Thus its positive native radical is **2√tau88**. Treating tau33 and tau88 as independent field generators, or interpreting a saved coefficient vector in a different field, would be incorrect.

The new immutable field decoder instead finds the unique basis square class containing each native squared contact t. If t/u=r² with positive r∈K and u one of 1, tau, sigma, tau·sigma, it encodes the acute native product as −r√u. Existence, uniqueness and physical positivity are checked. Field selection first replaces native parameters by fixed square-class representatives, then chooses at most two independent representatives for the wall pair. Harmless independent padding is used when fewer are needed.

Raw295 contributes tau295=3+√5. It is nonsquare in K: the rational candidates for the square of the rational coefficient of a putative square root are 5/2 and 1/2. Its ratios with the previous five representatives are also nonsquare, checked exactly. The exceptional rational ratio with tau87 is 2/3, which is neither a rational square nor five times one. The other nonsquare ratio norms give a direct obstruction. All six representatives and all fifteen pairwise ratio tests are saved through the exact domain assertions. A single wall pair still needs at most two representatives; no common six-radical coefficient representation is assumed.

Every self-projection, cross-projection and actual-wall test uses the field selected from that pair's native encodings. At the common finite vertex,

    g = p^T V^-1 q - sqrt(h_p h_q),
    h_p = p^T V^-1 p - 2 > 0.

Native omitted-root selfchecks fix the same outward future sign. Identical external covectors denote one actual wall and are skipped. The inherited exact routines test signs before squaring, fixed-K membership of the normalized two-cycle, algebraic integrality and the finite Coxeter alphabet where applicable. All unexcluded pairs remain compatible.

## Complete pair accounting and capacities

The baseline's 4,785,958 pair outcomes are transported by full native star keys. Its ordered first-owner rule is preserved recursively, including every inherited unknown outcome; no later ancestor overrides one. Transport checks the complete forest, occupied pairs, both heights and native omitted product. Saved field coefficients are never used as transport keys.

There are 203,025 new pair tests, with 197,993 exclusions:

| First applicable reason | New exclusions |
|---|---:|
|Positive product of distinct actual walls|66320|
|Product in the open interval (−1,0)|51850|
|Nonintegral native two-cycle|42147|
|Distinct full stars sharing a quadrilateral pair|22166|
|Squared native product outside K|15510|

The combined relation covers 4,988,983 pairs, with 4,884,920 exclusions. Among the 172,364 new-versus-old pairs, 167,578 are excluded and 4,786 remain unknown. The certificate's new exclusion list must be combined with the recursively ordered inherited domains; it is not the entire relation by itself.

Capacity-one orbit cliques and explicitly verified proper colorings give 145 and 207 new columns respectively, after removing 2,226 exact historical tuples. There are 208 mixed columns, 93 violated mixed columns and 55 additional vertex raises attributable to mixed columns. Unknown pairs remain edges of the conservative compatibility graph, so any actual simultaneous stars form a clique and every proper coloring gives a valid upper capacity. Optimal coloring and geometric realizability of allowed configurations are not claimed.

The largest allowance increase occurs at forest 982,

    3 A1 + A3 + H3 + H4.

Column 109 is a cost-one clique joining raw23/24/33/88/89/178/295:

    2 Q26328 + Q26402 + Q26448 + Q26563 + Q26637
      + Q27388 + Q45769 + Q45798 + Q60822 + Q68025.

It raises the vertex allowance from 272052307/10000000000 to 82294913/2500000000. Its wholeface coefficient two is two distinct actual corners, not an orientation multiplier. The literal clique and all other capacity witnesses are saved.

## Preservation and scope

The repair keeps all face masses fixed and takes the per-vertex maximum of the old allowance and each new load/capacity. Consequently all 438,266 inherited finite-model columns, including 96,176 zero-budget signed columns and 52 positive-cost signed columns, remain valid; the 443 baseline union columns also remain valid. All primal charge variables are nonnegative. Signed face coefficients do not turn a dual load inequality into an equality, so increasing its positive-cost vertex allowance is safe.

Older implicit local-packing covers and full joint scans were constructed at other masses and are not reused here. The selected 243-Q domain is still smaller than the full positive Q support. The package establishes new necessary geometric capacities and their exact fixed-witness effect. It does not establish complete field, forest, face or implicit-family coverage. An enlarged finite solve, if undertaken later, needs its own preserved-matrix and exact certificate checks.

## Verification

FINAL/FROZEN: generation 52402 and fresh saved replay 11884 both terminate exit 0/PASS in 74.82 and 84.44 seconds. The separate no-producer saved-capacity check 86982 also terminates exit 0/PASS. It expands all inherited relations with explicit pair ownership, verifies literal clique/coloring capacities and actual occurrence coefficients, and independently recomputes every new load, allowance maximum and ratio conditional on the saved pair certificates. It does not independently reconstruct the physical orientation or pair-geometry census. All 44 direct source bindings, four local links and source syntax pass.

Root's narrow complete source-delta and report read passes the positive-K square-class decoder, immutable field selection, retained 220-class domain, new provenance and six/eight-corner reconstruction, inherited repaired baseline, deduplication load guard, conservative capacities and monotone repair. The relevant inherited field routines were also read. This review does not claim an independent physical orientation census, pair census or saved-capacity replay. Any later review is recorded in a successor aggregate rather than changing these frozen files.

    python3 scripts/fresh_raw33_295_mixed_actual_corners_20261017.py
    python3 scripts/fresh_raw33_295_mixed_actual_corners_20261017.py --verify
