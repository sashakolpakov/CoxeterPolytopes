# Mixed actual-corner capacities across the two new cohorts

Joint testing of raw87/88/89 against raw178/296/297 gives263 new capacity columns after removing1,263 exact duplicates of previous exports. Forty-six new inequalities violate the saved witness. All are cross-cohort constraints. Keeping its4,410 face masses fixed and raising33 vertex allowances changes the exact ratio from12.301640818329087 to

    1342109395844304360 / 110008693388412467
        =12.200030329472785... >12.

This is an exact repair of a fixed partial-model witness. No optimizer is invoked, the enlarged finite optimum is not determined, and no dimension theorem follows. The [source](../scripts/fresh_mixed_raw87_large_core_corners_20261017.py) and [certificate](fresh_mixed_raw87_large_core_corners_20261017.json) retain every previous inequality through unchanged face masses and monotone vertex allowances.

## Complete represented domain

The baseline is the frozen [raw87/88/89 witness](fresh_raw87_88_89_actual_corners_20261017.md), built on the [latest finite model](fresh_raw178_296_297_full_model_20261017.md). The new union includes every current positive class in these six raw families with explicit or unique saved provenance. No ambiguous row is assigned a representative.

| Raw core | Core order | Spectator rank | Classes | Actual corners | Literal orientations | Forests |
|---:|---:|---:|---:|---:|---:|---:|
|87|8|7|36|288|9530|239|
|88|8|7|27|216|9248|176|
|89|8|7|26|208|1042|173|
|178|9|6|18|144|764|124|
|296|10|5|14|112|5233|105|
|297|10|5|11|88|324|77|

There are132 classes,91 with an explicit provenance map and41 with exactly one saved context. Their1,056 actual prism vertices yield26,141 distinct literal stars across465 complete rank13 vertex forests. The two cohorts share135 forests and396,080 literal cross pairs. Relative to the earlier47-class [large-core package](fresh_raw178_296_297_actual_corners_20261017.md),42 classes remain positive, face60847 is newly included, and five old classes have zero mass. The43 current large classes are reconstructed from the current support, not copied as the old47-class domain.

Every aggregate is a four-dimensional canonical dotted prism with eight actual corners. Each corner omits two roots from its complete core and includes all spectator components, giving thirteen incident walls. Every whole retained component is embedded into every whole matching component of the target forest through all weighted automorphisms and identical-component choices. Full unordered external contact pairs determine literal stars. Occupied quadrilateral pairs, native omitted-root products and origin multiplicities are checked when keys coincide. Complete orientation orbits are disjoint and belong to one face class. Coefficients count actual prism corners; orientation counts are not substituted for those multiplicities.

## One exact field for all six families

Let K=Q(√5). The three new large-core dotted roots lie in K:

    tau178=tau297=6+2√5=(1+√5)^2,
    tau296=(7+3√5)/2=((3+√5)/2)^2.

They therefore embed in the already proved field

    E=K(√tau87,√tau88),
    tau87=(9+3√5)/2, tau88=5+2√5.

The raw89 dotted root is also1+√5. The source checks all positive K roots and the three nonsquare conditions for tau87,tau88,tau87/tau88. No third radical and no independent sign choice are introduced by combining the cohorts.

In one complete finite vertex frame with Gram V, an external wall of norm2 and contact covector p has strictly positive squared timelike height h_p=pᵀV⁻¹p−2. Its outward side fixes the positive future coefficient. Two distinct walls therefore have the forced product

    pᵀV⁻¹q − sqrt(h_p h_q).

All reconstructed self-products reproduce the prescribed native omitted-root product. Identical full covectors represent one actual wall and are skipped, rather than rejected by their positive self-product. Components with no supported contacts contribute zero to the projection; every active component is retained.

The frozen exact biquadratic arithmetic determines signs before squaring and applies the fixed-ground-field two-cycle condition. If a height root is in E, its unique positive physical square-class representative is used; otherwise a nonzero E projection makes the squared product leave E and hence K. Algebraic integrality and the finite Coxeter alphabet are tested only after K membership. These are necessary conditions for the compact-arithmetic K model. Sharing an occupied quadrilateral port pair also excludes two distinct full stars, since that actual two-face has a unique full wall star. Unknown pairs remain allowed.

## Reuse of exact old pair domains

All current corners and orientation orbits are freshly reconstructed. For every retained face, its source context, provenance and complete corner record must agree with the frozen predecessor after canonical JSON serialization. At each forest an inherited literal star is reused only when its complete unordered contact pair, occupied pairs, both exact squared heights, native omitted product and full forest agree.

These identities allow the entire old pair relation to be transported into the current frame, independently of old face masses. An old unknown pair stays unknown. Direct hashes of both frozen source/data packages are checked before reuse. In total3,227,336 pair tests are reused, including3,188,803 exclusions. Only396,129 previously unseen pairs require new exact product tests:396,080 cross-cohort pairs and49 additional pairs involving the newly represented large class.

The new tests exclude390,562 pairs. Their first applicable reasons are:

| Reason | New exclusions |
|---|---:|
|Distinct full stars sharing a quadrilateral pair|42437|
|Positive product of distinct walls|74971|
|Product in the open gap(−1,0)|116664|
|Nonintegral normalized two-cycle|41426|
|Squared native product outside K|115064|

Across the whole union there are3,623,465 pairs,3,579,365 excluded and44,100 unknown. Among cross-cohort pairs,390,520 are excluded and5,560 remain unknown.

**Certificate format:** each forest's `incompatibility_certificates` contains only newly computed exclusions. Its `inherited_pair_domains` gives the bound predecessor and a complete new-to-old star-index map. A consumer must expand those inherited relations as well as the new list before checking a capacity or reusing this output. The new list alone is not the complete incompatibility relation; the source and saved schema state this explicitly. This avoids duplicating millions of frozen records.

## Capacities and a concrete new obstruction

Within each forest, mutually incompatible orientation groups have capacity1. Universally incompatible groups give further capacity-one cliques. Proper colorings of the conservative literal compatibility graph bound the size of simultaneous actual collections. Every color class is checked pairwise incompatible; no optimal-coloring or realizability assertion is needed.

The263 new columns comprise131 clique columns and132 coloring columns. Their costs are131 ones,115 twos,16 threes and one four. Exact tuples of vertex, cost and all face coefficients are compared with both predecessor exports;1,263 duplicate candidates are omitted, with their current loads directly checked against the baseline allowance. Of the new columns251 involve both cohorts. All46 violations and all33 required allowance increases come from those mixed columns; omitting them leaves the previous ratio exactly unchanged.

The strongest new increase is at vertex4008, with forest A1+E6+H4+I2(5). Raw88 face45685 has two literal orientations and raw297 face68119 has four. Each orbit already has capacity1. All eight cross-orientation pairs are now excluded: four force a positive wall product and four force the nonintegral two-cycle

    (80+32√5)/9.

Its field trace160/9 is not an integer, so it cannot be an algebraic integer. For the positive-product pairs, the projection p is positive and

    p^2−h = 7/3 + 2√5/5 >0,

so p−√h>0 exactly. Consequently the two face occurrences share capacity1. Their total load is18134463/5000000000, while the old allowance is283351/156250000. This literal incompatibility explains the largest new repair, without relying on an abstract port-packing relaxation.

## Repair, preservation and limits

At each vertex the allowance is replaced by the maximum of its previous value and every new load divided by its positive cost. All4,410 face masses remain unchanged. The objective numerator remains122903790828233/9900000000000, and the new normalizer is10000790308037497/9828000000000000. Their quotient gives the displayed ratio12.200030329472785.

This preserves all436,262 columns of the frozen finite model and all1,002 columns of the immediately preceding repair. The96,176 zero-budget signed transfers depend only on unchanged face masses. The52 signed positive-cost columns impose load≤cost·mu, since their variables are nonnegative; increasing mu preserves them. No old fixed-mass simplex or joint-port covers are silently reused. Complete implicit coverage at the current masses, other actual-corner families and broader field/domain coverage remain unfinished.

## Verification and review

FINAL/FROZEN: final generation70380 and separate fresh saved replay13362 both terminate exit0/PASS. They reconstruct the current represented domain and every newly required product, while checking the exact frozen pair-domain transports. The source/data predecessor hashes remain explicit. A separate no-producer consistency walk28027 also terminates exit0/PASS: it expands all inherited exclusions, checks1,056 occurrence records and all263 clique/coloring capacities, and recomputes the33 allowance increases and exact ratio. It treats saved pair-geometry certificates as inputs and does not independently repeat their projection/embedding census. All direct source hashes, source syntax and local report links pass.

Construction's narrow full source/delta review passes the current explicit/unique lifting, variable core orders, complete-key transport, fixed-field embedding, conservative unknown-pair retention, capacity coefficients and exact old-column deduplication. This is a mathematical/source review, not a duplicate pair census or physical reconstruction. Any later review is recorded in the successor aggregate rather than changing this frozen report. No predecessor is edited and no optimizer is run.

    python3 scripts/fresh_mixed_raw87_large_core_corners_20261017.py
    python3 scripts/fresh_mixed_raw87_large_core_corners_20261017.py --verify
