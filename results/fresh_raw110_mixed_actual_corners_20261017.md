# Raw110 mixed actual-corner capacities at the weighted witness

Five uniquely represented raw110 classes supply new exact geometric constraints. Together with the currently positive classes of the prior geometry domain, they lower the full fixed witness from **12.133893289851665** to **11.920281779865938 < 12**. Omitting the new cross-cohort columns would leave **12.124167626625663 > 12**. Thus the improvement depends on mixed actual-corner incompatibility.

The [source](../scripts/fresh_raw110_mixed_actual_corners_20261017.py) and [certificate](fresh_raw110_mixed_actual_corners_20261017.json) preserve every frozen predecessor. The baseline is the full repair from the [complete weighted graph-cone package](fresh_weighted_actual_corner_capacities_20261017.md), at the 4,713 face masses of the [current finite model](fresh_raw33_295_mixed_full_model_20261017.md). No optimizer is run, and no new dimension or global feasibility claim follows from this fixed-witness repair.

## Complete selected domain

All 184 positive Q classes belonging to the frozen 243-class geometry export are retained. The new classes are exactly faces 48161, 48241, 48261, 48279 and 48317, each with one unique saved source context and canonical prism core90 from raw source110. No ambiguous merged row receives an arbitrary provenance.

Raw110 has core order8, prism dimension3 and six distinct actual corners. Its spectator rank is7, so each complete finite corner frame has rank13. The five classes give 30 actual corner occurrences and 650 distinct orientation occurrences over29 forests. There are24 common forests with the old cohort and34,708 cross-cohort literal pairs.

| Complete combined domain | Count |
|---|---:|
|Represented Q classes|189|
|Explicit provenance classes|130|
|Unique saved-context classes|59|
|Actual corner occurrences|1612|
|Literal stars|30109|
|Finite forests|491|

All six raw110 corners are reconstructed, with complete active-component embeddings into whole matching finite components and every choice among identical spectator copies. The previously retained corner/provenance records agree literally with their frozen versions. Whole-face coefficients count actual corner occurrences, not embeddings or orientations. Equal complete external-wall pairs are deduplicated, and distinct complete orbit sets are checked disjoint.

## The genuinely new square class

In the common norm-two convention, raw110 has dotted square

    tau110 = 15 + 6 sqrt(5) = 3 tau88.

This is a new square class over K=Q(sqrt(5)): sqrt(3) is not in K. Its norm45 is not a rational square, so tau110 itself is not a square in K. Exact K square-root tests also verify independence from each of the existing representative square classes. This is not the earlier raw33 identity tau33=4tau88, which uses the coefficient2 in K.

Every external wall retains its native contact encoding. The field for a self or cross projection is selected from those native encodings and contains at most two independent radicals. Missing radicals are padded only by an independent square class. The positive K coefficient of each native contact is recovered exactly in that field. No coefficient vector from one field is reinterpreted in another.

At a finite vertex, let J be the inverse full incident Gram, and let p,q be complete covectors of two external norm-two walls. Their squared timelike heights are

    h_p = p^T J p - 2 > 0,
    h_q = q^T J q - 2 > 0.

The common future side fixes their actual product as

    p^T J q - sqrt(h_p h_q).

Every height is checked to lie in K and be strictly positive, and each literal prism's omitted-wall pair reproduces its prescribed native product with this sign. Identical covectors denote one actual wall and are skipped as an incompatibility test.

The inherited exact pair-field gate then excludes positive products, products strictly between minus one and zero, squared products outside the fixed ambient K, and nonintegral native two-cycles. Its square-class cases and physical signs are exact. The finite Coxeter alphabet remains an additional necessary test where applicable. Every pair not excluded remains compatible; there is no inference from numerical near-equality.

## Pair provenance and capacity certificates

The previous full geometry supplies complete mapped pair domains. The first predecessor containing both stars owns the outcome, including unknowns. Recursive transport checks the whole forest, complete native covectors, occupied square pairs, both squared heights and the native omitted product. A direct new exclusion is forbidden inside any inherited mapped domain, so it cannot override an inherited unknown.

| Pair result | Count |
|---|---:|
|Complete literal pairs|4254990|
|Inherited pairs reused|4182709|
|New pairs tested|72281|
|New exact exclusions|71939|
|Total exact exclusions|4170376|
|New cross-cohort pairs|34708|
|Cross-cohort exclusions|34620|
|Cross-cohort unknowns retained|88|

The new exclusion reasons are10,457 shared occupied-quadrilateral pairs,14,804 positive wall products,6,806 products in(-1,0),20,910 nonintegral two-cycles and18,962 two-cycles outside K. These total71,939. A shared occupied pair determines the same actual quadrilateral, so two distinct complete aggregate stars cannot both occupy it.

Internally incompatible complete orbits have capacity one; universally cross-incompatible orbit cliques retain capacity one. More general simultaneous sets form cliques in the conservative literal compatibility graph, so explicitly checked proper colorings give safe integer capacities. The code retains every unknown edge. It exports263 new columns after1,838 exact historical tuple duplicates, including the weighted predecessor's exported columns, are removed. Their coefficients preserve actual corner multiplicities.

Of the263 columns,45 mix raw110 with the old cohort. There are21 violated columns, of which17 are mixed, and14 vertex allowances increase. Eleven allowances increase further because of the mixed columns. These are the stated clique/coloring exports, not a claim to have exhausted all weighted capacities in the enlarged raw110 graph.

For example, forest2902 is `A1+A2+H3+H3+H4`. One new capacity-one clique proves

    alpha_26448 + 2 alpha_26460 + alpha_45769 + alpha_48317 <= mu_2902.

It mixes raw23, raw88 and raw110. Its exact required allowance exceeds the weighted baseline by40746569/10000000000, the full mass of the leading new raw110 class. The same increase is forced at forest3324, `A1+A5+H3+H4`, by another saved mixed clique. The certificate records every literal star group used in each inequality.

## Monotone full-witness repair

All face masses remain fixed. At each forest the allowance is raised to the maximum of its old value and every new load divided by its certified capacity. The objective stays12258439268597/990000000000. The final normalizer is112297812488946653/108108000000000000, giving

    1338621568130792400 / 112297812488946653
      = 11.920281779865938.

The complete 439,061-column inherited finite model remains satisfied. Its96,176 zero-budget signed columns depend only on the unchanged alpha. All ordinary and52 positive-cost signed columns impose load <= positive_cost*mu; increasing mu preserves them. The full old weighted orbit-capacity cone is also preserved: its exact mixtures were constructed at these same face masses, and their required allowances only become easier to satisfy.

No old-alpha implicit simplex cover or joint-port scan is reused. The new raw110 capacities are necessary local geometric constraints in the selected fixed-K domain. The repaired ratio is not a new optimized finite-model value, a primal target12 certificate, a realizable hyperbolic polytope or an ambient dimension bound.

## Verification

Domain generation95134, full generation37543 and separate fresh saved replay90771 all terminated exit0/PASS. The separate saved-capacity checker39386 also terminated exit0/PASS: it imports no geometry producer, recursively reconstructs the first-owned saved relations, checks all actual occurrence coefficients, clique/coloring capacities, duplicate handling, exact loads and both with/without-mixed repairs. It does not independently reconstruct native projections or the physical orientation census.

Root completed a full new-source/delta and mathematical report read, including the inherited capacity ownership and deduplication section. It found no gap in the weighted baseline at the same4713 face masses, retention of all184 old positive classes, the new tau110 square-class guards, first-owner unknown preservation, historical tuple load checks or monotone preservation of the complete old weighted cone. This was a source and mathematical review, not an independent physical pair, covector or orientation census.

All50 direct source bindings, syntax and local report links pass. Source, data and report are final/frozen. No optimizer was run for this package, and all frozen predecessors remain unchanged.
