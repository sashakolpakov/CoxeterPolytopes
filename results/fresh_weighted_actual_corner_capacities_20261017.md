# Complete weighted capacities in the represented actual-corner graphs

The complete nonnegative orbit-weight cone strengthens the saved actual-corner capacities. At the new 4,713 fixed face masses it requires **84 vertex-allowance increases**, lowering the full witness from 12.382248317539402 to **12.133893289851665 > 12**. The separate restriction to represented quadrilateral classes requires 86 increases and drops from 12.076081050570941 to **11.8338646495987 < 12**.

Every positive represented forest has an exact matching charge and covering-mixture certificate. Consequently, further nonnegative orbit weights on these same conservative compatibility graphs cannot eliminate the repaired full witness. This is a limitation of this graph domain at these fixed masses, not a limitation of new geometric compatibility, new face classes or other incidence constraints.

The [source](../scripts/fresh_weighted_actual_corner_capacities_20261017.py) and [certificate](fresh_weighted_actual_corner_capacities_20261017.json) bind the [243-class geometry package](fresh_raw33_295_mixed_actual_corners_20261017.md), the [new finite model](fresh_raw33_295_mixed_full_model_20261017.md), and its [support diagnostic](fresh_raw33_295_mixed_full_model_20261017_support.md). They do not repeat the physical projection or pair-field census.

## The actual graph domain

The current masses have 184 positive Q classes among the 243 exported class identifiers. Their complete saved corner orbits occupy 486 finite forests, with 1,582 orbits and 29,459 distinct literal stars. No forest has more than 14 positive orbits.

The frozen geometry supplies full common-frame covectors, occupied pairs, squared heights, native omitted-wall products and actual whole-face corner multiplicities. Inherited pair relations are expanded recursively. The first mapped predecessor containing both stars owns the pair, including an unknown outcome. Transport checks all these literal data, and a new direct exclusion is additionally required to lie outside every inherited mapped pair domain. Thus no later exclusion overrides an inherited unknown.

The compatibility graph joins every pair not excluded by the saved exact geometry. Unknown pairs remain edges. An actual simultaneous set of represented corner stars is therefore a clique in this graph. Each distinct literal star can occur at most once at the same vertex, as established in the bound geometry package. Positive-mass orbits retain their entire orientation sets; omitted classes receive zero new payment.

For orbit g, its demanded count is

    b_g = sum_face m_(g,face) alpha_face,

where m counts actual corner occurrences, not orientations. This is recomputed from the new alpha. The restricted witness has identical masses on all represented Q classes, so it has the same b_g. Its different initial vertex allowances are kept separate throughout.

## Exact capacity and the full nonnegative weight cone

For each forest, enumerate every maximal literal compatibility clique. Any actual clique extends to one, and all proposed orbit weights are nonnegative. The maximum weighted charge is therefore attained on a maximal clique. Replace each clique by its vector of counts in the complete orbits. Identical vectors may be merged; a coordinatewise dominated vector may be discarded. Every retained vector keeps a concrete compatible literal-clique representative.

For weights w_g >= 0, the normalized charge constraints are

    sum_g c_(j,g) w_g <= 1  for every retained count vector j.

The best charge against the current demands is the finite local problem

    maximize sum_g b_g w_g.

Its covering dual uses nonnegative mixture masses lambda_j with

    sum_j lambda_j c_(j,g) >= b_g  for every g,
    minimize sum_j lambda_j.

The saved charge and mixture satisfy both systems exactly and have equal objective values. Hence each local optimum is proved by rational arithmetic. Local numerical LPs only propose supports; exact linear systems recover both certificates. A saved replay enumerates the graph cliques and verifies the certificates without any optimizer.

The covering certificate also explains the universal statement directly. For any other nonnegative weights w, let h(w) be their maximum literal-clique charge. Then

    sum_g b_g w_g
      <= sum_j lambda_j sum_g c_(j,g) w_g
      <= (sum_j lambda_j) h(w).

Thus the mixture covers every nonnegative orbit-weight inequality, not only the particular optimizer output or subset-coloring cuts. The literal cliques are graph-compatible; they need not be simultaneously realizable by actual hyperbolic walls. This distinction makes the result a certified obstruction within the conservative graph relaxation.

## Complete results

| Exact domain or certificate count | Count |
|---|---:|
|Positive represented Q classes|184|
|Completed forests|486|
|Complete positive corner orbits|1582|
|Literal stars|29459|
|Maximal literal cliques enumerated|95989|
|Distinct nondominated orbit-count vectors|1477|
|Positive exact mixture terms|1299|
|Exported weighted columns|86|
|Full-witness vertex raises|84|
|Represented-Q restriction vertex raises|86|
|Unfinished forests|0|

The enumeration guards allow an unfinished forest to be recorded explicitly, but none remains. All 1,477 distinct vectors are nondominated in this run. The 86 exported columns are precisely the local optimal charges violated by at least one of the two baseline witnesses. Each has positive vertex cost one and nonnegative rational face coefficients. No unchanged global finite-model optimization was run.

For example, at forest 3137, `A1+A3+H3+H4+I2(5)`, the five orbit weights are `(1,1,1/2,1/2,1)`. Their actual occurrence coefficients give

    2 alpha_26474 + (alpha_45881 + alpha_26918)/2 + alpha_26448 <= mu_3137.

The exact demanded allowance is 411298683/20000000000, whereas the previous allowance was 188073201/10000000000. The saved five count vectors and literal representatives certify the inequality and an equal-mass covering mixture. This combines orbits with different simultaneous multiplicities, which the old single-orbit and all-orbit colorings did not fully exploit.

## Separate repaired witnesses and preservation

For each forest v let beta_v be its certified weighted-cone optimum. Set its new allowance to max(old_mu_v,beta_v), separately for the full and restricted witnesses. All face masses remain fixed.

| Witness | Exact repaired ratio | Decimal |
|---|---|---:|
|Full current finite-model dual|83663848008174525 / 6895053880039298|12.133893289851665|
|All S plus represented Q restriction|1305510807036181200 / 110319903572705863|11.8338646495987|

Every one of the 439,061 inherited exported model columns is preserved by this monotone change. The 96,176 zero-budget signed columns depend only on alpha. The 52 positive-cost signed columns, like all ordinary budgeted columns, impose load <= cost*mu; increasing mu preserves them. Their primal variables are nonnegative despite signed face coefficients. The already verified restricted alpha is also unchanged, with no simplex restoration.

No previous-mass implicit simplex-packing cover or joint-port scan is imported. The result does not prove target12 feasibility, a new optimal value for the enlarged global finite model, a globally realizable polytope or a dimension bound. In particular, the full repaired witness remains above12 even after all nonnegative orbit weights in this domain are admitted.

## Verification

FINAL/FROZEN: final generation 4169 and fresh separate-process saved replay 86895 terminate exit 0/PASS. Both reuse or verify only exactly matched charge/cover certificates while reconstructing the graph clique domain; neither invokes an optimizer. The generation that first produced these certificates used bounded local LP proposals, with numerical recovery failures resolved by objective normalization and exact rechecking. No forest remains open.

A separate checker importing no project producer passes all 486 charges by an independently written integer weighted-clique search, all 1,299 literal compatible mixture terms and their orbit marginals, both exact normalizers and all 86 occurrence-coefficient columns. The independent search visits 20,795 branches. It uses the bound saved geometric relations rather than reconstructing native projections or pair exclusions.

The parent read the complete source and report, finding no issue in the exact charge/cover argument, and requested the explicit direct-exclusion/inherited-unknown guard now present. Separately, the parent used no-producer Fraction arithmetic to check all 486 saved charges against all 1,477 count vectors, all 1,299 covering terms and demands, and both repaired normalizers and ratios. That check does not duplicate clique enumeration or physical pair geometry. Direct source bindings, report links and source syntax pass. All these are internal checks, not external mathematical refereeing.
