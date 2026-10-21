# Mixed raw24/raw23/raw25/raw39 actual-corner capacities

The new mixed capacities invalidate the current fixed witness above12. Keeping all5,059 face masses fixed and increasing only vertex allowances lowers the [within-raw24 witness](fresh_raw24_actual_corner_compatibility_20261017.md) from12.417934074952033 to **11.06349648741785**. If every genuinely raw24-mixed column is omitted, the corresponding repair remains12.417926597612544. Thus the crossing comes from the new mixed compatibility constraints, rather than just recertifying older families at the current masses.

This is an exact fixed-witness repair. It is not a reoptimized finite-model value, a feasible target12 rule or a dimension theorem. No optimizer is invoked.

The [source](../scripts/fresh_mixed_raw24_actual_corners_20261017.py) and [certificate](fresh_mixed_raw24_actual_corners_20261017.json) reconstruct the positive cohorts from the **new optimized** face-mass vector in the frozen [mixed-raw39 model](fresh_mixed_raw39_full_model_20261017.md). They do not reuse the earlier mixed package's support.

| Raw family | Canonical core | Positive classes | Actual corners | Whole forests | Literal orientations |
|---|---:|---:|---:|---:|---:|
|23|19|71|710|395|7634|
|24|20|45|360|219|4522|
|25|21|54|432|305|5414|
|39|29|11|44|33|318|

In total181 classes contribute1,546 actual corner occurrences,17,888 literal orientations and458 whole vertex forests. Ninety-three classes use the inherited explicit actual-provenance map. The remaining88 have exactly one saved source context each, checked against the complete context list before lifting. No ambiguous class is represented arbitrarily. The raw24 selection agrees exactly with the frozen within-family45-class census. The changed counts for raw23 and raw39 reflect the new optimized support, not deletion of positive records during this calculation.

## Common finite-vertex geometry and actual multiplicities

For each selected prism, the source reconstructs every finite corner from its canonical quadrilateral partition. Its rank8 orthogonal spectators complete the retained five-root corner to a rank13 finite vertex. Each whole active finite component embeds into every matching whole vertex component, including all identical copies and all internal weighted automorphisms. Inactive components have zero covectors and no occupied ports and can be quotiented. Distinct orientations are keyed by the unordered pair of **complete** external contact covectors in one labelled full vertex frame. Repeated keys must agree in native omitted product and occupied square-port pairs. The complete corner orbits are disjoint, and one orbit belongs to one actual face class.

If G is the positive definite norm2 incident Gram matrix and p is a complete external covector, its positive timelike height is h_p=p^T G^-1 p−2. The common future outward sign fixes the product of two actual external walls:

    g(p,q)=p^T G^-1 q-sqrt(h_p*h_q).

Every height is checked strictly positive and every selected star's prescribed native omitted product is recovered exactly. Identical covectors are the same actual wall and are skipped before an incompatibility test. A positive product or a norm2 product in(−1,0) excludes two distinct walls. Distinct full stars sharing an occupied two-port pair are also incompatible: the pair determines a unique actual quadrilateral and therefore its unique full wall star. Coefficients in the final inequalities count actual prism-corner occurrences over the whole face, not orientation counts.

## Exact arithmetic uses one field per wall pair

The fixed ambient ground field is K=Q(sqrt5). The dotted squares are

    tau23=(7+sqrt5)/2,       Norm(tau23)=11,
    tau24=(23+9*sqrt5)/2,    Norm(tau24)=31,
    tau25=6+2*sqrt5=(1+sqrt5)^2,
    tau39=(41+17*sqrt5)/2,   Norm(tau39)=59.

The three nonsquare parameters have pairwise independent K-squareclasses: their norms and pairwise norm ratios are nonsquares in Q. One external-wall pair contains at most two of these radicals. The source therefore constructs a separate biquadratic field F=K(sqrt(tau),sqrt(sigma)) for that pair. If fewer than two radicals occur, it pads with independent unused parameters; this changes no physical value or exact membership decision. Every saved expression includes its two bases. Cached field operations belong to the corresponding field object; no coefficients or cache entries are transferred between differently based fields.

Multiplication uses the four-element basis1,sqrt(tau),sqrt(sigma),sqrt(tau*sigma). Nested exact sign comparisons select the positive physical square roots before squaring. The inverse finite-vertex projections and each individual squared height are checked to lie in K. For g=p−sqrt(h), with p in F and h>0 in K, the four K-squareclasses represented by F are exactly1,tau,sigma,tau*sigma. If h belongs to one, the source subtracts its positive root in that basis, squares g and requires all three nonconstant coefficients to vanish. If sqrt(h) lies outside F, then p≠0 implies g²=p²+h−2p*sqrt(h) lies outside F and hence outside K; p=0 leaves g²=h. This exhausts the cancellation cases, including harmless padded fields.

After fixed-K membership, the native two-cycle must be integral in O_K. A finite acute product must also satisfy the fixed-K Coxeter angle alphabet. Unknown or unexcluded pairs remain compatible; no sufficiency of the arithmetic tests is asserted. The pair-field algebra generalizes the frozen [mixed-raw39 engine](fresh_mixed_raw39_actual_corners_20261017.md), without requiring a global field containing all three radicals at once.

## Complete stated pair domain and capacity certificates

The source tests950,774 unordered literal-star pairs in the458 forests. It excludes933,174 and leaves17,600 conservatively compatible. Of540,810 cross-family pairs,528,980 are excluded. Raw24 has common forests with another selected family in212 cases; all339,222 such cross pairs are tested,335,420 are excluded and3,802 remain allowed.

The first recorded exclusion reasons are157,319 shared occupied square pairs,191,594 positive distinct-wall products,162,080 products in(−1,0),401,793 nonintegral native two-cycles and20,388 squares outside K. There are1,134 internally incompatible corner orbits and1,679 universal incompatibility edges between them. Maximal orbit cliques give capacity-one inequalities. Verified proper colorings of the conservative literal compatibility graphs give further capacities, including whole-forest union columns. Actual simultaneous stars form a clique in that graph, so each proper coloring is a valid upper capacity. Optimality of the coloring is unnecessary and is not claimed.

The1,793 exported columns include357 involving raw24 and another family. In the baseline witness148 columns violate their capacities;139 of those are genuinely raw24-mixed. Taking the simultaneous maximum allowance at each forest raises120 allowances,115 beyond the repair omitting raw24-mixed columns.

For example, the largest new capacity-one load at vertex1505 is

    alpha26331 + 2*alpha26392 + alpha26410
      + 2*alpha26627 + 2*alpha26645 <= mu1505.

Its actual-star orbit clique is saved with the literal pair exclusions and whole-face multiplicities. Every exported load is recomputed exactly before and after repair.

The unchanged objective is69217670225393/5544000000000. The repaired normalizer is127083144831498288000922893661/112612500000000000000000000000, giving

    1405983926453295312500000000000 / 127083144831498288000922893661
    =11.06349648741785... .

All old and new inequalities remain valid for this repaired dual certificate: face masses are fixed and every vertex allowance increases or stays unchanged. This preserves all433,753 already verified finite-model columns, including the52 signed columns with positive vertex cost, as well as the newly reconstructed complete represented simplex-packing covers and all17 selected joint families. No cover is reused at changed face masses.

## Verification and remaining scope

Primary generation98701 and fresh saved replay93854 both terminate exit0/PASS. The replay recomputes the full stated corner/orbit domain, exact pair gates, capacities and repair without an optimizer.

A separate narrow source/algebra read by construction_audit passes the PairField routines, nested sign rule, four squareclass test, independent padding and strict finite-angle test. That review does not claim a duplicate projection, orbit, provenance or complete pair census. Root's separate full source/mathematical read passes the current explicit/unique cohort, common-frame embeddings and projections, pair-field signs and native self-checks, wall identity, orbit coefficients, conservative capacity extraction and exact with/without-mixed repair; it does not duplicate the full pair census. The inherited geometric principles and older source reviews are not presented as a fresh independent census of this package.

Only the stated positive explicit/unique cohorts receive these new capacities. Other merged provenance classes, other local diagrams, missing face and vertex domains, and global chamber compatibility remain unresolved. The exact fixed-witness crossing justifies studying an enlarged model in a subsequent step; it does not supply its optimum or feasibility.
