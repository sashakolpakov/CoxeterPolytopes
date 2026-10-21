# Exact actual-corner capacities for raw178, raw296 and raw297

The new larger-core cohort removes the current fixed dual obstruction above12. Its469 necessary capacity columns include91 violated inequalities; raising65 vertex allowances with all4,647 face masses unchanged lowers the exact ratio from12.457543530716451 to **2877693523508930/247081190427751 =11.646752707185115**. This is a monotone repair of one witness, not an optimized feasible rule or a dimension theorem.

The [source](../scripts/fresh_raw178_296_297_actual_corners_20261017.py) and [certificate](fresh_raw178_296_297_actual_corners_20261017.json) use the new [mixed-raw24 finite-model dual](fresh_mixed_raw24_full_model_20261017.md), whose [support diagnostic](fresh_mixed_raw24_full_model_20261017_support.md) identified these47 uncovered explicit/unique classes. No optimizer or unchanged large packing scan is run.

## Complete represented cohort

| Raw core | Canonical core | Core order | Spectator rank | Classes | Actual corners | Complete literal orientations | Finite forests |
|---:|---:|---:|---:|---:|---:|---:|---:|
|178|158|9|6|18|144|770|124|
|296|267|10|5|16|128|6256|120|
|297|268|10|5|13|104|385|91|

There are41 already explicit provenance rows and6 rows with exactly one saved source context. Every positive row of these three families with explicit or uniquely saved provenance is included; the reconstructed set is compared with the frozen support diagnostic. No ambiguously merged row of these families occurs in this support, and no arbitrary representative is chosen.

All three selected aggregates are canonical four-dimensional dotted prisms, with8 actual vertices and3 occupied quadrilateral port pairs at each vertex. Their whole-face corner multiplicities are reconstructed from the frozen [full-star quadrilateral aggregation](fresh_quad_product_aggregation_20261017.md). Raw178 has9 roots and6 orthogonal spectator roots; raw296/297 have10 and5. Deleting each corner's two omitted roots leaves exactly13 incident walls after adjoining the spectators. The implementation uses these actual core orders, rather than the seven-root assumption of earlier families.

Every retained raw-core component is mapped into a whole matching component of the complete finite vertex forest, through every weighted graph isomorphism and every choice among repeated component copies. Different raw components use different target components. The code rechecks every labelled edge/nonedge and that the image consists of whole target components. Automorphisms acting only on unused spectator components have zero effect on the external covectors and occupied pairs and may be quotiented. The inherited rank-independent embedding and positive-definite inverse helpers are bound as source inputs.

Full literal stars are keyed by their unordered pair of external contact covectors in this common13-wall frame. Identical keys are merged before pair counting; their occupied pairs and native omitted-root products must agree. Distinct complete orientation orbits are verified disjoint, and each orbit belongs to exactly one Q-face class. Global column coefficients count actual prism-corner occurrences, never the number of orientations.

The resulting domain has376 corner occurrences,7,411 distinct literal stars and236 finite forests. The three families share84 forests; there are11,670 cross-family literal pairs.

## Exact geometry in one finite-vertex frame

All finite contacts have labels2,3 or5. The dotted squared products and their positive square roots are

    tau178=tau297=6+2 sqrt(5),       sqrt(tau)=1+sqrt(5),
    tau296=(7+3 sqrt(5))/2,         sqrt(tau)=(3+sqrt(5))/2.

Thus every external contact covector lies in K=Q(sqrt5); no independent contact radical is needed. Passive B/F/other components of the ambient finite forest are orthogonal to these covectors. Every component carrying a contact is checked to have a positive-definite normalized Gram with labels3/5, and its exact K inverse is reconstructed.

Let V be the complete positive-definite incident Gram, normalized to diagonal2, and let p be an external wall's full contact covector. Its squared timelike height is

    h_p=p^T V^(-1)p−2 >0.

Strict positivity is checked for every literal wall. The common finite vertex fixes the future sign of every genuinely external wall. Hence two actual external normals with covectors p,q have the unique forced product

    g=p^T V^(-1)q−sqrt(h_p h_q).

Every reconstructed star first reproduces its own prescribed omitted-root product exactly with this sign. Identical external covectors are the same actual outward wall and are skipped in pair exclusions.

Write a=p^T V^(-1)q in K and h=h_p h_q>0. The sign of a−sqrt(h) is negative when a<=0 and otherwise is the sign of a²−h. Before squaring, these exact signs exclude positive native products and the open Coxeter gap(−1,0).

For the additional fixed-ground-field gate,

    (a−sqrt(h))²=a²+h−2a sqrt(h).

If h is not a K-square, this belongs to K precisely when a=0, giving h. Otherwise the positive physical K square root gives the product directly. The code applies algebraic-integrality and finite-angle checks only after this membership test. Integrality of each normalized two-cycle and membership in the fixed ambient ground field K are necessary arithmetic conditions; this is explicitly a compact-arithmetic Q(sqrt5) model. The finite squared-product alphabet used for −2<g<0 is that of labels2,3,4,5,6,10. No numerical root or sign decision is used.

A separate combinatorial exclusion applies when distinct complete stars occupy the same two-port pair: that pair determines a unique actual quadrilateral, whose full wall star is unique. They cannot both occur at one actual vertex. These arguments are the same actual-normal/unique-face principles as the [first dominant-corner certificate](fresh_dominant_cap_compatibility_20261017.md), applied to the larger complete cores.

## Exhaustive pair checks and conservative capacities

All761,872 unordered literal pairs within their common finite forest are tested. Exactly760,379 are excluded, with the first applicable certificate classified as follows:

| Exclusion | Pairs |
|---|---:|
|Distinct full stars sharing a quadrilateral port pair|101546|
|Positive product of distinct actual walls|220673|
|Actual-wall product in(−1,0)|286772|
|Nonintegral normalized two-cycle|144468|
|Normalized two-cycle outside K|6920|

The1,493 unexcluded pairs remain compatible. Among cross-family pairs,10,683 are excluded and987 remain compatible. Lack of an exclusion is never treated as a realizability proof.

An orientation orbit whose distinct literal stars are all incompatible has capacity1. Universally incompatible such orbits form capacity1 cliques. More generally, in the conservative compatibility graph every simultaneous actual collection is a clique. A saved proper coloring therefore bounds its size by the number of colors. Every color class is checked pairwise incompatible. The code exports both individual-orbit colorings and whole-forest union colorings, together with capacity1 orbit cliques; it does not claim an optimal coloring or complete compatibility classification.

The469 distinct exported coefficient vectors consist of268 orbit-clique columns and201 coloring columns. Their costs are316 ones and153 twos. Of101 mixed-family columns,39 violate the initial witness; mixed columns require26 additional vertex increases beyond the within-family repair. Within-family columns alone already lower this witness to

    1438846761754465/121170865226688 =11.874527420949352.

The largest new loads occur in raw296 corners. For example, at vertex1385 with forest2A1+2A2+A3+A4, the two actual face classes68056 and68064 cannot coexist within the selected orbit family. Their capacity1 load is54310239/10^10, versus the old allowance7524359/10^10. This illustrates a geometric incompatibility absent from the previous abstract local-simplex and selected smaller-core constraints.

## Exact witness repair and remaining scope

For each affected vertex v set

    mu'_v=max(mu_v, max_new_columns_at_v load/vertex_cost).

All coefficients, costs and face masses are rational. The source directly checks every new column at the resulting mu' and verifies mu'>=mu coordinatewise. Exactly65 allowances rise; the support remains807 vertices. The objective numerator is unchanged at41109907478699/3300000000000 and the new normalizer is247081190427751/231000000000000, giving the ratio stated above.

This monotone change preserves every435,793 exported column of the frozen parent model, including its96,176 zero-budget signed transfers and52 signed columns with positive vertex cost. The signed columns are nonnegative charge variables with some negative face coefficients, not free primal variables. The [column constructor](../scripts/fresh_actual_corner_enlarged_model_20261017.py) adds them with positive vertex cost; the [inherited solver/verifier](../scripts/fresh_sqrt5_labelled_packing_full_model_20261017.py) applies bounds(0,None) to every active variable and checks load_j<=cost_j*mu_v. These are inequalities, not free-variable equalities. Increasing the vertex allowance therefore only relaxes them; zero-budget constraints depend only on the unchanged face masses. No independent rebuild of that entire matrix is claimed here: its exact generation and fresh replay are bound inputs.

The full implicit local-simplex packing and joint families have not been reconstructed at the current4,647 face masses. Previous complete covers at older masses are not reused. The present result proves that these new necessary actual-corner capacities defeat the saved finite-model witness; it does not settle a strengthened optimum, complete face/forest/field coverage or a dimension bound. Further optimization, if undertaken, must retain the inherited obligations and all new exported capacities.

## Verification scope

FINAL/FROZEN: final generation80828 and separate fresh saved replay54932 both terminated exit0/PASS. The initial generation63803/replay62818 also passed; the final runs additionally bind the nonnegative-variable engine and signed-column constructor directly. Construction_audit independently confirmed only the simpler fixed-K squared-product algebra; that statement is not a source or census review. The producer exactly reconstructs the full selected provenance domain, projections, native self-products, all literal pair tests, exported capacities and rational repair. The root full mathematical/source review passes, including the core orders, spectator ranks, exact K-height gates, actual provenance/orbits, multiplicity invariants, capacity extraction and the nonnegative-variable signed-inequality convention bound directly above. This is a written/source review, not an independent projection/pair census or full-matrix rebuild. Any later broader review belongs in a successor aggregate. No frozen predecessor was changed.
