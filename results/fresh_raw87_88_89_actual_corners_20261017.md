# Joint actual-corner capacities for raw87, raw88 and raw89

The new joint domain yields1,002 necessary capacity columns, including113 violated inequalities at the current finite-model dual. Keeping all4,410 face masses fixed and raising75 vertex allowances repairs these inequalities and lowers the exact ratio from12.41451686839335 to

    1342109395844304360 / 109100031098664527
        =12.301640818329087.

The ratio remains above12. This is an exact repair of the saved witness, not an optimized bound for an enlarged model or a dimension theorem. No optimizer is invoked. The [source](../scripts/fresh_raw87_88_89_actual_corners_20261017.py) and [saved certificate](fresh_raw87_88_89_actual_corners_20261017.json) bind the [new large-core finite model](fresh_raw178_296_297_full_model_20261017.md) and its [support diagnostic](fresh_raw178_296_297_full_model_20261017_support.md).

## Represented domain and full vertex frames

| Raw core | Canonical core | Classes | Actual corners | Literal orientations | Finite forests |
|---:|---:|---:|---:|---:|---:|
|87|76|36|288|9530|239|
|88|77|27|216|9248|176|
|89|78|26|208|1042|173|

All three cores have8 roots and7 orthogonal spectator roots. The selected aggregates are four-dimensional canonical dotted prisms with8 actual vertices and3 occupied quadrilateral port pairs per vertex. Deleting the two omitted roots at each corner and adjoining the spectators gives exactly13 incident walls.

There are51 explicitly split provenance classes and38 classes with exactly one saved source context. Every current positive class of these three families with unambiguous provenance is included, and the reconstructed set is compared exactly with the frozen support diagnostic. No ambiguous row of these families occurs in this support. No merged row is assigned an arbitrary representative.

The712 actual corner occurrences give19,820 distinct literal stars across380 complete finite vertex forests. The families share154 forests, with435,936 mixed-family literal pairs. Every whole retained raw component is embedded into a whole matching component of the target forest, through all weighted internal automorphisms and all choices among identical spectator copies. Different retained components have disjoint target images. The code rechecks all edges, nonedges and complete-component boundaries. Only permutations with no effect on external covectors or occupied pairs are quotiented.

Full stars are keyed by the unordered pair of their complete external contact covectors. Identical keys are merged, with equality checks for occupied pairs and native omitted-root products. Complete orientation orbits are verified disjoint and each belongs to one Q-face class. Column coefficients count actual prism vertices, including repeated isomorphism types at distinct vertices; they do not count orientations as additional actual occurrences. These are the frozen [full-star quadrilateral aggregation](fresh_quad_product_aggregation_20261017.md) conventions.

## Exact field and common future signs

Write K=Q(√5). All finite contacts have labels2,3,5. The dotted squared products are

    tau87=(9+3√5)/2,
    tau88=5+2√5,
    tau89=6+2√5=(1+√5)^2.

The first two give independent square classes overK. For tau87 the norm is9, but the rational square-root equations require15/4 or3/4 to be a rational square; neither is. Tau88 has nonsquare rational norm5, and tau87/tau88 has nonsquare norm9/5. The source also checks these three nonsquare tests exactly. Thus every selected external covector lies in the common biquadratic field

    E=K(√tau87,√tau88)

with basis1,√tau87,√tau88,√(tau87 tau88). The raw89 contact uses its positive K root1+√5. The implementation reuses the frozen immutable `PairField` arithmetic with these new parameters; coefficients are never interpreted in a different field.

Let V be the full positive-definite incident Gram, normalized to diagonal2. Every component carrying a selected external contact has labels3/5 and its exact inverse is computed overK. Other finite components are passive: all selected covectors vanish there, so they contribute zero to every projection. Each external wall with contact vector p has squared timelike height

    h_p=pᵀV⁻¹p−2.

For every literal wall the source verifies that this quantity is inK and strictly positive. It does not assume that arbitrary E-valued covectors have K-valued squared height. At the common finite vertex, the outward side fixes the positive future coefficient, so two genuinely external walls with covectors p,q have the unique forced product

    g=pᵀV⁻¹q−√(h_p h_q).

Every reconstructed star reproduces its own prescribed omitted-root product with this sign. Identical contact covectors determine the same outward wall and are skipped in incompatibility tests; a wall's positive self-product is not treated as an exclusion.

For a=pᵀV⁻¹q inE and h=h_p h_q inK, exact nested radical comparisons determine the signs ofa−√h anda+1−√h before any squaring. They exclude positive products of distinct actual walls and the open Coxeter gap(−1,0).

For the arithmetic gate, √h lies inE exactly when the K-square class ofh is one of1,tau87,tau88,tau87 tau88. In that case the positive physical root is placed in the corresponding basis coordinate, the productg is squared exactly, and membership inK requires all three nonconstant coefficients to vanish. If √h is outsideE and a is nonzero, then

    g²=a²+h−2a√h

is outsideE and hence outsideK. If a=0 its square ish inK. These cases exhaust the field-membership test. Algebraic integrality of the normalized two-cycle is tested only after membership in the fixed ambient ground fieldK. For−2<g<0 the additional finite squared-product alphabet is that of labels2,3,4,5,6,10. These are necessary conditions for the compact-arithmetic fixed-K model; there is no claim for unrestricted ground fields.

There is also a combinatorial exclusion: two distinct complete stars cannot occupy the same quadrilateral port pair. That pair determines a unique actual two-face, whose complete wall star is unique. This uses the same actual-face principle as the [first dominant-corner certificate](fresh_dominant_cap_compatibility_20261017.md).

## Exhaustive pair tests and conservative capacities

All2,567,854 unordered literal pairs within a common finite forest are tested. Exactly2,530,570 are excluded;37,284 remain allowed. The first applicable exclusion certificates have counts

| Exclusion | Pairs |
|---|---:|
|Distinct complete stars sharing a quadrilateral port pair|389179|
|Positive product of distinct actual walls|741367|
|Product in(−1,0)|495111|
|Nonintegral normalized two-cycle|574369|
|Normalized two-cycle outsideK|330544|

Of435,936 mixed-family pairs,432,944 are excluded and2,992 remain allowed. Unknown pairs are retained in the conservative compatibility graph; they are not asserted realizable.

There are712 complete corner orbits, of which616 have pairwise incompatible distinct orientations and therefore capacity1. The graph of universally incompatible capacity1 orbits has383 edges. Its maximal cliques give440 capacity columns. Proper colorings of the conservative literal compatibility graph give another562 columns, covering individual orbits and complete represented forest unions. Every color class is checked pairwise incompatible. A simultaneous actual collection is a clique of this compatibility graph, so its size is bounded by the number of colors. No optimal-coloring or complete-realizability claim is needed.

The1,002 exported columns have costs776 ones,171 twos,39 threes,10 fours and6 fives. There are286 mixed-family columns,41 of which violate the initial witness. Mixed columns require21 additional vertex increases beyond the within-family repair. Using only the within-family columns gives the larger ratio

    1342109395844304360 / 108679674535598087
        =12.349221706628278.

A concrete mixed violation is column255 at vertex668, with forest4A1+D4+H3+I2(5). Its three selected raw87/raw88 face classes45429,45480,45676 have total capacity1. Their exact load is7627641/5000000000, versus the previous allowance3119741/5000000000. This forces the largest new vertex increase. The saved orbit and pair witnesses identify the actual contact configurations, rather than relying on abstract simplex-port disjointness alone.

## Repair and inherited obligations

For every vertex, the new allowance is the maximum of its previous allowance and every new column's load divided by its cost. All face masses remain unchanged. The source checks all1,002 inequalities at these rational allowances and verifies their coordinatewise domination of the old allowances. Seventy-five allowances rise; the positive support remains984 vertices. The unchanged objective numerator is122903790828233/9900000000000 and the new normalizer is109100031098664527/108108000000000.

This monotone change preserves every436,262 exported column of the frozen parent model. Its96,176 zero-budget signed transfers depend only on the unchanged face masses. Its52 signed columns with positive vertex cost, like all other positive-cost columns, impose load≤cost·mu. Their variables are nonnegative; signed refers to face coefficients, not free variables or dual equalities. The bound [engine](../scripts/fresh_sqrt5_labelled_packing_full_model_20261017.py) and [signed-column constructor](../scripts/fresh_actual_corner_enlarged_model_20261017.py) establish that convention. No duplicate complete matrix rebuild is claimed in this geometric package.

The complete implicit local-simplex and joint-port families have not been rebuilt at the current4,410 masses. Older complete covers are not reused. This package covers exactly the stated actual raw87/88/89 cohort; it does not test mixed pairs with every previously treated family or resolve other face, forest, provenance and ground-field gaps. Since the repaired witness remains above12, this result does not remove the current obstruction and does not warrant a dimension claim.

## Verification

FINAL/FROZEN: domain run60460, full generation57067 and separate fresh saved replay74134 all terminate exit0/PASS. Generation and replay each reconstruct every selected provenance, literal orientation, positive height, native self-product, pair exclusion, capacity and exact repair in about50seconds. All27 direct source hashes, Python syntax and8 local links pass.

Construction's narrow complete new-source/delta read passes the fixed two-radical setup, exact nonsquare guards,51/38 provenance domain, order8/spectator7 corner reconstruction, complete-component embeddings, native self-products and positive K heights, actual occurrence coefficients, conservative capacities and monotone signed-column preservation. That review reuses the previously reviewed `PairField` and helper mathematics; it does not independently rebuild orientations, rerun the2.568-million pair census or verify the numerical output. Any later broader review is to be recorded in a successor aggregate. No frozen predecessor is changed and no optimizer is called.
