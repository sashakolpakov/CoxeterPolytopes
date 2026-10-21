# Mixed raw39/raw23/raw25 actual-corner capacities

The new raw39 mixed constraints remove the current fixed witness above12. The exact repair changes its ratio from **12.382070564084598** to **11.585825626632172**. Without columns combining raw39 with another family, the same reconstruction only lowers it to **12.379327350872302**. All4,661 face masses remain fixed; only86 vertex allowances increase. This is a fixed-witness repair, not an optimized value, a feasible target12 rule, or a dimension theorem.

The [source](../scripts/fresh_mixed_raw39_actual_corners_20261017.py) and [certificate](fresh_mixed_raw39_actual_corners_20261017.json) start from the frozen [raw39 repair](fresh_raw39_actual_corner_compatibility_20261017.md). They reconstruct every positive raw23, raw25 and raw39 class whose provenance is either explicitly mapped in the current model or uniquely specified by its entire saved source-context list. No arbitrary representative is selected from a merged row. There are78 explicit classes and70 unique-context classes, with no ambiguous positive row from these three families omitted.

| Raw core | Prism dimension | Actual face classes | Actual finite corners | Vertex forests | Literal orientations |
|---|---:|---:|---:|---:|---:|
|23|5|63|630|340|7,094|
|25|4|54|432|307|5,892|
|39|2|31|124|88|1,246|

There are148 face classes,1,186 actual corner occurrences,409 whole rank13 vertex forests and14,232 literal orientations. Raw39 shares87 forests with at least one other family. Every corner retains five roots from a seven-root core and the complete rank8 spectator forest. All active components are embedded into whole matching vertex components, over every identical component copy and internal weighted automorphism. The unused components are orthogonal to both external covectors, so their permutations add no new literal stars.

Each literal star is identified by its unordered pair of complete external contact covectors in the same labelled vertex frame. Identical keys must agree in occupied port pairs and native omitted product. Distinct automorphism orbits are checked disjoint, and each orbit belongs to exactly one actual face class. Coefficients count actual prism-corner occurrences at vertices, not orientations or presentations. The raw39 selected classes and their source contexts are also checked against the frozen within-family package.

## Common-frame product and exact arithmetic

At a finite vertex, its incident norm2 Gram matrix G is positive definite. An external actual wall with complete covector p has positive timelike height

\[
h_p=p^TG^{-1}p-2>0.
\]

The future outward sign is fixed by the common vertex, so two such wall normals have forced product

\[
g(p,q)=p^TG^{-1}q-\sqrt{h_ph_q}.
\]

The source checks each star's prescribed native product with this formula. Identical covectors denote one wall and are skipped in distinct-wall tests. Different actual walls cannot have positive product or a norm2 product in the open interval(−1,0). Two distinct complete aggregate stars cannot occupy the same two-port pair: that pair determines a unique actual quadrilateral and its full wall star. Unknown pairs are retained as compatible.

All arithmetic exclusions use the fixed ambient ground field K=Q(sqrt5). Write

\[
\tau=(7+\sqrt5)/2,\qquad\sigma=(41+17\sqrt5)/2.
\]

The dotted raw25 square is6+2sqrt5=(1+sqrt5)². The norms of tau and sigma are11 and59; their product has norm649. Thus tau, sigma and tau·sigma are nonsquares in K and the common projection field is the biquadratic extension

\[
F=K(\sqrt\tau,\sqrt\sigma),
\qquad p=a+b\sqrt\tau+c\sqrt\sigma+d\sqrt{\tau\sigma}.
\]

Multiplication is exact in this four-term basis. Physical signs are computed by two nested quadratic comparisons. For g=p−sqrt(h), with h>0 in K, if p≤0 then g<0; otherwise sign(g)=sign(p²−h). No numerical field or sign decision is used.

There are exactly four ways sqrt(h) can lie in F: h is a K-square times one of1,tau,sigma,tau·sigma. The source tests all four using the complete exact K-square-root routine, chooses the physically positive root, subtracts it from p and squares in the four-term basis. The native two-cycle g² lies in K exactly when its three nonconstant coefficients vanish. If none of the four tests succeeds, sqrt(h) is outside F. In that quadratic extension, the coefficient−2p of sqrt(h) shows g² can lie in F, hence in K, only when p=0; its square is then h. These cases exhaust all cancellations.

A surviving two-cycle must be an algebraic integer in O_K. Finite acute products are additionally checked against the complete fixed-K squared alphabet0,1,2,(3+sqrt5)/2,3,(5+sqrt5)/2. These are necessary exclusions; no converse realization claim is made. The algebraic four-squareclass argument received an independent conceptual check by construction_audit, without an implementation or census replay.

## Certified capacities and repair

All627,684 unordered literal-star pairs are tested;608,207 are excluded and19,477 remain conservatively compatible. Across raw39 and the other two families specifically, all92,676 pairs are tested:90,792 are excluded and1,884 remain allowed. The latter exclusions have first reasons10,072 shared quadrilateral port pairs,24,746 positive wall products,13,236 products in(−1,0),30,496 nonintegral two-cycles and12,242 squares outside K.

There are784 internally incompatible corner orbits and816 universal incompatibility edges between such orbits. Maximal cliques of these orbits have capacity1. For each single orbit and each whole forest union, the source also verifies proper colorings of the conservative compatibility graph. Simultaneous actual stars form a clique in that graph, so the number of colors bounds their count. This remains valid when the coloring is not optimal and when some unknown pairs cannot actually coexist.

After duplicate coefficient tuples are removed,1,440 necessary columns remain. Of these,100 violate the input witness and require86 vertex increases. There are168 columns combining raw39 with raw23 or raw25;74 violate, and65 forests require an additional increase beyond all columns that do not mix raw39. For example, at forest3A1+A3+A4+H3 (vertex944), the capacity-one column is

\[
2\alpha_{26328}+\alpha_{26394}+\alpha_{26410}
+\alpha_{28678}+\alpha_{28760}\le\mu_{944}.
\]

Its exact load is223214567/10^10. The first three face classes are raw23 and the last two raw39. This is an actual corner-occurrence inequality; the coefficient2 represents two distinct corners of the same prism at distinct vertices in the global count.

The repair takes the maximum of the old allowance and every new load divided by its positive vertex cost. It preserves all face masses and directly verifies all1,440 new inequalities. Because only vertex allowances increase, all432,021 inherited model columns, including the52 signed columns with positive vertex cost, remain satisfied. The preceding complete local simplex packing covers, all17 selected joint families, cap constraints and incidence inequalities also remain satisfied by the same monotonicity. Their previously stated coverage limitations are unchanged.

The exact repaired ratio is

\[
\frac{247898037929594052000000000000000}
{21396665711915633070810190213789}
=11.585825626632172\ldots.
\]

Excluding the genuinely raw39 mixed columns gives instead

\[
\frac{247898037929594052000000000000000}
{20025162183963579836128012769579}
=12.379327350872302\ldots>12.
\]

No optimizer is invoked. The result removes this fixed obstruction through new mixed geometry; optimization over an enlarged finite domain and complete global face/field coverage are separate questions.

## Verification status

Final generation49584 and fresh saved replay53813 both terminated exit0 with `PASS_EXACT_MIXED_RAW39_RAW23_RAW25_FINITE_VERTEX_CORNER_CAPACITIES`. The replay reconstructs the cohort, full contact frames, exact projections, pair certificates, clique/coloring capacities and rational repair. Direct source hashes, Python syntax and local report links pass. The frozen predecessors are unchanged.

Root independently read the complete producer source and mathematical report and found no gap in the four-squareclass and nested-sign arguments, common-height formula, full provenance and active embeddings, orbit disjointness and occurrence coefficients, conservative capacities or monotone repair. This is a mathematical/source review, not an independent census replay. Construction_audit separately checked the new field algebra; that narrower check makes no source or census claim.
