# Raw-25 actual-corner compatibility on the new dual

The exact calculation reconstructs all **55 positive raw-25 classes** in the selected cohort, their **440 finite prism corners**, and **6,260 distinct literal orientations** in 296 complete rank-13 vertex forests. It produces **500 necessary incidence inequalities**, of which 28 violate the current fixed witness. Raising 25 vertex allowances gives

\[
 \frac{67685455440946980000000000000000}
 {5466425626097138727892456269673}
 =12.38203170968089\ldots,
\]

from the preceding 12.43082950275823. All face masses remain unchanged. The repaired witness is still above 12. This is an exact fixed-witness repair, not an optimized value, a feasible target-12 weighting, or a dimension theorem. No optimizer is used.

## Domain and provenance

The baseline is the [enlarged-model joint-family certificate](fresh_actual_corner_enlarged_model_20261017_joint_ports.md), which already includes the complete local simplex-packing covers and its stated 17 joint families at these face masses. The [support report](fresh_actual_corner_enlarged_model_20261017_support.md) identified 19 positive raw-25 rows in the existing explicit-provenance map and 36 additional positive rows with exactly one saved source context.

Every additional row is admitted only after verifying that its entire original `source_contexts` list has length one. Its context, spectator forest, and canonical core are checked against the frozen producer data. No representative is chosen from a row with several saved provenances. The selected existing rows retain their original explicit provenance indices; newly admitted rows retain the original face index and unique context, without inventing an index in the old map.

All 55 classes have canonical core 21, raw source 25, and prism dimension four. Their aggregate is a four-prism with six square faces and eight vertices. Its seven-root core has squared normalized dotted product

\[
 \tau=6+2\sqrt5=(1+\sqrt5)^2.
\]

Eight spectator ranks complete each finite vertex to rank 13. The program reconstructs each of the eight retained five-root corner diagrams and recognizes its whole finite components. The additional positive-definite two-deletion of this raw core is not a vertex of the selected maximal prism; the corner list is obtained from the six actual square partitions, as in the [full-star aggregation theorem](fresh_quad_product_aggregation_20261017.md).

Each retained whole component is embedded into every matching whole component of the complete vertex forest, including identical spectator copies and every internal weighted automorphism. Passive component permutations have zero external contacts and no occupied square pairs, so quotienting them loses no literal data. A literal star is keyed by its unordered pair of complete external-wall contact covectors. Duplicate keys are identified, and agreement of occupied square pairs and native external product is asserted. Distinct complete orientation orbits are disjoint. Each such orbit has only one Q-face class among its origins.

The global coefficient of an orbit counts actual prism-corner occurrences, including different vertices of one prism. It does not count the number of automorphism presentations or literal orientations of one corner.

## Exact pair exclusions

Use normalized squared wall norm 2. At a finite vertex let \(G\) be the positive-definite Gram matrix of its 13 incident walls, and let \(v\) be its future unit timelike vector. For an external wall with complete incident contact covector \(p\), its spatial projection is \(s_p\), where

\[
 s_p^2=p^TG^{-1}p,\qquad
 r_p=s_p+\sqrt{p^TG^{-1}p-2}\,v.
\]

The square root is strictly positive because this is a genuinely external wall. Its sign is fixed by the common inward half-space convention: \(r_p\cdot v<0\). Consequently two distinct external walls have the forced product

\[
 r_p\cdot r_q=p^TG^{-1}q-
 \sqrt{(p^TG^{-1}p-2)(q^TG^{-1}q-2)}.
\]

This is the same finite-vertex argument used in the [raw-23 family](fresh_dominant_cap_compatibility_family_20261017.md). Here every active corner component has labels 3 and 5. Other spectator components are orthogonal to all external covectors. Their normalized Gram entries therefore need not be placed in the active \(K=\mathbb Q(\sqrt5)\) calculation.

Two different literal stars are excluded if they occupy the same two-port pair: the corresponding actual quadrilateral determines its full wall star. Otherwise all four external-wall pairs are examined. Identical complete covectors denote the same actual outward wall and are skipped. A positive forced product contradicts acuteness; a product in \((-1,0)\) is not a Coxeter product of distinct norm-2 walls.

There are also fixed-ground-field arithmetic gates. Because \(\sqrt\tau=1+\sqrt5\), every forced product has the form \(g=p-\sqrt h\), with \(p,h\in K\) and \(h>0\). If \(p\ne0\) and \(h\) is not a square in \(K\), then

\[
 g^2=p^2+h-2p\sqrt h\notin K.
\]

Otherwise \(g^2\) is computed exactly and must be an algebraic integer in \(\mathcal O_K\). These are the native two-cycle conditions for the fixed compact-arithmetic ground field, as in the [independent ground-field review](fresh_sqrt5_spectator_groundfield_gate_20261017.md). This step is not an assertion for arbitrary arithmetic ground fields. The code additionally checks the complete fixed-\(K\) finite-angle alphabet when \(-2<g<0\); this yields no further exclusion in the saved census.

Every one of the 149,510 unordered distinct literal-star pairs is tested. The first applicable exact exclusion is recorded:

| Reason | Pairs |
|---|---:|
| Shared occupied square pair | 24,666 |
| Positive actual-wall product | 86,932 |
| Product in \((-1,0)\) | 9,416 |
| Nonintegral native two-cycle | 19,828 |
| Native two-cycle outside fixed \(K\) | 1,008 |
| Total excluded | 141,850 |
| Unexcluded and retained as compatible | 7,660 |

An unexcluded pair is not claimed geometrically realizable. Pair certificates are compact: literal-star indices reference a reason code and, where needed, a saved exact radical-expression table. Shared square pairs are reconstructed from the saved literal stars.

## Capacities and repair

The conservative compatibility graph has one vertex for each distinct literal star and an edge for every unexcluded pair. Any simultaneous actual collection is a clique in this graph. Thus an explicitly verified proper coloring with \(c\) colors proves that at most \(c\) stars in that selected collection can occur at one vertex.

The source retains the earlier capacity-one orbit-clique construction and also colors each complete corner orbit and the union of all cohort orbits at a forest. Every saved color class is checked pairwise against the exact exclusion certificates. It uses several deterministic greedy colorings and retains the smallest certified one; no chromatic optimality is assumed. Identical coefficient vectors retain the strongest available capacity.

There are 35 capacity-one orbit-clique columns and 465 coloring columns. Their vertex costs range from one to nine. For example, the two-color certificate at forest \(4A_1+3H_3\), vertex index 690, gives

\[
 f_{26832}\le2v_{690}.
\]

The current face mass is \(28713807/5000000000\), so the allowance rises from \(19892209/10000000000\) to \(28713807/10000000000\). The union certificate at \(3A_1+A_3+D_4+H_3\), vertex 975, gives

\[
 f_{26798}+2f_{26870}+f_{26911}\le2v_{975}.
\]

All 500 exact loads are compared with their stated vertex cost. The new allowance is the maximum of the old allowance and each load divided by its cost. Only these nonnegative vertex allowances increase; no S or Q face mass changes. Therefore every inherited finite-model inequality, signed incidence column, cap, complete local packing cover, and previously certified joint constraint remains valid by monotonicity. The 36 new uniquely represented rows gain the present actual-corner constraints; they do not retrospectively gain all joint-family coefficients omitted from the preceding scan.

## Verification scope

The primary [source](../scripts/fresh_raw25_actual_corner_compatibility_20261017.py) and [saved exact certificate](fresh_raw25_actual_corner_compatibility_20261017.json) contain the complete selected provenance map, literal corner orbits, all pair decisions, proper colorings, exact coefficients, and repaired masses. Final generation 48947 and separate fresh saved replay 41747 both terminated with exit code 0 and PASS. They reconstruct these data without an optimizer. The final source and data are frozen:

- Source SHA256: `e14f9f9af4ded877905e6eb38b208cac58d74603b9a9a3d6ed2b01da425a75fc`.
- JSON SHA256: `eb78ebd28ba2adca06e8cb1af61593e2a1af61dd3202843d31cc2f4819551d90`.
- All 19 direct source hashes, seven local links, and source syntax checks pass.

An independent internal source/mathematical read passed the unique-context lifting, complete active-component embeddings, identical-wall handling, actual-corner multiplicities, arithmetic gates, proper-coloring capacities, and monotone repair. This was a scoped source/proof review, not a duplicate raw projection or orbit census, external refereeing, or formalization.
