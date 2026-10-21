# Independent review of the raw23 finite-vertex clique family

FINAL/FROZEN. Source/proof review PASS. The separate saved-certificate checker passed generation and a fresh saved replay. No mathematical or counting issue was found in the stated family.

This reviews the [primary family](fresh_dominant_cap_compatibility_family_20261017.md), which produces 294 necessary incidence inequalities from 57 positive raw23 prism classes. Its repaired fixed-witness ratio is approximately 8.804856730635848. This is not an optimized feasible rule, a complete face catalogue, or a new dimension theorem.

## Whole-face coefficients and repeated presentations

At a fixed simple finite vertex, its complete labelled incident frame and an unordered pair of external contact covectors determine the full wall star. Every genuinely external actual wall has strictly negative product with the future vertex, so its timelike sign is fixed. Identical covectors therefore mean one actual wall. Identical unordered pairs mean one local full star, regardless of the number of isomorphism presentations.

The source deduplicates these full-star keys before exclusion tests. Complete orbit sets are merged when equal and asserted disjoint otherwise; a merged orbit is asserted to belong to only one selected Q class. It also checks that duplicate keys have the same occupied quadrilateral pairs. Thus no local incidence is counted twice merely because it has two presentations.

The whole-face coefficient instead counts actual finite vertices of an abstract prism belonging to the selected orbit types. The ten distinct omitted-root pairs describe its ten actual vertices. Equivalent corner types at different prism vertices still contribute their multiplicity globally; they are not simultaneous occurrences at one vertex. Summing the local capacity-one inequality over all vertices therefore gives precisely the saved integer face coefficients, including coefficient two when two distinct prism vertices qualify.

## Orbit completeness, products and clique capacities

The componentwise mapping algorithm uses every whole matching target component, including each choice among identical components and every internal weighted graph isomorphism. These are exactly the finite-forest automorphisms relevant to the five retained raw roots. Unmapped spectator permutations change neither complete external covectors nor occupied port pairs. The orthogonality and complete-component assertions support that quotient.

Passive components with normalized entries outside Q(sqrt5) are harmless here: every exported covector is explicitly supported on the inverted 3/5-labelled components. The positive definite incident Gram is block diagonal, so the other inverse blocks make zero contribution to every projection.

For complete covectors p,q, the product is

\[
p^T V^{-1}q-\sqrt{(p^T V^{-1}p-2)(q^T V^{-1}q-2)}.
\]

Both squared heights are strictly positive. The source's radical comparison first determines the sign of the expression before squaring, so it introduces no false root or future-side ambiguity. Identical single-wall covectors are skipped before applying the positive-product or (-1,0) gap tests. The internal product of each literal star is also checked against its prescribed omitted-root product in the primary.

Distinct full stars occupying the same quadrilateral port pair cannot coexist: that port pair determines one two-face, whose full supporting-and-side wall set is unique. The other exclusions use either a strictly positive actual product or the universal norm-two Coxeter gap (-1,0). A complete orbit internally excluded pairwise has capacity one. A clique of such orbits has capacity one when every cross-orientation pair is excluded. Unknown pairs remain allowed; omitting their orbits or edges is conservative.

## Independent saved-certificate checks

The [review source](../scripts/fresh_dominant_cap_compatibility_family_20261017_review.py) imports no producer arithmetic or graph helpers. Its [certificate](fresh_dominant_cap_compatibility_family_20261017_review.json) verifies:

- all 570 actual corner records, their unique ownership and their whole-face multiplicities;
- 6,010 distinct literal full-star keys and their disjoint orbit partition;
- 560 internal capacity-one orbits and all 378 universal cross-orbit edges from the saved pair witnesses;
- all 294 clique columns, maximality in the saved graph, integer coefficients and unit vertex costs;
- all 22,397 shared-quadrilateral exclusions, with literal occupied-pair intersections;
- all 112,038 geometric exclusions through separate Fraction arithmetic on their 82 distinct radical expressions;
- all 120 maximum allowance increases, the exact normalizer and the exact repaired ratio.

The 134,435 proved exclusions and ten unresolved tested pairs agree with the primary. Separate rational arithmetic evaluates signs in Q(sqrt5), then the nested physical radicals with their positive branches. No floating-point sign or angle decision is used in this review.

This does **not** independently reconstruct the 6,010 projection covectors, rerun the component-embedding census, or establish completeness of the underlying raw23 geometric catalogue. Those steps received the source/proof review above. The checker validates the complete exported counting/clique/sign certificate; it does not turn the selected family into full geometric coverage. No optimizer or large historical replay was run.

All face masses remain fixed and vertex allowances only increase, so the primary's inherited constraints are preserved by monotonicity. Syntax, direct manifest bindings and local links pass. These are internal checks, not external refereeing or formal proof-assistant verification.
