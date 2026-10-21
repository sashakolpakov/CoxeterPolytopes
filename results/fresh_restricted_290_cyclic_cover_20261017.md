# The restricted witness survives every local cyclic-product test

The stronger hypothetical witness from [the complete unowned-pair frontier](fresh_restricted_290_unowned_pair_frontier_20261017.md) also survives all Vinberg cyclic-product membership and integrality necessities within its saved common finite-vertex configurations. Its exact ratio remains

\[
\frac{18670750401222348}{1549726013368729}
=12.047775051950415\ldots>12.
\]

Every one of its1,415 positive-mass cover terms is checked, including singleton-star terms. No new incompatibility or actual capacity column is asserted. The actual full and restricted witnesses remain unchanged. This is a certificate that these local arithmetic necessities cannot eliminate the surviving relaxation witness, not a proof that its configurations assemble into a physical polytope.

The [source](../scripts/fresh_restricted_290_cyclic_cover_20261017.py), [saved exact certificate](fresh_restricted_290_cyclic_cover_20261017.json), and [separate standard-library saved-data checker](../scripts/check_restricted_290_cyclic_saved_20261017.py) are provided. No optimizer is used in generation or verification.

## Current covers and complete physical frames

The predecessor reconstructed all571 forests of the current290-class domain, retained every historical and explicit owner outcome, and hypothetically deleted all570,398 genuinely unowned pairs. The two recent owned unknowns at forest1763 remain compatible. This successor rebuilds the same ownership and graph states; it does not retest an owned unknown and silently change its outcome.

For each forest the program takes either the predecessor's attaining exact cone cover or its explicit feasible cover. Zero-demand stars are removed from each term. This cannot destroy compatibility, increase the cost, or reduce a positive-demand marginal. The saved terms are merged only when their literal star sets agree, and their positive marginals and exact total cost are verified. The multiplier is `max(old_mu, cover_mass)`, exactly as in the predecessor.

Each term is reconstructed in its complete13-root positive-definite finite-vertex frame. Every distinct external wall is identified by its full native contact covector. Its projection square and strictly positive timelike height square are recomputed. With the outward future sign fixed, a contact covector determines a unique external normal. Identical covectors therefore represent one wall and are deduplicated before any cycle test.

For each star, the reconstructed pair's actual omitted product is checked against its frozen native value. For distinct external walls with covectors p,q and height squares h_p,h_q, the actual product is

\[
g_{pq}=p^T G^{-1}q-\sqrt{h_ph_q}.
\]

Every such product is recomputed in an immutable exact field selected from that pair's native covectors. No coefficient array from a previous field is reinterpreted. The existing acute-sign, forbidden-angle-gap, fixed-K squared-product, algebraic-integrality and finite-angle-alphabet tests all pass in these chosen configurations. Here K=Q(sqrt(5)), and all roots have norm2, so the symmetric Gram entries are the normalized Vinberg entries.

The graph used for the cycle certificate includes **all** nonzero inner products among the13 basis roots and the external walls: basis-basis edges, every native wall-basis contact, and every external-external edge. Passive basis components are retained even though their external contacts vanish. Thus the check includes cycles through the finite frame, not merely triangles formed by three external walls.

## Exact all-cycle certificate

Pair-square tests alone do not imply all cyclic-product conditions in general. This package proves the conditions separately for each saved configuration.

Write s_ij=g_ij² for each nonzero edge. The exact physical checks give s_ij in O_K. In each connected nonzero-edge component, choose a spanning tree, root it, and define t_i as the product of edge squares along the tree path from the root to i; the root has t=1. The certificate stores the tree, all t_i, and for every edge a positive element c_ij in K satisfying

\[
c_{ij}^2=t_i t_j s_{ij}.
\]

For a tree edge this relation is automatic. For a chord, repeated initial tree-path edges occur twice, so the relation says precisely that the squareclass product on its fundamental cycle is trivial. Fundamental cycles generate the cycle space. Equivalently, modulo K-squares each edge has class[t_i][t_j], whose product around any cycle is trivial because every vertex occurs twice. Hence every cyclic Gram product lies in K. Signs do not affect membership, and zero-edge cycles have product0.

Each g_ij is integral over O_K because it satisfies the monic polynomial X²-s_ij. A product of such entries is integral over O_K. Since the cyclic product lies in K and O_K is integrally closed, it belongs to O_K. This proves all cycle lengths at once; it is not a finite check restricted to triangles.

The program also checks the displayed c_ij are algebraic integers. That is a redundant exact consistency check: their squares are products of integral edge squares, and c_ij already lies in K.

## Complete result and preserved witness

|Checked domain|Count|
|---|---:|
|Current forests|571|
|Positive cover terms|1,415|
|Singleton-star terms|637|
|Terms with multiple stars|778|
|Distinct external walls, summed over forests|3,319|
|External-wall occurrences across cover terms|3,887|
|Native star omitted-product checks|2,270|
|Distinct external-wall pair checks, summed over forests|3,731|
|Nonzero-edge squareclass certificates|18,157|
|Sum of cycle-space dimensions|1,866|

Every check passes. The covers still dominate all required actual-corner marginals at the fixed restricted face masses. Their multiplier sum remains 1549726013368729/1544400000000000. All earlier model columns and capacities remain satisfied because the face masses are fixed and these hypothetical multipliers dominate the actual predecessor multipliers.

Consequently adding local incompatibilities solely because a common-frame configuration fails a cyclic-product necessity cannot remove this above12 witness: the displayed covers use configurations passing every such necessity. This conclusion concerns the current local corner-capacity relaxation. It does not claim universal sufficiency of pair tests, physical realization of a cover term, consistency of wall identities between different finite vertices, global cycles involving walls from different frames, complete implicit packing or face coverage, target feasibility, a finite-model optimum, or a dimension theorem.

## Verification and review

Generation session84150 completed with exit0/PASS in86.33 seconds. Separate fresh native saved replay session72623 completed with exit0/PASS in77.01 seconds. It rebuilds the complete current frames and every stated physical arithmetic check and compares the resulting certificate exactly; it runs no optimizer.

The independent standard-library saved-data checker imports no producer or project arithmetic. Its final-location run passes all18 direct bindings,1,415 cover/marginal records,3,731 saved pair squares,2,270 native omitted-square identities,18,157 exact K-squareclass equations, all tree-parent and cycle-rank identities, and the multiplier/ratio reconstruction. It independently rebuilds each term's edge list from the saved native contacts, cross-wall squares and common saved basis edge list. This is certificate/accounting verification; it does not independently reconstruct the finite basis or the physical projections and field-sign tests.

The root and construction agents read the complete producer source and final report, including chosen covers, native identity and positive-height reconstruction, immutable pair fields, all-edge squareclass certificates and multiplier preservation. The root also read the separate saved-checker source. These source/mathematical/report reviews do not duplicate the571-frame arithmetic replay or projection reconstruction. The construction review does not claim execution of the saved checker.

All predecessor files remain unchanged. A subsequent structural task can inspect cross-vertex consequences of a concrete surviving configuration: forest1325=3A1+2H3+H4, stars12 and46, a raw25/raw23 pair of mass62840451/5000000000 using three distinct external walls. They are face26918/provenance541 and face26448/provenance490, sharing one wall. The common passive H3+H4 has rank7; the remaining nine active normals lie in Lorentz dimension7. A full finite vertex must retain all passive roots, so only84 six-root active subsets need examination. These may force vertices beyond the already-counted two prism stars. That census is not performed or claimed here.
