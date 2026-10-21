# Role-specific spectator extensions preserve two exceptional payments

Two finite-component payments in the affine-tree argument survive arbitrary acute positive-semidefinite spectator extensions, even without excluding small affine cusp factors. The relevant E7 component cannot extend to finite E8 through its allowed attachment position. Its only nontrivial extension is affine E7, where the marked pair is still short. The analogous D8 fork pair extends only to affine E8 and also remains short.

These are partial payment lemmas, not a complete nonsimple-cusp quadrilateral theorem or a new dimension bound. The [checker](../scripts/fresh_cusp_spectator_role_payments_20261017.py) and [data](fresh_cusp_spectator_role_payments_20261017.json) independently calculate the inverse entries, affine kernels and small graph costs used below.

## The precise E7 role

In the D6-to-E7 configurations of the [frozen affine-tree proof](fresh_simple_cusp_affine_tree_completion_20261017.md), the new cap d attaches by a simple bond to one old fork leaf. The resulting E7 tree has arms of lengths1,2,3, with d at the end of the length-two arm. The other marked root is either the branching vertex or the leaf on the length-one arm. Their distances from d are2 and3.

Every spectator under consideration is orthogonal to every old affine-factor root. In this corner it can therefore meet this E7 only at d. With norm-two normalization, the E7 inverse entry is

` (G_E7^-1)_(d,d)=2. `

For any spectator with nonzero contact magnitude t, the Coxeter condition gives t>=1, while positive semidefiniteness of its principal extension gives

` 2−2t^2 >= 0. `

Hence t=1. The extension is the affine E7 diagram with arms1,3,3. Its kernel is strictly positive. In a positive-semidefinite full corner Gram, a kernel vector of a principal block is also annihilated by the full matrix. Since all cross products are nonpositive, that positive kernel forces every further spectator to be orthogonal to the whole affine component. Thus the component is exactly E7 or this affine E7, with arbitrary orthogonal components elsewhere. This exhausts every number of spectators and every label; it is not a bounded attachment search.

In the finite E7 case the original finite weight pays the pair. In the affine case its unchanged distance2 or3 pays it under the ideal rule below. Extending the length-two arm does not produce finite E8, whose arms would have lengths1,2,4. A generic abstract E7-to-E8 embedding is therefore irrelevant to this attachment role.

## The D8 sibling pair

In the fork-leaf P6 exception, the opposite finite component is D8. Its marked roots c,d are the two fork leaves, and every other node of this component is an old affine-factor root. Spectators may meet this D8 only at c and d. The relevant inverse block is

` [[2,3/2],[3/2,2]]. `

If both contact magnitudes are nonzero, they are at least1 and their projection cost is at least7, exceeding the available norm2. If exactly one is nonzero, its cost is2t^2, so positive semidefiniteness again forces t=1. Adding that one simple-bond spectator gives affine E8, with arms1,2,5. Its positive kernel then forces all other spectators to be orthogonal to the entire component. The complete possibilities are D8 itself or one such affine E8 extension.

The marked pair remains at distance2. In finite D8 it lies in the common fork-invariant end segment; in affine E8 it is paid by the short-distance ideal rule. No exclusion of affine E8 is used.

The construction agent independently checked the role geometry and both positive-semidefinite arguments conceptually. That read is separate from, and does not duplicate, the producer matrix calculations here.

## A factor rule with cost at most three times its order

For an irreducible affine factor of order m, use weight1 on:

* every pair if m<=7;
* pairs at distance at most3 if m>=8;
* additionally, for a classical factor, pairs in a fork-invariant end segment of order at most6.

All other pairs have weight0. The total cost is at most3m for every irreducible affine diagram. For m<=7, binom(m,2)<=3m. The classical cycles, paths, single forks and double forks of orders8 through13 each have cost exactly3m; the checker counts all24 literal graphs. For m>=14 the all-rank short-pair and end-segment formulas are those proved in the [frozen affine weight interface](fresh_simple_cusp_vinberg_weight_frontier_20261017.md), giving exactly3m. The exceptional affine E6,E7,E8 diagrams have respective costs21,20,23 and respective orders7,8,9. The affine F4 and G2 diagrams have order at most5 and are already covered by the small-order bound. The rank-one affine factor is also covered.

Consequently, if this rule and the finite weights eventually satisfy every face-payment obligation, their ideal-vertex cost would be at most3s_v V_v by the proven factor-fiber multiplicity. This sentence is conditional: these two spectator lemmas do not establish all face payments. New affine components after replacing a central root, small-affine triangle cores, and general all-finite-corner transfers in the presence of small cusps remain separate cases.

## Verification and status

Generation and a fresh saved replay pass. The checks include the literal E7 inverse diagonal and D8 inverse block, three affine extensions with positive kernels and every single-node deletion positive definite,24 small classical graphs, three exceptional affine costs, and the348 inherited graph rows. The all-spectator conclusion is the positive-kernel proof above. The all-rank cost statement uses explicit classical formulas and the finite affine classification, not extrapolation from the graph rows. Syntax, all five direct source hashes and local report links pass. The root independently reviewed the mathematical/source scope, including the two inverse-entry arguments, the positive-kernel isolation step and the all-rank cost accounting, and found no gap. That review did not duplicate the matrix calculations or graph census. This package is final and frozen; no global payment or ambient dimension conclusion is asserted.
