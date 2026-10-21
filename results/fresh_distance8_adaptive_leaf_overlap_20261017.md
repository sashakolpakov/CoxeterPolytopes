# Nine adaptive finite-leaf payments for distance-eight overlap cores

This local lemma closes nine explicitly listed high-core overlap contexts. It does **not** constitute a complete distance-eight overlap classification, and gives no new dimension bound. Generation and a separate fresh saved replay both passed. Root and `exact_hub_counts` completed separate full source/mathematical reads with no gap identified; neither review claims a duplicate nine-core arithmetic run or full overlap census.

The input is an actual quadrilateral face star, with opposite marked pairs, no marked path of length at most eight, and one of the nine literal core identifications in the [source](../scripts/fresh_distance8_adaptive_leaf_overlap_20261017.py) and [data](fresh_distance8_adaptive_leaf_overlap_20261017.json). The core has exactly its displayed edges: exclusion of further cross-core bonds is a separate obligation of the full overlap classification. Any additional supporting wall is a white root. Every edge principal diagram is positive definite, and every corner principal diagram is positive semidefinite. Normalize roots to square two.

The conclusion in every case is a **finite** corner paid by distance at most eight, or by a marked leaf at distance nine, or by a whole connected component of order at most sixteen. Thus these cases use the finite rule already budgeted in the [finite8 coupon package](fresh_finite8_a17_d11_coupon_budget_20261017.md). No new ideal pair or coupon is introduced.

## Complete one-root contact domains

All nine cores are nonsingular Lorentzian with rational, simply-laced Grams. Their four literal corner Grams are positive definite. Let `J` be the inverse core Gram and let `x` contain the nonnegative contact magnitudes of an additional white root. Its orthogonal residual is spacelike or zero, hence

`x^T J x <= 2`.

For each finite corner principal Gram `K`, the corresponding necessary bound is `x_K^T K^-1 x_K <= 2`. Every inverse `K^-1` is entrywise nonnegative. Therefore setting every nonzero contact to one gives a necessary test for each support, with no upper-label cutoff assumption. The source enumerates **every subset** of the union of eligible core nodes, discarding either a length-at-most-eight one-root shortcut or one of these monotone corner violations.

For all remaining supports, the inverse-core submatrix is entrywise positive except for one two-coordinate pattern and its renumberings. In the positive cases its unit quadratic value is a lower bound. A value above two excludes the support for all labels; a two-contact value exactly two forces both contacts to be one. In the exceptional negative-entry two-coordinate case, completing the square gives a lower bound of four times the square of one nonzero coordinate, hence at least four. This excludes it for **all real** contact magnitudes at least one; it is not a unit-label inference.

Every surviving singleton has inverse diagonal at least `4/5`. Its contact square is at most `5/2`, while a Coxeter label at least five has square at least `(3+sqrt(5))/2 > 5/2`. Thus only labels three and four remain, and their squares one and two are checked exactly. The saved domain records include all support rejection counts and the complete-square certificates. The nine graphs are constructed afresh from the frozen high-rank minimal-diagram catalogue, rather than loaded from an exploratory overlap JSON.

## Two-root exclusion and protected leaves

For additional roots with core-projection product `a < 0` and residual squares `rho,sigma >= 0`, the orthogonal complement of the Lorentz core is positive definite. If `a^2 > rho*sigma`, residual Cauchy implies their actual inner product is negative. They therefore have a native bond. Whenever the two contact paths and this bond give a marked path of length at most eight, the star is already paid. The saved exact tests use squared quantities, so contacts of magnitude `sqrt(2)` cause no approximation. Shortest paths to opposite marked sets avoid other marks internally and survive in the corresponding actual corner.

A simple marked leaf whose neighbor remains of degree two cannot belong to a classical affine component: affine A has no leaves, affine C has doubled terminal bonds, and affine B/D simple leaves meet a branch. Every chosen base component has at least ten nodes, excluding exceptional affine diagrams. The actual factor-fiber lemma says the marked component at an ideal corner is the whole varying affine factor, so unlisted ideal walls cannot repair this obstruction. Consequently such a protected leaf produces a finite corner, and its distance-nine pair receives the finite leaf payment. Core-orthogonal whites cannot alter a core leaf or its neighbor.

Four of the nine contexts have only two possible contact locations on the remote tail. Their opposite marked simple leaf and degree-two neighbor are untouched in every extension. This proves their finite distance-nine payment directly.

Three more contexts have a corner whose two marked nodes are simple leaves with already protected degree-two neighbors. The literal two-root inequalities exclude simultaneous contacts at the two leaves in an unpaid star. A single root meeting both leaves was already excluded by the one-root shortcut test. At least one marked leaf therefore stays untouched. This is an **adaptive** choice of protected leaf; neither end is asserted to be untouched in every realization.

## The two reversed five-node contexts

These are opposite orientations of the same remaining configuration. In the first numbering the marked corner is `(0,11)`, its finite base component has twelve nodes, and the remote support set is `{1,2,3}`. If any supporting root contacts that set, the residual inequalities forbid every root contacting marked node11, unless a short paid path appears. Node11 and its degree-two neighbor are then protected, giving the finite leaf-nine payment.

Suppose instead that the remote support set is untouched. The opposite endpoint3 and its neighbor2 remain a simple leaf and degree-two neighbor, so this entire marked component is finite. Its only relevant contacting root `G` has a unit contact at node11, projection square `5/4`, and residual square `3/4`. The only other possible profile that can coexist with `G` is a unit contact `B` at node13, outside the chosen base component. Its projection square is one and the projection product is `1/2`. All other profiles are excluded by positive projection and residual Cauchy. The native product `G.B` must be zero: any value at most minus one exceeds the residual Cauchy bound. Their residual Gram and inverse are

```
S     = [[3/4, -1/2], [-1/2, 1]]
S^-1  = [[2, 1], [1, 3/2]].
```

There is at most one `G`. If `G` is present there is at most one `B`: two distinct B roots would have opposite unit residuals, while G would need product `-1/2` with both. A doubled contact at the outside node has zero residual and is incompatible with G.

If B is absent, any whole connected component of core-orthogonal whites contacting G has total rank at most three. This is the established endpoint-inverse bound `m/(m+1)`: the available residual norm is `3/4`. It controls indirect attached nodes, not just G's immediate neighbors. If B is present, any core-orthogonal root contacting G costs at least the inverse entry two; equality forces a unique unit G-contact and zero B-contact. Such a root lies in the span of the two residuals. Every further root contacting B would then have positive product with it, and every root right to both is also right to it. Thus no indirect bridge reaches an arbitrary B-tail. This is the same [two-residual no-bridge argument](fresh_distance9_supporting_extensions_review_20261017.md), applied to the displayed literal matrix.

The whole finite marked component therefore has order at most `12+1+3=16`. It receives the existing whole-component payment. The reversed orientation is checked separately with its literal renumbering.

## Scope and verification

The exact source reconstructs nine cores and all36 corner Grams, enumerates the complete support domain before deriving its label bound, and verifies each adaptive pair inequality and both residual matrices. The methods are four universally protected marked leaves, three adaptive leaf choices, and two adaptive leaf-or-bounded-tail choices. The proof retains the actual-edge/corner hypotheses, the absence of short marked paths, and the separate no-extra-cross-edge premise. It does not promote the exploratory2614 count to a complete theorem, nor does it address the ideal costs of the other overlap contexts.

The final verification scope is nine reconstructed cores,36 rational corner checks, complete powerset support filtering, all surviving label3/4 profiles, and the stated adaptive/residual certificates. Root additionally inspected the inherited positive-semidefinite/corner helpers and two-residual proof. The source/data are stable; only this verification metadata was finalized after replay.
