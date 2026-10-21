# Arithmetic field bounds and a stronger unconditional estimate

Two proof candidates now have independent mathematical and exact-arithmetic reviews: **compact arithmetic dimension at most20 over Q(cos(2π/7))**, and **unrestricted finite-volume Coxeter dimension at most166**. These are repository research results, not externally validated or priority-certified theorems. No new finite-volume polytope or new example dimension has been established.

| Class | Result at this checkpoint | Status |
|---|---|---|
|Compact arithmetic, ground field Q(cos(2π/7))|n≤20|New internally reviewed proof candidate|
|Compact, all Lannér subdiagrams of order≤3|n≤28|New internally reviewed fallback; includes arithmetic fields outside Q(√2,√5)|
|Compact arithmetic over Q(√2) or Q(√5)|Overall ceiling remains29|Local arithmetic restrictions sharpened; global improvement open|
|Noncompact arithmetic|n≤21, n≠20;21 is attained|Published restriction, unchanged|
|General compact Coxeter|n≤29|Published restriction, unchanged|
|General finite-sided finite-volume Coxeter|n≤166|New internally reviewed proof candidate; no arithmetic hypothesis|

The distinction between the first and third rows is essential: the overall arithmetic bound has **not** yet fallen below29. Vinberg's high-dimensional field theorem, together with the cubic-field candidate, leaves only the two quadratic fields in dimensions22–29. The [scope report](fresh_frontier_bound_scopes_20261015.md) records the published hypotheses and algorithm scope; [Nikulin's survey of the ground-field restrictions](https://arxiv.org/pdf/0708.3991) states Vinberg's field list. Allcock's published edgewalking algorithm is formulated for Lorentzian lattices over the integers; using its geometric idea over other fields or without an arithmetic lattice requires additional completeness arguments. [Allcock's paper](https://web.ma.utexas.edu/users/allcock/research/edgewalking3.pdf).

## The cubic field: a complete small-face proof

The [primary proof](fresh_compact_global_bound_20261015.md) and [independent review](fresh_compact_global_bound_20261015_review.md) give n≤20. Arithmeticity forces each minimal indefinite diagram to have at most three nodes and every such triangle to contain a7/14 bond. A quadrilateral with no corner whose marked nodes are within graph distance two is then tightly constrained.

The decisive step uses the rank of the **entire** quadrilateral star, including all white supporting walls. Rational Schur complements and the degree-three field exclude the six-node obstruction. The only seven-node exception has two corners that are entire D5 components. Assigning a half-unit to each of those two exceptional angles repairs the distance-two weighting. Every vertex has total weight at most `2n−2`, while the strict h-vector bound exceeds40 at n=21. This excludes every n≥21.

The independent checker uses a separate rational multivariate-polynomial implementation and a different complete graph enumeration. It verifies35,904 relaxed graphs, both determinant identities and their cubic remainders, the complete white-contact reduction, entire-component condition, and all-rank weight accounting. It does not import producer arithmetic. Generation and saved replay pass.

The [order-three fallback](fresh_arithmetic_highdim_20261015.md) and its [independent review](fresh_arithmetic_highdim_20261015_review.md) establish n≤28 under the broader condition that all Lannér diagrams have order at most three. This follows for arithmetic ground fields not contained in Q(√2,√5), by Vinberg's determinant-field proposition. Its distance-three weighting remains useful independently of the sharper cubic-specific argument.

## The remaining quadratic fields

The [Q(√2) triangular-face census](fresh_arithmetic_sqrt2_local_20261015.md) gives10 possible marked pairs out of21 in E7 and9 out of28 in E8. The finite angle alphabet is derived from the field, and every allowed contact profile is enumerated exactly. The report also treats contacts to other orthogonal components, which could otherwise invalidate an isolated-component exclusion.

Its [independent review](fresh_arithmetic_sqrt2_local_20261015_review.md) reconstructs all9516 profiles,32 accepted full determinants and237 replacement determinants with separate exact field arithmetic. Generation and saved replay pass; the review also verifies the ambient-component argument.

The [Q(√5) analysis](fresh_arithmetic_quadratic_bound_20261015.md) permits10 of21 E7 marked pairs and10 of28 E8 marked pairs on triangles, including arbitrary orthogonal vertex components. It separately rejects all23 frozen geometric quadrilateral witnesses:20 fail the cyclic-product field condition and3 fail integrality. Its [independent review](fresh_arithmetic_quadratic_bound_20261015_review.md) records the exact arithmetic and ambient-component checks. This rejects the saved witnesses, not every quadrilateral completion.

A further [Q(√2) quadrilateral rank lemma](fresh_arithmetic_sqrt2_quad_rank_20261015.md), with [independent review](fresh_arithmetic_sqrt2_quad_rank_20261015_review.md), excludes two disjoint Lannér cycles of orders3–5, at least one of order4 or5, each with exactly one4 bond, joined by a white bridge. Attached spherical white trees are allowed. Longer connector paths already violate the Lorentz negative-index constraint. The direct-placement census has4032 labelled cases and28 survivors of the finite-corner tests, all removed by the rank identity. Shared whites, extra contacts, other triangles and two order-three cycles remain outside scope; no global quadratic-field dimension bound follows.

The next dimension-bound problem is to control those quadrilateral completions and couple different face stars at one vertex. A triangle exclusion alone does not lower Vinberg's vertex-weight bound. No finite scout, restricted profile family or failed construction has been promoted to a global theorem.

## General finite volume:166

The [dimension-166 case binder](fresh_nikulin_coupled_charge_20261015_dimension166.md) binds the [new local proof](fresh_nikulin_coupled_charge_20261015.md), its [independent review](fresh_nikulin_coupled_charge_20261015_review.md), and the unchanged geometric dependency chain from the169 checkpoint.

Write `N=n−1` and let `M≥1` be the common cusp-height ratio. When `M≤21/20`, the local contact-cost coefficient improves from4/3 to7/6; the other branch keeps4/3. The resulting norm-two pairing lower bounds are `M(N−40)/9` and `M(N−43)/9`. At `N≥166` both contradict the strict narrow-facet threshold14. The compact-facet alternative remains143. All affine types, all ranks and every prescribed forbidden-node choice are covered.

The independent review reconstructs40 exceptional records,54 single-contact cases, the sole surviving multiple-contact support, all small-cycle cases and the symbolic classical tails. Generation, saved replay and recursive dependency binding pass. The proof uses no arithmeticity or simple-cusp assumption. Prokhorov's published995 and Nikulin's announced bound below300 are kept distinct from this repository candidate in the binder.

## Preserved construction work and archive

The [preceding checkpoint](contact_cost_bound_and_cusp_links_2026_10_15.md) remains frozen. The compact H9 catalogue still excludes164 of240 recipes; its T10 continuation stays at9022 visits. The33-wall H9 chamber and47-wall H20 chamber retain their exact infinite-volume escapes. Their finite vertices, finite edges and bounded cusp links are unchanged; none is a completed new example.

The [archive audit](arithmetic_field_bounds_and_coupled_charge_2026_10_15_audit.json), produced by the [archive checker](../scripts/audit_arithmetic_field_coupled_charge_checkpoint_20261015.py), preserves the previous source hashes, checks the new recursive manifests and local links, and binds updated README and research-plan snapshots. Archive integrity is separate from mathematical proof review.

```sh
python scripts/fresh_compact_global_bound_20261015_review.py --verify
python scripts/fresh_nikulin_coupled_charge_20261015_dimension166.py --verify
python scripts/audit_arithmetic_field_coupled_charge_checkpoint_20261015.py --verify
```

Fresh prefixes frozen by this checkpoint are:

* `fresh_nikulin_coupled_charge_20261015`
* `fresh_arithmetic_highdim_20261015`
* `fresh_compact_global_bound_20261015`
* `fresh_arithmetic_sqrt2_local_20261015`
* `fresh_arithmetic_quadratic_bound_20261015`
* `fresh_arithmetic_sqrt2_quad_rank_20261015`

Use new prefixes for subsequent work. Historical pending-review sentences are retained where a frozen producer artifact predates its independent review.
