# Direct rank-five Gale candidates for compact H9

The bounded search produced **1,063 saved realizable fifteen-facet nerves outside every standard single-face-truncation family**. All 1,063 have exact inherited diagram obstructions: 902 literal local-pattern applications and 161 two-Lannér-block applications. No candidate reached a full ambient label search or metric fit. This closes the saved proposal set, not the fifteen-facet class.

The primary data are [the terminal search ledger](fresh_h9_rank5_gale_20261014.json), [the separate arithmetic audit](fresh_h9_rank5_gale_20261014_independent.json), and [the final certificate with an explicit rational realization](fresh_h9_rank5_gale_20261014_certificate.json). The resumable [state file](fresh_h9_rank5_gale_20261014_state.json) is terminal; no job remains running.

The deterministic seed was `2026101406`. The completed budget was 3,000 balanced integral initialization proposals followed by four walks of 10,000 mutations each. Floating-point annealing selected the next state; every retained geometry and exclusion uses exact arithmetic. Temporary states with a large minimal nonface were permitted as stepping stones and were never treated as compact Coxeter candidates.

| Proposal outcome | Count |
|---|---:|
| Minimal nonface of order greater than five | 30,233 |
| A label does not occur in an actual facet | 410 |
| A vertex link has a simplex-boundary join factor | 885 |
| Passes these filters, including repeated nerves | 11,472 |
| Total proposals | 43,000 |

The saved representatives comprise one initialization case and 230, 344, 259, and 229 new cases from the four walks. Each has fifteen dual facets and between 228 and 366 actual vertices. Their literal minimal-nonface sets are distinct, and their keys under the existing colored-clutter `pynauty` canonicalization are distinct. The separate audit verifies every recorded relabeling and key; it does not independently implement graph-isomorphism canonicalization. The exclusion result does not require the representatives to be pairwise nonisomorphic.

For realizability, each saved matrix has fifteen rows and rank five. All 3,003 five-row minors are nonzero. Alternating cofactors enumerate all 5,005 six-row circuits; the positive circuits cover every label. Summing their positive dependences gives an explicit strictly positive integral vector `w` with `Gᵀw=0`. Thus `diag(w)G` has zero row sum. Its transpose kernel has dimension ten and contains the constant vector, giving fifteen rational points in affine dimension nine. Gale duality identifies precisely the complements of positive six-circuits as facets. They have nine vertices and are eight-dimensional simplices. Every label occurs on a facet, and every ridge has exactly two incident facets. The complete minimal-nonface clutter is checked, with every member of order two through five.

The final certificate also constructs the actual fifteen rational points for case 0, together with all 310 supporting affine functionals. Each functional is zero on exactly its nine-vertex facet and positive on the other six points; each facet has affine rank eight. This is a Euclidean polytopal realization, not a hyperbolic Coxeter realization.

The exclusion of standard truncation families uses a sufficient structural test rather than a comparison only against the old 45,679 nerves. A face truncation of a simple polytope is dual to stellar subdivision at a simplex `F`. Its new vertex has link

`boundary(F) * link(F)`.

In the minimal-nonface clutter of this link, `F` is a single member disjoint from every other member. None of the fifteen vertex links of any saved nerve has such a member. Consequently no saved nerve is isomorphic to the result of one standard face truncation of any fourteen-facet parent. This does not exclude other face surgeries or imply that the nerve admits a Coxeter metric.

The exact diagram gates are independent of the ambient facet count. Each local pattern is used only after an explicit relabeling identifies its entire induced minimal-nonface clutter. The [inherited local proof package](fresh_h9_newfacet_frontier_20261014_local_independent.json), [seven-node relation](fresh_h9_newfacet_frontier_20261014_relations.json), and [additional fifteen-facet patterns](fresh_h9_fifteen_facet_20261014_patterns.json) exhaust the relevant necessary Lannér and spherical constraints. The [two-block proof](fresh_h9_bugaenko_remaining_20261014_block_gate.json) applies only to disjoint blocks of orders `(3,4)`, `(3,5)`, `(4,4)`, `(4,5)`, or `(5,5)` whose union has exactly those two minimal nonfaces. It supplies a second negative direction, contradicting Lorentzian inertia. No `(3,3)` claim is made. These inherited proofs include every finite integer label: labels two through six are literal, and category seven covers all finite labels at least seven, with the corresponding interval proof where a metric value is required. Missing triangles are allowed; no fourteen-facet-specific restriction is imposed.

The fresh independent audit uses integer Bareiss determinants rather than the proposal engine's bounded `int64` Leibniz sums. It reconstructs all positive circuits, all downward face sets, every minimal nonface, and every vertex link. Its link test uses minimalized ideals obtained from the whole nerve, independently of the producer's explicit facet-link construction. It checked **3,192,189 determinants, 5,320,315 circuits, 342,574 listed polytope vertices, and 15,945 vertex links**, and rechecked every literal obstruction application. It inherits the earlier obstruction-library proofs; it does not repeat their label searches or determinant/Bernstein arithmetic. There were no new full-label searches to replay, no node-limited saved cases, and no unresolved saved frontier.

The final ledger and explicit realization were generated and freshly replayed with:

```sh
python scripts/fresh_h9_rank5_gale_20261014_certificate.py --verify
```

The completed separate all-case replay is [the independent checker](../scripts/fresh_h9_rank5_gale_20261014_independent.py). It ran with two workers while the bounded search was producing immutable cases. The short final command binds that finished audit, every primary case, all per-case audit records, and the terminal state; it is not a second complete determinant replay. Sources are [the proposal driver](../scripts/fresh_h9_rank5_gale_20261014.py), [exact geometry and gate helpers](../scripts/fresh_h9_rank5_gale_20261014_core.py), and [the standalone final certificate](../scripts/fresh_h9_rank5_gale_20261014_certificate.py).

This provides a reproducible family genuinely outside the previous truncation construction, but no compact H9 example. A useful next search would penalize the actual induced obstruction supports during generation, rather than continue the same walk using only missing-face size and truncation filters. Any changed proposal must still pass the full realization and diagram checks; no claim about unsampled Gale configurations follows from this batch.
