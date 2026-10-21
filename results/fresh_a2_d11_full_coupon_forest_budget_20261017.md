# Conditional finite budget for full A2/D11 short-leaf coupons

Upgrading the existing A2/D11 short-fork-leaf half-payments to units leaves **every finite-forest maximum through rank75 unchanged** from the [frozen D10/large-D budget](fresh_d10_coupon_large_d_forest_budget_20261017.md). In particular, the maximum at rank66 remains495, attained by D66, with finite curvature slack0.

This is an exact conditional counting theorem. It does not prove the proposed coupon transfer or any ideal-pair reduction that uses it, and does not assert an ambient dimension bound.

[Independent source](../scripts/fresh_a2_d11_full_coupon_forest_budget_20261017.py) · [saved certificate](fresh_a2_d11_full_coupon_forest_budget_20261017.json)

## The precise change

A whole A2 has two possible marked roots, and a whole D11 has two short fork leaves. The four existing half-payments therefore have total cost2. Making each a unit raises that total to4, an increment of2 for each A2/D11 component pair. In doubled units this is exactly4a_2d_11.

All existing base costs and cross coupons remain. In particular the new D10 coupons still replace the dominated B10 contribution without adding it twice; the A1/D11 full upgrade is already present and receives no further increase. With the predecessor's notation the doubled representative cross cost becomes

\[
a(2d_{10}+3t+2u)+2a_1(d_{11}+d_{12})+2a_4d_{13}+4a_2d_{11}.
\]

The rank-preserving type-domination comparisons are repeated with this changed coupon function. Replacing a small rank2 component by A2 can add the new coupon and cannot lower the cost. Replacing an A11 or other rank11 component by D11 also remains valid. Thus the same representative types and complete1,608 large-component multisets apply.

## Complete revised certificate

The source uses the predecessor's independently reviewed standard-library method and imports no project arithmetic. It generates the full large domain again, then compares each multiset, used rank, base cost and original three interaction coordinates with the frozen predecessor. The additional fourth coordinate d_11 distinguishes the new coupon; there are357 distinct interaction keys.

For each key, the residual costs of A1,…,A9 and D17,…,D75 are recounted from the original component-pair coupon function. Residual components still have no mutual coupons. Item-first unbounded-knapsack generation is verified by every last-item upper-bound inequality and a saved attaining predecessor. Every global candidate is bounded, and every reported maximum has a typed forest whose full component-pair cost attains it. Fresh replay checks the saved certificates without rerunning maximization.

| Certificate part | Count |
|---|---:|
| Large component multisets |1,608|
| Revised interaction tables |357|
| Direct interaction-formula comparisons |109,344|
| Bellman upper-bound alternatives |860,013|
| Attaining Bellman predecessors |26,775|
| Global upper-bound alternatives |22,504|
| Typed global attainers |76|

The checker additionally asserts equality of every global maximum with the predecessor, not merely that adding coupons cannot reduce it. Representative examples are:

| Rank | Maximum | Typed attainer | Finite curvature slack |
|---:|---:|---|---:|
|65|487|D65|−7|
|66|495|D66|0|
|67|503|D67|17/2|
|68|511|D68|16|
|69|520|A5+4D16|24|
|74|560|2A5+4D16|69|
|75|567|D75|161/2|

The original same-component graph formulas and typed-cost exactness have not changed. The predecessor's all-rank D count, and the direct D18–75 graph checks repeated here, give D66 cost493+2=495.

## Validation and scope

Final generation and a separate fresh saved-certificate replay both terminated directly with exit0/PASS. Logs are `/private/tmp/a2_d11_budget_gen.log` and `/private/tmp/a2_d11_budget_verify.log`. The geometric inputs to the predecessor are frozen and hash-bound. The proposed additional full A2/D11 transfer is deliberately an explicit hypothesis of this package; it is not silently treated as proved by the cost calculation.

The order16 reviewer read the entire successor source and report, including the final equality checks for all76 maxima. The precise four half-to-unit leaf increments, unchanged fork-center coupons, fourth interaction coordinate, rank2 residual cost, repeated type domination, Bellman/global upper bounds and attainment, predecessor comparisons and conditional scope all passed. This was a source/mathematical review, without numerical replay.

The parent reviewer read the complete successor delta and report against the already reviewed predecessor source. The doubled4a_2d_11 increment, fourth key, changed rank2 residual cost, repeated partner domination, complete global typed certificate and equality of all76 maxima passed. This review did not repeat the arithmetic. Twenty-two direct bindings, Python syntax and all three local report links pass. Source, data and report are final/frozen.

No ideal budget, new local face-payment theorem or global dimension conclusion is established here. A later composition may combine this exact budget with a separately completed coupon-transfer theorem and ideal-vertex comparison.
