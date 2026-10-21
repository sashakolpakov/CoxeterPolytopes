# Exact finite-forest budget with D10 coupons and large-D half payments

For the specified finite payment rule, the exact maximum at total rank66 is **495**, attained by the whole finite diagram D66. This equals the finite curvature coefficient at that rank. The computation covers every rank from0 through75 and retains all existing cross-component coupons.

This is a finite counting theorem conditional on the separately proved geometric payment rules. It does not recompute ideal-vertex costs or assert a new ambient dimension bound.

[Independent source](../scripts/fresh_d10_coupon_large_d_forest_budget_20261017.py) · [exact certificate](fresh_d10_coupon_large_d_forest_budget_20261017.json)

## Payment rule and type domination

The input [large-D transfer](fresh_large_finite_d_leaf_half_transfer_20261017.md) changes the whole finite D_n cost to8n−33 for n≥20 and adds a half-payment between every mark of a whole finite component of rank at most5 and each short fork leaf of a separate whole D10. D18 and D19 keep8n−32. The [rank17 theorem](fresh_finite_d17_long_pair_support_20261017.md) supplies A17=114 and D17=127; the weaker upper bound B17≤127 is enough here. Other same-component costs and existing coupons are unchanged.

The finite Coxeter classification and the unchanged base costs are inherited inputs. The checker separately compares every same-rank representative and every coupon partner among228 type representatives, including I2(*) for all labels at least4. Its51,984 comparisons establish the following domination:

- Ranks1 through9 are represented by A of that rank.
- Ranks10 through16 are represented by D of that rank.
- Rank17 retains A17, which has cross coupons, and D17, which has none; B17 is dominated by D17.
- Ranks18 through75 are represented by D of that rank.

Each replacement preserves rank and can only increase both the base cost and every coupon interaction. Applying the replacements one component at a time therefore gives an upper bound for every finite forest.

In particular, B10 has base45 and its old A2 coupon has total cost1. D10 also has base45, but the new A2/D10 coupons have total cost2. Thus D10 dominates B10. The representative formula includes the new D10 contribution and **does not add the old B10 coupon a second time**. The old A1/A11 unit and A2/A11 half coupons are explicitly included before type domination; D11 dominates them as well.

Let d_r count whole D_r components, u count A17 components, t=d_11+⋯+d_16, and let a_s count residual A_s components for1≤s≤5. Write a=Σs a_s. In doubled units the representative cross cost is exactly

\[
a(2d_{10}+3t+2u)+2a_1(d_{11}+d_{12})+2a_4d_{13}.
\]

The same-component costs are doubled too, so the entire calculation uses integers.

## Complete domain and exact certificates

The implementation imports no project or producer arithmetic. It enumerates large components by their number and a nondecreasing rank tuple, rather than using the earlier producers' consecutive-rank count recursion. There are exactly1,608 multisets of D10,D11,…,D16,A17 with total rank at most75. This is exhaustive because every such component has rank at least10, so there are at most7.

Their residual interaction depends only on

\[
(L,b,c)=(2d_{10}+3t+2u,\ d_{11}+d_{12},\ d_{13}),
\]

giving152 distinct tables. For each original multiset and every residual component, the checker compares this compressed formula against the explicit type-pair coupon function. No dummy component is used to represent the total rank.

The remaining components are A1,…,A9 and D17,…,D75. They have no mutual cross coupons. Once a large multiset is fixed, their contributions are additive. An item-first unbounded-knapsack calculation produces each residual table. A separate verifier checks, for every fitting rank r and total n,

\[
V(n)\ge V(n-r)+c(r),
\]

and checks a saved predecessor attaining equality. Induction from V(0)=0 proves both the upper bound for every residual forest and attainment. Fresh replay uses these saved certificates; it does not rerun the maximizing algorithm.

At each total rank the global certificate dominates every fitting large-multiset/residual-table combination. An explicit typed forest is independently recounted using the original component-pair coupon function and attains the reported maximum. Every attaining component is an actual representative with its stated exact cost, so using B17≤127 as an upper bound does not turn these maxima into merely unattained bounds.

| Certificate part | Count |
|---|---:|
| Large component multisets |1,608|
| Compressed interaction tables |152|
| Direct interaction-formula comparisons |109,344|
| Bellman upper-bound alternatives |366,168|
| Attaining Bellman predecessors |11,400|
| Global upper-bound alternatives |22,504|
| Typed global attainers, including rank0 |76|

## Maxima and direct large-D count

The finite curvature coefficient is β_n=n(n−6)/8 for even n and (n−1)(n−5)/8 for odd n. The following values are exact.

| Rank | Maximum | Typed attainer | β_n minus maximum |
|---:|---:|---|---:|
|64|480|4D16|−16|
|65|487|D65|−7|
|66|495|D66|0|
|67|503|D67|17/2|
|68|511|D68|16|
|69|520|A5+4D16|24|
|70|527|D70|33|
|71|535|D71|85/2|
|72|543|D72|51|
|73|551|D73|61|
|74|560|2A5+4D16|69|
|75|567|D75|161/2|

The D_n formula also has a direct graph proof. Removing one short fork leaf leaves a path of order n−1. For n≥18 that path has8n−44 pairs at distance at most8; the restored leaf adds9, giving8n−35. Exactly three leaf-marked pairs have distance9. The long-end pair keeps weight1. For n≥20 the two short-fork-leaf pairs have weight1/2 each, giving8n−33 in total; for n=18,19 all three still have weight1, giving8n−32. The checker builds each D graph from18 through75 and independently counts these pair classes. In particular D66 has493 short pairs and total leaf-nine weight2, hence495. The argument, rather than the finite graph check alone, proves the stated all-rank formula.

## Validation and scope

Generation and a separate fresh saved-certificate replay both terminated directly with exit0/PASS. The logs are `/private/tmp/d10_large_d_budget_gen.log` and `/private/tmp/d10_large_d_budget_verify.log`. The source binds the final local large-D theorem, the rank17 costs, the earlier finite-budget certificate and relevant coupon inputs. Verification does not invoke an optimizer or import their arithmetic functions.

The parent reviewer read the entire source and report. The inherited finite-classification scope, all228 typed partner comparisons, separate A17 role, preservation of A11 coupons, D10 domination without B10 double counting, complete1,608-vector domain and152-key compression, Bellman upper bounds and attainment, global typed witnesses and the direct D66 count493+2=495 all passed. This was a full source/mathematical read, without duplicating the certificate arithmetic.

Two other reviewers separately read the full source and report. Both passed the coupon accounting, distinct rank17 types, complete multiset/residual decomposition, interaction compression, all-last-item inequalities with attaining predecessors, typed global exactness and all-rank D count. They did not repeat the arithmetic, Bellman verification or geometric census.

Nineteen direct bindings, Python syntax and all four local report links pass. Source, data and report are final/frozen. The proof of geometric face coverage, the inherited finite classification and earlier coupon-transfer theorems remain separately bound inputs. No ideal scalar comparison or new dimension conclusion is asserted by this package.
