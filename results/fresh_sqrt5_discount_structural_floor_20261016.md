# Exact limits of the current payment family

The [small exact certificate](fresh_sqrt5_discount_structural_floor_20261016.json) and its [checker](../scripts/fresh_sqrt5_discount_structural_floor_20261016.py) show why further context splits of the present rule cannot reach the low teens. These are obstructions to a specified payment model, not existence results or universal limitations of the weighted two-face method.

For a simple convex d-polytope, write its nonnegative symmetric h-vector as h_i. The identities f0=Σh_i and f2=Σ binom(d−i,2)h_i give

`f2/f0 > min_i (binom(i,2)+binom(d−i,2))/2 = floor((d−1)^2/4)/2`.

The inequality is strict because h0=hd=1 and the endpoint coefficient exceeds the minimum. Requiring two units on every triangle, one on every quadrilateral, and nonnegative payments on longer faces gives total payment at least `5f2−binom(d,2)f0`. Thus the vertex-budget threshold is

`T(d) = 5 floor((d−1)^2/4)/2 − binom(d,2)`.

The checker independently constructs A1, A2, D4 and D5 Gram matrices, proves them positive definite by exact elimination, and counts every unordered root pair at distances one and two.

| Dimension | Formal finite vertex forest | Unchanged baseline | Required vertex budget |
|---|---|---:|---:|
| 17 | 3D5 + A2 | 25 | 24 |
| 14 | 2D5 + D4 | 22 | 14 |
| 13 | D5 + 2D4 | 20 | 12 |

None of these forests has a component of rank six. Every existing long-component discount is therefore zero, and every available compensation is nonnegative. No adjustment of those compensation contexts or coefficients can bring these vertex costs within the required budgets. This rules out using that family to exclude dimensions 17, 14 or 13; it does not rule out the dimensions geometrically.

Allowing discounts on small D4/D5 adjacent pairs is a substantive change. Even if every adjacent pair in the displayed forests received one uniform discount, dimensions 14 and 13 would respectively require discounts at least 8/11 and 4/5. A cap of one half, useful for the current quadrilateral pruning, cannot support that goal.

There is a further limit to an adjacent-only strategy. At dimension 11, the forest 2D5+A1 retains eight distance-two payments even if every adjacent payment is reduced to zero; T(11)=15/2. Excluding dimension 11 through this per-vertex budget therefore requires changing distance-two payments or imposing additional geometric/incidence restrictions. This does not prevent an upper bound of 11–13 by other means.

The existing [996-state triangle catalogue](fresh_sqrt5_triangle_extension_complete_20261016.md) is more reusable than its original discount range suggests. Its hereditary pruning depends only on the baseline and the number of literal black label-three edges, not on the rank-six eligibility threshold. A compact triangle has at most two such edges. Consequently its three stopping rules remain safe for any adjacent-only deduction of at most one per corner: baseline at least four; baseline three with at most one such edge; or baseline two with none. The checker compares these conditions with the frozen source's actual pruning function. A fresh literal payment evaluator is still required. This short algebraic extension does not alter the historical catalogue's claimed parameter range.

Quadrilateral coverage changes substantially for these larger discounts. The current baseline-one closure suffices only when every discounted corner loses at most one half and all added payments are nonnegative. Larger small-component discounts require additional baseline-two, baseline-three and potentially baseline-four domains. No coverage or successful weight assignment for those enlarged rules is asserted here.

Generation and a fresh saved replay pass. The certificate gives a concrete reason to pursue payments on small D components, global short-face compatibility, or higher-dimensional face information, rather than more context splits inside the unchanged family.
