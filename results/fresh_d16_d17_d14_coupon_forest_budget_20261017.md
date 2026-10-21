# Exact finite budgets after the rank16/17 and D14-role changes

With the proposed finite costs A16=116, B16=104, D16=118 and D17=125,
and the added A1/A2--D14 role coupon, the exact rank65 forest maximum
remains487. If every irreducible component has rank at most17, however,
the exact maximum is479, attained by D17+3D16. The curvature allowance
is480. The previous small-component obstruction is therefore removed
by this local rule; large finite components still obstruct the unrestricted
comparison.

The [source](../scripts/fresh_d16_d17_d14_coupon_forest_budget_20261017.py)
and [certificate](fresh_d16_d17_d14_coupon_forest_budget_20261017.json)
prove a conditional counting theorem. The geometric justification is the
separate [local transfer](fresh_finite_d16_d17_long_pair_half_transfer_20261017.md).
No ambient theorem is asserted by this budget package alone.

## The exact rule counted

Retain every payment in the frozen
[E6/D11 budget](fresh_e6_d11_quarter_coupon_forest_budget_20261017.md),
except for the four stated same-component reductions. The inherited A17
and B17 costs remain114 and105. Add a half between every mark in a whole
A1 or A2 component and the unique long-arm neighbor of the fork center
of a separate whole D14 component. That single D14 role gives total
increments1/2 and1 per component pair. It is distinct from the existing
fork-center and short-leaf half-coupons.

All finite types are checked for domination by the chosen representatives.
In particular D16 has cost118, which dominates A16=116 and B16=104;
retaining their old costs120 would invalidate this step. D17 cost125
dominates the inherited B17=105, while A17 remains separate because of
its coupons. E6 continues to represent rank6 and preserves its D11
quarter-coupon interaction. Arbitrary I2 labels retain their uniform
rank2 representative bound.

The complete1,608 large vectors consist of D10 through D16 and A17.
The interaction key has five coordinates:

    (2*d10+3*(d11+...+d16)+2*a17,
     d11+d12, d13, d11, d14).

Relative to the predecessor, the last coordinate adds d14 to the doubled
rank1 residual cost and2*d14 to the doubled rank2 cost. Rank6 still gets
the existing2*d11. Each D16 large component has doubled base236, and
the residual D17 item has doubled base250. There are626 distinct keys.

Two exact domains are computed. The unrestricted domain allows residual
ranks1 through9 and17 through75. The small-component domain allows only
1 through9 and17 as residual ranks. Together with the large-vector list,
this is precisely the domain of all forests whose components have rank
at most17. No new coupon joins two residual components.

## Upper and attaining certificates

Each of the two domains uses item-first integer knapsack generation.
Replay reads the saved tables and checks every last-item upper inequality,
an equality predecessor at every table entry, and literal component-pair
interaction formulas. It reconstructs an actual typed forest attaining
every global maximum through75 and checks every large-vector alternative.

|Combined check across both domains|Count|
|---|---:|
|Literal residual interaction formulas|125424|
|Bellman upper inequalities|1944982|
|Attaining Bellman predecessors|93900|
|Global alternatives|45008|
|Typed global attainers|152|

The unrestricted maxima never exceed the frozen predecessor's maxima.
At65 the value487 is still attained by D65, so lowering only the small
components cannot prove a general ceiling64. With component ranks at
most17 the exact value is125+3*118=479. This is a complete upper-and-attainment
calculation, not a search over selected partitions.

## Verification

Generation70315 and separate saved replay16089 both terminated exit0/PASS.
A separate full source mathematical review checked domination, the new
role's exact cost, the five-coordinate compression, both complete domains,
Bellman upper/equality verification and typed global attainment. It did
not independently replay the arithmetic. The same reviewer completed the
final report review and confirmed its conditional counting scope.

The local geometric transfer and composition with ideal margins remain
separate proof obligations for any resulting global corollary. This
package neither repeats the geometric censuses nor asserts a new chamber.

    python3 scripts/fresh_d16_d17_d14_coupon_forest_budget_20261017.py
    python3 scripts/fresh_d16_d17_d14_coupon_forest_budget_20261017.py --verify
