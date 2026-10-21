# The depth17 half rule has exact finite maximum479 at rank65

Suppose the current finite payment rule is changed only as follows: a unit
at marked distance8 in a whole finite A/B/D component becomes a half when
both marks are at distance at least17 from every endpoint, branch and
endpoint of a nonsimple bond. Every other payment and coupon stays fixed.
Then the **exact maximum over all finite Coxeter forests of rank65 is479**,
attained by D17+3D16. The finite curvature allowance is480.

This is a conditional counting theorem. It does not prove face coverage
for the modified weights or assert a new ambient dimension bound. The
[source](../scripts/fresh_depth17_rank65_conditional_budget_20261018.py) and
[certificate](fresh_depth17_rank65_conditional_budget_20261018.json) use
the completed [current finite budget](fresh_d16_d17_d14_coupon_forest_budget_20261017.md).

## Strip one large component

If every component has rank at most17, no weight changes. The frozen
budget gives maximum479, with the stated attaining forest.

Otherwise select a component of rank r>=18. Finite Coxeter classification
leaves only A_r, B_r or D_r. None participates in a current cross-component
coupon: every designated large coupon component has rank at most17 and
every small partner has rank at most6. The source checks the literal
inherited coupon function against all finite types through rank65,
including exceptional and uniform dihedral representatives.

For A_r, B_r and D_r, respectively, the numbers of changed pairs are

    max(0,r-42), max(0,r-43), max(0,r-44).

These follow by counting starting positions of an8-edge segment at the
specified distance from the path endpoints, doubled bond or fork.
The source independently constructs all144 relevant finite graphs for
18<=r<=65 and counts every eligible pair by breadth-first distances.
It checks the formulas and saves every counted pair.

The inherited doubled same-component bounds are16r-68 for A/B,
16r-64 for D18/D19, and16r-66 for D_r with r>=20. Halving one unit
subtracts one from the doubled cost. The new D bound dominates the
new A/B bounds at every relevant rank, as the source checks explicitly.

Let M(s) be the completed old unrestricted maximum at rank s, taken from
the frozen certificate. The remaining modified forest costs no more
than M(65-r), since the modification only decreases weights. Thus a
forest with a component of rank at least18 has doubled cost at most

    new_doubled_D_cost(r) + 2*M(65-r).

All48 values are computed exactly. Their maximum is953, uniquely at
r=65. Therefore this branch costs at most953/2=476.5. The other branch
has exact maximum479, and D17+3D16 still attains it. This proves the
complete rank65 result without another knapsack optimization or a
claim that the remaining forest contains only small components.

## Scope and checks

The argument retains all existing coupons, including E6/D11 and the
A1/A2--D14 role addition. It introduces none. The rank65 result would
combine with the already positive dimension65 ideal margins only after
a complete local payment theorem supplies its geometric hypothesis.

Generation52c6ab passed144 graph counts,48 large-rank bounds and28,512
zero-cross-coupon comparisons. Fresh saved replay104c66 terminated exit0/PASS.
An initial
graph-construction run detected a missing D-fork edge before any certificate
was written; the final graph has both fork leaves and asserts connectivity.
The old finite-budget certificates are bound inputs; this package does not
rerun their Bellman tables or geometric payment proofs. A separate full
source/report mathematical read passed the literal pair-count formulas,
absence of large-component cross-coupons, domination and one-component
stripping argument. That reviewer did not repeat the graph counts or the
inherited Bellman calculation. No gap was identified within this conditional
counting scope. Source, certificate and report are final/frozen.

    python3 scripts/fresh_depth17_rank65_conditional_budget_20261018.py
    python3 scripts/fresh_depth17_rank65_conditional_budget_20261018.py --verify
