# Unrestricted finite-volume dimension at most65

**Every finite-sided, finite-volume real hyperbolic Coxeter polytope has
dimension at most65.** The new argument excludes dimension66; the completed
[ceiling66 theorem](fresh_ideal5_weighted_dimension66_20261017.md) excludes
all higher dimensions. Arithmeticity, compactness, simple cusps, a bound
on the number of cusp factors and a bound on the number of facets are
not additional hypotheses.

The [source](../scripts/fresh_ideal5_weighted_dimension65_20261017.py) and
[certificate](fresh_ideal5_weighted_dimension65_20261017.json) assemble
the local payment theorems, complete exact finite budget and all1,659
ideal comparisons. Final verification status is recorded below.

## Simultaneous local payment rules

Start with the completed ceiling66 rule, including the finite rank17
costs A114/B105/D127 and all its existing cross-component coupons.
Apply these three refinements:

1. The [large finite D theorem](fresh_large_finite_d_leaf_half_transfer_20261017.md)
   halves the two short-fork-leaf distance9 payments in a whole D_n,
   n>=20, saving1. It adds a half-coupon between any mark of a whole
   finite component of order at most5 and a short fork leaf of a separate
   whole D10.
2. The [order14 B theorem](fresh_order14_affine_b_fork_projection_transfer_20261017.md)
   deletes two further redundant ideal pairs, lowering B86 to84. Every
   affected actual factor would require a new doubled-terminal root
   whose norm-two projection exceeds2. No finite payment changes.
3. The [order12 B/C theorem](fresh_order12_bc_unit_coupon_transfer_20261017.md)
   removes B's two remaining halves and C's two remaining nonbasic units,
   giving B64 and C63. Upgrade the existing whole A2/D11 short-fork-leaf
   halves to units. This adds2 per A2/D11 component pair; the already
   full A1/D11 payments remain full.

These refinements have compatible complete provenance. All finite
short8 units and all protected ideal units at distances at most9 remain.
The large-D theorem's affected high-overlap routes are fully partitioned:
2,531 whole components of order at most17,18 protected simple leaves,
nine adaptive simple-leaf cases,40 cases with two protected halves,
and16 cases with the new small-component/D10 coupon. Exact residual
bounds control whole supporting components, including indirect bridges.
Its universal cycle D leaf-nine alternatives have order at most19 and
are unchanged.

The order12 B/C proof reconstructs all758 eligible high patterns, all28
affected B maps and four C maps,14 literal B geometric interfaces and16
finite-payment variants. Radical or zero-residual sealing makes the
required finite components whole. The strengthened A2/D11 and existing
A1/D11 payments give finite totals1 or2 even after the ideal B half is
removed. Both affected C contexts already have whole finite B11 components;
their missing-terminal projections also give the direct contradiction16>2.
The universal immediate B alternative is excluded by the inherited strict
determinant gaps. Low, short-path, protected and all later interfaces are
checked at their actual outputs.

Neither B/C proof consumes the new D10 coupon or a large-D half. Conversely,
the large-D proof never needs a removed distance10 B/C ideal pair. The
order14 projection proof introduces no replacement payment. Thus the
three theorems apply simultaneously: every triangle receives at least2,
every quadrilateral at least1, and all weights remain nonnegative.

## Exact finite maximum495, including both coupon changes

The [D10/large-D budget](fresh_d10_coupon_large_d_forest_budget_20261017.md)
and its [full A2/D11 successor](fresh_a2_d11_full_coupon_forest_budget_20261017.md)
count the actual simultaneous rule. The successor is an abstract counting
theorem with an explicit geometric hypothesis; the completed B/C theorem
above supplies that hypothesis here.

Finite type domination retains D10 through D16 and coupon-bearing A17
as large components, with A1 through A9 and D17 onward as residual types.
D10 replaces B10: its new rank2 interaction is already larger than the
old A2/B10 coupon, which is not added again. A17 remains separate from
D17. Existing A11 coupons and the A1/D11, A1/D12 and A4/D13 upgrades
are preserved.

Write a=sum_(s=1..5) s*a_s, let d_r count D_r components, let u count
A17 components, and put t=d_11+...+d_16. The exact doubled cross cost is

    a*(2*d_10+3*t+2*u)
      +2*a_1*(d_11+d_12)+2*a_4*d_13+4*a_2*d_11.

The last term is precisely the new A2/D11 increment. All1,608 large
multisets through rank75 are covered. Their357 distinct residual tables
have860,013 Bellman upper alternatives and26,775 attaining predecessors.
All22,504 global alternatives and76 typed attaining forests are checked.
Generation uses an item-first algorithm; replay checks every last-item
upper inequality and a predecessor attaining equality. This proves
the exact maxima without importing the geometric producers' arithmetic.

Every maximum through75 remains unchanged after the A2/D11 upgrade.
At66 the maximum is495, attained by D66. Its direct graph count is493
short8 pairs plus total leaf-nine weight2. The finite curvature coefficient
is66*(66−6)/8=495. Zero finite slack suffices because the global inequality
has a strictly positive additive constant.

## Actual cusp weights and corrected complete concentration

For an affine factor of rank r, hence order r+1, the normalized cost is
bounded above by c(r):

|Rank r|c(r)|
|---|---:|
|1 through10|r/2|
|11|16/3|
|12|75/13|
|13|6|
|14|32/5|
|15|107/16|
|16|120/17|
|17 and larger|163/23|

Only the order12 and order14 maxima change: their A/B/C/D costs are
60/64/63/62 and70/84/84/82, respectively. Other envelopes are inherited.
For a cusp of link dimension d and factor ranks r_i, set

    H(t)=product_i(1+t+...+t^r_i), V=H(1),
    G=sum_j [t^j]H(t)*min(j,d−j).

The [actual factor-fiber theorem](fresh_cusp_factor_payment_framework_20261017.md)
gives V/(r_i+1) copies of each factor pair, including nonsimple cusps.
The normalized ideal margin is

    (5G−dV)/(2V) − sum_i c(r_i).

Use the inherited [all-rank convex-order theorem](fresh_distance14_unrestricted_dimension156_20261017.md):
concentrating two ranks decreases normalized curvature. For a margin
reduction, the total cost must also stay the same or increase.

Low ranks concentrate only up to10, preserving their linear cost. Two
ranks at least18 concentrate until one reaches17, preserving their
constant cost. Two rank13 parts may become12 and14: the cost increase
is11/65. Every move increases the bounded integer sum of squared ranks,
so the process terminates.

The old duplicate12 move is invalid: its cost increment is now−8/39.
Consequently the complete terminal domain retains arbitrary copies of
10,11,12,14,15,16,17, with at most one13, one part1 through9, and one
part at least18. No old terminal list is reused after its proof ceases
to apply.

At ambient66, d=65. All1,659 terminal profiles have strictly positive
margin. The exact minimum is

    13367/495924,

attained at(10,11,12,16,16). The primary computes every integer polynomial,
V and G and every rational margin. It also compares each row with the
previous [conditional scalar diagnostic](fresh_order14_b84_scalar_frontier_20261017.md);
the new B/C theorem and finite budget now supply that diagnostic's missing
geometric and finite-cost obligations.

## Strict contradiction and scope of the bound

Let D5=sum_(two-faces)(5−number of sides), and let v_f count finite
vertices. Face coverage gives

    D5 <= 495*v_f + sum_(ideal v) V_v*sum_i c(r_vi).

The [shape inequality](fresh_cusp_link_shape_curvature_20261017.md) gives

    D5 >= 495*v_f + sum_(ideal v) (5G_v−65V_v)/2 + 325.

Every ideal difference is positive, the finite difference is zero, and
the global constant325 is positive. These inequalities contradict each
other, including the all-finite and all-ideal cases. Dimension66 is
therefore impossible. The completed ceiling66 excludes every dimension
at least67, proving the ceiling65.

At ambient65 the new comparison still has two obstructions. Its finite
maximum is487 at D65, against curvature480. Among1,528 terminal ideal
profiles,16 have negative margin; the minimum is−78444/485537 at
(12,12,12,12,16). These are formal comparison obstructions, not realized
polytopes or an optimality claim. No ceiling64 is asserted.

Compact arithmetic remains at most21 and arbitrary compact Coxeter at
most29. The effective compact-facet ambient ceiling becomes65, while
the older auxiliary argument remains79. No new finite-volume example
or construction dimension follows from the proof.

## Verification status

The local theorems and both finite-budget packages have completed
generation, fresh replay and separate full internal mathematical/source
reviews. The large-D local proof also has an independent16-core rational
arithmetic/interface check. Their inherited all-label classifications and
arbitrary-support geometry retain their explicitly recorded scopes.

Primary generation and a separate fresh saved replay both pass. Full
internal source and mathematical composition reviews find no unresolved
step. An [independently implemented verifier](fresh_ideal5_weighted_dimension65_20261017_review.md)
agrees with all1,659 primary scalar rows and all1,528 lower-boundary
profiles. It also reconstructs the full1,608-vector/357-table finite
upper-and-attainment certificate and checks290 supplied B/C maps,
16 finite-payment variants and four C projection bounds. It imports
no project arithmetic; original geometric-search completeness and
all-rank structural theorems retain their separate recorded reviews.
Its saved certificate and replay status are recorded in that report.

This is a complete, internally checked computer-assisted proof. The
internal reviews were performed within the agent workflow. No external
refereeing, proof-assistant formalization or publication-priority claim
is made.

    python3 scripts/fresh_ideal5_weighted_dimension65_20261017.py
    python3 scripts/fresh_ideal5_weighted_dimension65_20261017.py --verify
