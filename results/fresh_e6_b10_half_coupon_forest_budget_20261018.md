# Enlarged E6/D11 and A1/B10 coupons retain the finite maximum464

At rank64 the exact maximum finite cost remains **464**, attained by the
forest of four D16 components, after both of these additions:

- Upgrade the typed whole E6/D11 coupon from 1/4 to 1/2 per eligible marked
  pair, increasing its component-pair cost from 1 to 2.
- Add a half coupon between a whole A1 and the doubled-terminal mark of a
  whole B10, adding 1/2 per component pair.

The first supports the [order17 D118 refinement](fresh_order17_e6_half_coupon_transfer_20261018.md).
The second is included explicitly in this counting rule; this report does
not assert the separate ideal order11 theorem that may use it. All other
coupons and the proved depth11 zero rule and finite D16/D17 exchange are
retained from the [finite closure](fresh_ambient64_finite_closure_20261018.md).
The [source](../scripts/fresh_e6_b10_half_coupon_forest_budget_20261018.py)
and [saved certificate](fresh_e6_b10_half_coupon_forest_budget_20261018.json)
use standard-library exact arithmetic. This is a finite counting theorem,
not an ambient dimension bound.

## Complete finite-type reduction

Work in doubled costs. Same-component costs are unchanged: A16/D16 have
232, B16 has 208, A17 has 228, B17 has 210 and D17 has 254. For ranks at most 15
the conservative same-component bound is r(r−1). For A/B/D of rank r≥18,
start from 16r−68 for A/B, 16r−64 for D18/D19 and 16r−66 for D≥20. The new
whole-component deep distance8 zero rule subtracts twice the number of
eligible pairs, respectively

    max(0,r−30), max(0,r−31), max(0,r−32).

These are the same proved all-rank formulas used in the frozen finite
closure; their geometry is not reproved here.

Keep A17 separately because of its cross-component coupons. Every other
finite component is dominated by A_r for r≤9, except that E6 is retained
at rank6, and by D_r for r≥10. The inherited finite-type classification
includes A/B/D, E6/E7/E8/F4/H3/H4 and arbitrary dihedral labels. The source
checks same-component dominance and every one of 38,025 type/partner
comparisons at the enlarged coupon weights.

The changed B10 interface is explicit. With A1, doubled cross cost for
D10 is 2 and for B10 is 1; with A2 the costs are 4 and 2. Thus replacing B10
by D10 still dominates both its new A1 half and old A2 coupon. Other
partners give no additional B10 cost. E6 continues to dominate other
rank6 components, now with doubled E6/D11 cross cost 4. No old B10 coupon
is added again after replacing it by D10.

## Complete forest domain and exact certificates

Separate all D10 through D16 components and all A17 components. Their
multiplicities form 753 complete count vectors of total rank at most 64.
All remaining components have ranks1–9 or at least17, with E6 as the
rank6 representative and D_r as the representative for r≥17.

For counts d10,...,d16,a17, the five-coordinate interaction key is

    (2d10+3(d11+...+d16)+2a17, d11+d12, d13, d11, d14).

There are 348 keys. Let the coordinates be (b,u,a,d,e). The residual item
of rank r has its doubled same-component cost plus these cross terms:

- rb when r≤5;
- an additional 2u+e for r=1;
- an additional 4d+2e for r=2;
- an additional 2a for r=4;
- an additional 4d for r=6.

The last term is the upgraded E6/D11 cost. The new B10 coupon disappears
only after its explicit D10 domination check. These formulas also retain
the A1/D11 and A1/D12 unit upgrades, A4/D13 neighbor coupon, A1/A2–D14
neighbor halves, A17 coupons and all broad small-component coupons.
The source compares the compressed formula with literal typed component
interactions 42,921 times.

For each key, generation uses item-first unbounded dynamic programming.
The saved verifier does not rerun that optimization. For every residual
rank n it checks all last-item upper inequalities

    V[n] ≥ V[n−r]+w[r],

and one saved attaining predecessor with equality. Induction on n proves
both the upper bound for every residual forest and attainment by a typed
forest. There are 597,168 inequalities and 22,272 predecessors.

Each of the 753 large-component vectors is then combined with the exact
residual value at its remaining rank. Its attaining typed forest is
reconstructed and its full same-component plus cross-component cost is
recounted literally. This covers every rank64 forest after domination.
The largest doubled value is 928, hence cost 464, and four D16 components
attain it. They carry no cross coupons, so this attainment is unchanged
by either addition. It is attainment in the formal finite-forest domain,
not an assertion that such a vertex occurs in an actual chamber.

## Geometric and global scope

The baseline finite rules are already proved. Enlarging a finite coupon
cannot invalidate any old face-payment lower bound. The new E6/D11 half
is exactly the coupon required by the order17 refinement. Thus the finite
budget 464 remains available when that refinement is composed with the
baseline rule, including the optional A1/B10 half.

The new order11 ideal application and the complete ideal curvature
comparison are separate obligations. No conclusion about dimension63,
physical realizability, or a new example follows from this budget alone.
The certificate does not repeat the inherited geometric classification,
complete overlap census or all-rank deep-corridor proof.

## Validation

Final source-bound generation dfea64 and separate saved verification
abf19c both terminated exit0/PASS. They passed all complete count-domain,
type-domination, Bellman and typed-attainment checks. The saved verification
uses only integer/Fraction arithmetic and saved Bellman tables; no global
or local LP is invoked, and it does not regenerate the optimizing tables.

The parent and two other reviewers each completed a full mathematical
source/report read. They checked typed coupon accounting and domination,
E6/A17 distinctions, complete large/residual decomposition, interaction
compression, Bellman induction and typed attainment. Those source reviews
do not duplicate the arithmetic enumeration, the inherited finite-type
classification or geometric payment censuses.

All ten direct source bindings, the Python syntax and all four local
report links pass. Source, certificate and report are final and frozen;
there are no pending computations.
