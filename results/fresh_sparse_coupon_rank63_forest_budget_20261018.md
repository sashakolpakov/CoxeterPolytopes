# All new sparse coupons fit finite rank63 maximum449

For the explicit finite rule below, every rank63 forest has cost at most
**449**, and the formal forest2D17+D29 has actual cost449. This is below
the finite curvature899/2 by1/2. The count includes every preceding
coupon and all nine new typed component-pair increments. It is a finite
counting theorem; simultaneous geometric face coverage and the global
ideal comparison are separate obligations.

The [source](../scripts/fresh_sparse_coupon_rank63_forest_budget_20261018.py)
and [certificate](fresh_sparse_coupon_rank63_forest_budget_20261018.json)
bind the [sparse finite transfer](fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018.md),
the [B12 transfer](fresh_order12_b_center_terminal_transfer_20261018.md),
the [D13 transfer](fresh_order13_d_center_pair_half_transfer_20261018.md),
and the [previous complete finite budget](fresh_e6_b10_half_coupon_forest_budget_20261018.md).
The typed marked roles and exact cost increments are compared directly
with those local certificates before the forest count runs.

## Actual rule and the conservative representative allowance

Same-component costs are A16/D16=114, B16=104, A17=114, B17=105 and
D17=125. All other preceding costs remain, including deep11 finite8
zeros and the large-D fork-leaf halves. Every old cross-component coupon
is preserved. The new increments are:

|Actual component pair|Added total cost|Marked roles|
|---|---:|---|
|A4,A12|4|endpoint/endpoint unit|
|A4,A13|4|endpoint/second-node unit|
|D5,A12|2|long-endpoint/second-node unit|
|A2,A13|2|either A2 node/A13 endpoint half|
|A2,A14|2|either A2 node/A14 second node half|
|E6,A11|4|E6 length2-arm endpoint/A11 second node unit|
|A2,B9|1|either A2 node/doubled terminal half|
|A2,D10|1|either A2 node/fork center half|
|A3,D10|1/2|unique A3 middle/fork center half|

The last two coupons are distinct from the old small-component/D10
short-leaf halves. The new A2/A13 coupon is distinct from its old
middle-node coupon. Thus the totals above are genuine increments, without
overlap or suppression of an existing payment.

For a uniform upper bound, replace both types of every component pair
simultaneously. Use E6 at rank6, B9 at rank9, retain A17 separately, and
otherwise use A_r for ranks at most8 and D_r for ranks at least10.
Use an **upper allowance4 for E6/D11**, although its actual coupon still
costs2. This absorbs the new actual E6/A11 cost4 when A11 becomes D11.
It is not an added E6/D11 payment rule.

For every pair of finite types a,b the source verifies

    upper_cross(R(a),R(b)) >= actual_cross(a,b),
    base(R(a)) >= base(a).

This is deliberately a simultaneous two-end comparison. A5 alone need
not dominate the new D5/A12 coupon against an unchanged A12. After the
second component is also replaced by D12, its existing small-component
allowance absorbs that cost. Summing these pointwise inequalities over
all components and all unordered component pairs proves a bound for
every finite forest, including repeated types. No sequential one-end
domination is assumed.

All192 finite type representatives through rank63 and all36,864 ordered
type pairs are checked. Arbitrary dihedral labels use the inherited
uniform rank2 representative. Finite/affine classification and the old
same-component formulas remain explicit inherited inputs.

## Complete count domain and interaction key

Enumerate every multiplicity vector of B9,D10,...,D16,A17 whose rank sum
is at most63. There are1,329 vectors. All other components have ranks1
through8, represented by A except E6, or ranks at least17 represented
by D. In particular A17 and D17 are different choices with different
coupons; B9 is not silently folded into the additive small-rank part.

Write b9,d10,...,d16,a17 for the large multiplicities and set

    b = 2*d10 +3*(d11+...+d16) +2*a17.

The exact doubled residual-item interaction key is

    (b, d11+d12, d13, d11, d14, b9, d10).

It has993 distinct values. For residual rank r, start with the doubled
same-component cost. Add r*b for r<=5, and additionally

    r=1: 2*(d11+d12)+d14;
    r=2: 4*d11+2*d14+2*b9+2*d10;
    r=3: d10;
    r=4: 2*d13;
    r=6: 8*d11.

The rank6 coefficient8*d11 is the doubled conservative E6/D11 allowance4.
The rank2 B9 and D10 terms are the two new component-pair increments1;
the rank3 D10 term is the new middle-node increment1/2. All existing
special A1/A2/A4 and D14 coupons are retained. Large components have no
mutual cross coupon, and residual components have no mutual cross coupon,
so these additive items cover every interaction. The source independently
compares each item with the literal typed upper-cross formula for every
large vector.

## Exact upper certificates and actual attainment

Item-first dynamic programming produces each residual table through
rank63. Saved verification instead checks every possible last-item upper
inequality and an attaining predecessor at every rank. Induction on rank
therefore proves both the complete upper bound and attainment in the
representative upper-cost problem; no sampled cutoff is used.

The993 tables contain1,592,772 upper inequalities and62,559 attaining
predecessors. There are73,095 literal interaction identities and1,329
global alternatives. Every alternative reconstructs a typed forest and
recounts both its upper cost and its actual cost, asserting actual<=upper.
The two costs need not agree on a forest containing E6 and D11, and such
a forest is not claimed to attain its conservative allowance physically
or even under the actual counting rule.

The conservative maximum is449. The explicit witness2D17+D29 has no
cross coupon, and its actual costs are

    125+125+199 =449.

The D29 value is the inherited8*29−33=199; no deep11 distance8 pair exists
in D29. Thus actual and conservative costs agree on this witness, proving
the exact actual-rule maximum449 as well as its upper bound. Formal
attainment makes no claim that this forest is realized by a chamber.

At ambient63, the finite curvature is

    (5*sum_(j=0)^62 min(j,62−j) −63*62)/2 =899/2.

Consequently the complete finite comparison has slack1/2, with every
new coupon included. This report supplies no ideal margin, complete
simultaneous face-payment theorem, dimension62 conclusion, or example.

## Verification status

Final source-bound generation (session81636, terminal chunkb7b486) and
separate fresh saved-certificate replay (session88456, terminal
chunk6c791d) both exited0/PASS. Both bind the final frozen D13 source,
data and report. The source imports only the standard library. It retains
separate actual and conservative cost functions and checks the complete
coupon-role records from the local inputs. No numerical optimizer is used.

The root reviewer and a separate peer each completed a full source/report
mathematical read. Their scopes include the nine typed coupon guards,
simultaneous two-end type domination, distinct B9/E6/A17 representatives,
the conservative E6/D11 allowance, complete multiplicity/residual split,
interaction compression, Bellman upper-bound and attainment argument,
and the actual-cost forest449. Neither review independently reran the
1,329-vector enumeration or Bellman arithmetic. The generation and fresh
saved-certificate replay above supply that exact computation; the local
geometric classifications and transfers remain explicitly bound inputs.
