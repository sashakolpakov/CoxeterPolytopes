# The E6/D11 quarter coupon leaves every finite maximum through75 unchanged

Adding a quarter payment between either long-arm endpoint of a whole E6
component and either short fork leaf of a separate whole D11 component
leaves every finite-forest maximum through rank75 unchanged. In particular,
the maximum at ambient65 remains487 against curvature480, and at66 it
remains495 against495. The new local ideal transfer is therefore affordable
within the existing finite budget, but it does not remove the finite
obstruction at65.

The [source](../scripts/fresh_e6_d11_quarter_coupon_forest_budget_20261017.py)
and [certificate](fresh_e6_d11_quarter_coupon_forest_budget_20261017.json)
prove this conditional counting statement. The geometric justification of
the coupon is the separate
[order17 D transfer](fresh_order17_affine_d_three_quarter_transfer_20261017.md).
There are two eligible roots of each type, so the added total per component
pair is2*2*(1/4)=1. No old coupon has E6 as its small component.

## Complete finite counting

Start with the frozen
[full A2/D11 finite budget](fresh_a2_d11_full_coupon_forest_budget_20261017.md).
Retain its same-component payments and every existing cross-component
coupon. The only change is the additional E6/D11 cost just stated.

For type domination, represent every rank6 component by E6, whose
same-component upper cost is unchanged and whose new interaction is
nonnegative. All other representatives stay the same. The source checks
every actual-type/partner comparison, including exceptional types and
the full I2(m) family represented by its inherited uniform cost.

The1,608 complete large-component vectors and357 interaction keys remain
valid. Large components are D10 through D16 and coupon-bearing A17.
Residual components have ranks1 through9 and17 onward. There is no coupon
joining two residual components; the new E6/D11 interaction has its D11
member in the large vector. For each fixed vector, add exactly2*d11 to
the doubled residual cost of rank6. Other residual costs remain unchanged.

Generation uses item-first unbounded knapsack tables. The verifier checks
every last-item upper inequality and an attaining predecessor for every
entry. It also checks the compressed interaction formula against literal
component-pair costs, reconstructs actual typed forests, and checks every
global alternative at every rank through75. The resulting certificate has:

|Check|Count|
|---|---:|
|Large component vectors|1608|
|Residual interaction tables|357|
|Literal interaction formula checks|109344|
|Bellman upper inequalities|860013|
|Attaining Bellman predecessors|26775|
|Global forest alternatives|22504|
|Typed global attainers|76|

Every resulting maximum agrees exactly with the frozen predecessor.
This is an upper-and-attainment proof, not a search over selected forests.
It reuses the predecessor's complete large-vector domain and payment
formulas; it does not independently re-prove the geometric payment rules.

## Verification

Generation79676 and separate saved replay40871 both terminated exit0/PASS.
A separate internal full source/report mathematical review, including the
imported base, cross-cost, representative, domain and residual-cost helpers,
found no unresolved step. It did not independently replay the arithmetic.
The replay reads the saved tables
and checks all upper and attaining relations rather than rerunning their
item-first generator. No new ambient theorem or example is asserted here.

    python3 scripts/fresh_e6_d11_quarter_coupon_forest_budget_20261017.py
    python3 scripts/fresh_e6_d11_quarter_coupon_forest_budget_20261017.py --verify
