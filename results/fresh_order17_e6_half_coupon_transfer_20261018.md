# Order17 affine D cost118 with a typed finite half coupon

The whole affine D16 factor (17 roots) can use total weight **118**, with
all earlier weights preserved except the four distance10 pairs
(3,15), (3,16), (9,13), (9,14). Their weights fall from 3/4 to 1/2.
The required finite change is equally specific: raise the coupon from 1/4
to 1/2 between a whole E6, marked at an endpoint of either length 2 arm,
and a whole D11, marked at a short fork leaf. There are four eligible
marked pairs per pair of components, so the additional finite cost is
**one per E6–D11 component pair**.

The [source](../scripts/fresh_order17_e6_half_coupon_transfer_20261018.py)
and [certificate](fresh_order17_e6_half_coupon_transfer_20261018.json)
refine the [frozen quarter-coupon theorem](fresh_order17_affine_d_three_quarter_transfer_20261017.md).
This is a local face-payment theorem. It does not by itself prove a finite
forest budget or an ambient dimension bound.

## Why the same complete geometric proof applies

The frozen theorem reconstructs all 140 eligible selected patterns and 96
affected maps. Its 64 proper completions are impossible, using complete
necessary contact control, projection obstructions and Lorentz index.
These obstructions do not depend on the amount of the finite coupon.
Only four whole-affine contexts remain: 1981,2374,2549,2551, each with 8 maps.

Each has a 19-root core of signature (17,1,1). The positive affine kernel
forces every further supporting white to be orthogonal to the 17-root
ideal component. The two remaining possible core contacts have opposite
coefficients in the actual zero radical relation. Consequently both
contacts vanish, or the retained length 2 path pays. In an otherwise
unpaid face every extra white is orthogonal to the entire core. The
finite E6 and D11 components at the alternative corner are therefore
whole components, with no indirect supporting tail.

The present certificate checks the saved radical and positive affine
kernel exactly, and reconstructs the two finite component trees directly
from the saved Gram matrix. Their fork-arm lengths are respectively
(1,2,2) and (1,1,8). The actual marks are an E6 length 2 endpoint and a D11
short leaf. In each of the four configurations, the old proof used exactly
one ideal 3/4 payment and one finite 1/4 coupon, at different corners.
The new proof uses 1/2+1/2=1 at the same corners. No other finite payment
is consumed. In particular, the displayed finite D16 or D17 units remain
explicitly unused.

The four changed ideal pairs form a complete automorphism orbit and all
have distance10. Every basic pair, short8 pair and retained triangle pair
stays a unit. Completeness of the old selected, low-core, disjoint,
immediate and later interfaces is inherited from the frozen theorem;
this small refinement does not replace that all-label classification.
Its geometric alternatives are independent of the numerical split of 1
between the two specified corners.

## Simultaneous scope

All other ideal weights and finite payments are retained. The new finite
half coupon dominates its old quarter value, so faces not using the
changed ideal orbit cannot lose payment. When a changed ideal occurrence
is needed, the complete preceding argument gives either a retained short
alternative or the sealed finite half, retaining the quadrilateral total 1.
The changed pairs have distance10, and the prior triangle argument remains
unchanged with total at least 2.

The finite coupon lies in components of ranks 6 and 11. It is therefore
unaffected by the depth11 finite-distance8 zero rule, or by the separate
finite D16/D17 two-unit exchange. No existing same-component D16/D17
payment is used. The only interaction with newer order11 rules is through
a short8 alternative in an inherited branch. Among the proposed B/C/D changes,
only the opposite-fork pairs of order11 D have distance8; they retain
weight 1/2. At a distinct alternative corner that half already combines
with the original retained order17 half to give 1. More generally the
independent complete order11 transfer supplies its own other half using
another order11 half, a retained finite payment or a typed finite coupon;
it never uses the changed order17 distance10 orbit. A simultaneous binder
must use that independent transfer, rather than assume an old short8 unit
is unchanged. The frozen cycle11 short5 halves can also occur in a
short8 branch: a distinct such half and the original order17 half suffice.
A deep finite8 zero alternative must instead be handled through its
independent short7 replacement theorem. Complete simultaneous branch
bookkeeping belongs to the global binder; this packet alone does not
certify that final composition or the pending order11 theorem.

The order17 family totals become A85/B111/C108/D118, giving envelope118/17.
The usual cusp factor-fiber multiplicity is unchanged. The added finite
component-pair cost must be included in a separate complete forest budget;
no global conclusion is asserted in this report.

## Verification scope

Generation a4fd87 and separate fresh replay f68e65 both terminated with
exit 0/PASS. The parent completed a full producer/report mathematical read,
including the frozen predecessor's complete written argument. A separate
reviewer completed the full new source/report and predecessor-proof read,
passing the role counts, inherited proper exclusions, weight refinement and
explicit global-composition caveat. Neither reviewer repeated the geometric
census or arithmetic execution. The implementation uses only the standard
library and imports no producer arithmetic. It
checks all 136 saved pair weights, the changed orbit, exact radical/kernel
identities, four literal coupon-role trees and four payment totals. It
binds the frozen complete geometric theorem and current dimension64
package. It does not rerun the 140-pattern/96-map census or independently
prove its arbitrary-support sealing and completeness arguments.

All seven direct source bindings, the Python syntax and all three local
links pass. Source, data and report are final and frozen.
