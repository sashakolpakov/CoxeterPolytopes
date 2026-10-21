# A contextual payment removes the rank13 finite obstruction

Two changes lower the actual costs of both A2/A13 and A2/D13 component
pairs from5 to39/8. Their replacements preserve every triangle payment2
and quadrilateral payment1. The complete rank59 finite count includes a
new payment conditioned on other whole components at the same vertex.
Both formal forests10A2+3A13 and10A2+3D13 fall394 to1561/4=390.25,
below curvature783/2=391.5. Every representative forest containing a
rank13 component costs at most1563/4=390.75.

The overall exact finite maximum remains795/2=397.5, with the same three
actual typed attainers D16+D17+D26, A17+D16+D26 and2D16+D27.
There is no new dimension bound or constructed finite-volume polytope.

The [producer](../scripts/fresh_rank13_contextual_coupon_rebalance_20261018.py)
and [certificate](fresh_rank13_contextual_coupon_rebalance_20261018.json)
start from the [completed59 rule](fresh_depth6_dimension59_20261018.md)
with the [finite B14/D14 refinements](fresh_finite_rank14_existing_payment_transfer_20261018.md).
The earlier frozen packages remain unchanged.

## Exact new rule

Put delta=1/16. All conditions below refer to whole finite components
at the actual corner, not embedded subdiagrams.

* Between A2 and D13, reduce the fork-center weight1/2 to7/16.
  Either A2 mark is allowed; the D13 center is canonical node10.
  The four A2/short-fork-leaf unit roles remain units.
* Between A2 and A13, reduce the unique-middle weight1/2 to7/16.
  Either A2 mark is allowed; the A13 middle is canonical node6.
  The four A2/A13-endpoint unit roles remain units.
* Between A1 and A13, increase the unique-middle weight1/2 to9/16.
  Its other roles are unchanged.
* Between distinct A2 components, add1/8 for the marked pair if the
  same finite diagram contains a whole A10 component and at least two
  whole A1 components. All these additional components are unmarked.
  This is a single payment determined by their presence; it is not
  multiplied by the number of possible A10 or A1 witnesses.

All finite same-component and ideal weights stay unchanged. In particular,
actual D15 remains201/2. The rule is intrinsic under all component and
diagram automorphisms. It does not require choosing an orientation or
distinguishing one of several identical witness components.

## The three geometric interfaces that use a decreased coupon

In the [universal cycle14 transfer](fresh_cycle14_ideal5_transfer_20261017.md),
the only remaining unpaid equality is the sealed14-cycle with two opposite
simple leaves and an isolated spectator attached to one new leaf. Its two
adjacent finite corners have whole A2/A13-middle and A1/A13-middle roles.
They now pay7/16+9/16=1. The proof's all-label profile analysis, strict
spectator-energy hypotheses, height-first equality and arbitrary-support
orthogonality remain unchanged. The source reconstructs the equality Gram,
signature(15,1,1), positive-definite white block and four white-plus-black
edge diagrams, and checks the exact whole marked roles. It does not repeat
the predecessor's energy-spectrum classification.

In the [original order16 transfer](fresh_order16_finite_cross_payments_20261017.md),
only literal case384 uses the A2/D13-center role. Its two adjacent finite
corners each lose1/16. The opposite finite corner has exactly the displayed
whole components2A1+2A2+A10, with the marks in the two A2 components.
It gains1/8, so the four-corner total remains1.

The full18-root core has signature(17,1). Its actual affine16-node factor
has a strictly positive kernel, forcing every extra supporting white
orthogonal to that factor. The only possible remaining core contacts are
at roots3 and13. Both contacts supply an independently paid path of length2.
A single nonzero contact of magnitude t>=1 has projection square4t^2>=4,
exceeding norm2, since the inverse contact block is(4,-4;-4,4). Thus in the
remaining branch every additional supporting white is orthogonal to the
whole core. Neither the marked A2 factors nor the unmarked A1/A10 factors
can be enlarged or joined by a hidden spectator tail. The context condition
therefore holds in the full actual diagram. The source freshly reconstructs
the complete old nine-case order16 seal and checks the new corner roles.

Finally, the [three-singleton D15 theorem](fresh_finite_d15_three_singleton_half_transfer_20261018.md)
uses the center coupon in cases1147 and1269, paired with the selected
D15(1,12) half. Both cases also have the same opposite whole2A1+2A2+A10
corner. The old whole-D15 proof seals every other core component: allowed
remaining inverse diagonals are4,4,14, with strictly positive same-end
inverse blocks; a contact at both ends supplies the independent short8
route. The new coupon therefore applies. Each face now receives
1/2+7/16+1/8=17/16. No restoration of D15(1,12) is necessary.
The source rechecks the Lorentz signature, white positivity, remaining
inverse diagonals and all marked roles against both frozen seals.

## Complete output inventory and simultaneous composition

The [current coupon-output audit](fresh_small_d16_coupon_and_ideal17_rebalance_20261018.md)
divides the complete face-payment proof into universal short-path routes,
original sealed literals, later selected transfers and finite exchanges.
We refine that inventory by exact marked role, since the broad budget
allowances alone do not specify which geometric coupons are needed.

* The universal cycle/feature-free tree routes use A1/A15, A1/A17,
  A1/A11 or A2/A11, apart from the cycle14 middle pair just repaired.
  Whole connected alternatives and independently retained short paths
  are unaffected.
* The universal nine-pair affine-D13 theorem has exactly A5/D11,
  A2/D11, A4/D13 and A1/D13 cross outputs. Its description's separate
  lists of possible small and large factors must not be combined into
  an additional A2/D13 output. The affine-D14 theorem uses D12/D14.
  Neither consumes either decreased role.
* Original order16/17 sealed literal outputs are covered by the complete
  current216-configuration reconstruction. The only changed record is
  order16 case384, with two center debits and the contextual credit.
* Later ideal transfers use the separately sealed A4/D13, D5/D12,
  A4/D12, A5/D11, E6/D11, D5/D11 and A4/D11 roles, or smaller B/D
  factors. None is an A2/D13-center or A2/A13-middle role.
* The finite A15 transfer upgrades the A2/A13 **endpoint**, and the
  sparse rank15 exchange upgrades the A2/D13 **short fork leaves**.
  Those units remain unchanged. The earlier finite D15 alternative
  also uses a short fork leaf, not its center. The later three-singleton
  D15 center uses are the two seals repaired above. The finite14 and
  rank16/17 transfers use A4/D13 or their named same-component weights.

This is an audit of all outputs of the inherited exhaustive proof, not
a claim that every actual quadrilateral occurs in the216 saved cores.
Replay that proof with the same actual-corner priorities. Every unchanged
output retains its payment; the three changed interfaces now have the
replacement totals proved above. Short8-to-short7-to-short6 routes retain
their complete protected alternatives, including the repaired cycle14
case, and do not rely on a decreased coupon without this replacement.
Triangle payments use unchanged same-component or ideal elementary
weights, so retain2. Quadrilateral totals retain1.

The subsequent coupon-free B14 and C14 ideal deletions are compatible:
their first-entry projection exclusions use independent short8 escapes
and no new finite coupon. Conversely the present three seals use finite
roles or an independent short2 escape. Neither depends on their removed
ideal distance9 pairs.

As an exact supplementary interface check, all216 current literal
configurations are reconstructed from their Grams and actual marks,
starting from the finite14 packet's cornerwise lower payments. Exactly
two A2/D13-center debits and one contextual credit occur, all in case384.
Every corner lower bound is nonnegative and every total remains at least1.
The universal cycle equality and the later D15 cases are checked separately
as above; they are not silently assumed part of those216 literals.

## Complete finite count with the contextual charge

The two reduced component-pair costs are

    4*1 +2*(7/16) =39/8.

A1/A13 becomes25/16, while A1/D13 stays3/2. If a finite diagram has k
whole A2 components, a whole A10 and at least two A1, its new contextual
total is

    4 * binomial(k,2) * (1/8) = k*(k-1)/4.

There are four choices of the two marks for each unordered A2-component
pair. The charge is counted once for each marked pair regardless of how
many witness components exist. Without the required witnesses it is zero.

The rank59 certificate retains all961 fixed multiplicity vectors for
ranks9–17, with A17 fixed separately and D17 among the residual ranks.
There are960 interaction keys. All32,400 ordered actual type pairs satisfy
simultaneous representative domination, including both A2/A13 and A2/D13
reductions. Same-component bases and49,011 residual item identities are
checked exactly. The representative treats any rank10 component as an
A10 witness; this can overcount the contextual charge, so gives a safe
upper bound. It maps other rank2 types to A2 and preserves A1 counts.

For each key generate two exact tables: the ordinary residual table and
a table excluding rank2. The latter is necessary because the contextual
charge is quadratic in the number of A2 components. When rank10 is
present, enumerate the total A2 count k, reserve two A1 components, and
use the no-rank2 table for the remaining rank. Additional A1 components
are allowed in that table. This covers every diagram where the charge
is active. Compare all such options with the ordinary table, which
also bounds every diagram where the trigger is absent.

Fresh replay checks2,613,120 last-item upper inequalities and113,280
attaining predecessors across the1,920 tables, and all2,691 context
options. Each chosen row reconstructs its forest and exact contextual
charge. The maximum397.5 has actual typed attainers, so is exact rather
than a loose representative upper bound. Fourteen representative rows
remain above391.5, down from fifteen before this transfer. The maximum
among rows with rank13 is390.75; among rows with rank10 it is386.
Actual D15 remains100.5 while its representative allowance is101.

The certificate also includes a separate complete counting preflight:
increase the whole A1/D13 component-pair cost by1/2, as would result
from a new half at its unique long-arm neighbor of the center. Both
the overall397.5 maximum and the390.75/386 subdomain bounds remain
unchanged. This optional field asserts no geometric theorem for that
extra role. Any ideal improvement using it needs its own full proof.

All ideal weights are unchanged by the main rule, so the previous ten
negative terminal ideal comparisons are not repaired here. The other
finite obstructions also remain. No dimension58 conclusion follows.

## Verification scope

Final generation session97323/chunke0c336 and a separate fresh saved
replay session52827/chunk471153 both terminated exit0/PASS. The earlier
development replay revealed that both D15 seals also receive the new
context payment; the final rule therefore needs no D15 restoration.
The final source explicitly checks actual attaining forests and binds
the frozen B14/C14 interfaces.

A separate full source/report mathematical review reread the order16
projection seal, cycle14 equality, both D15 seals, complete coupon-output
inventory, exact affine-D13 output table and finite counting helpers.
It found no gap in arbitrary-support whole-component protection, the
once-per-marked-pair charge, noncircular payment composition, simultaneous
representative domination, the reserve-two-A1/no-A2-table enumeration,
or actual attainment. That review did not repeat arithmetic, native maps,
Bellman checks or convolutions. The root supplied the composition and
full source/report read; no external referee review is implied.

The producer rechecks the stated local Grams/seals, current literal
interfaces and all complete new forest certificates. It does not rerun
the inherited all-label census, old spectator-energy spectrum, cusp
convolutions, or construction searches. No external refereeing or
proof-assistant formalization is asserted.
