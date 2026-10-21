# Sparse D16 coupons and a fractional order17 reversal

Starting from the [completed407.5 rule](fresh_d15_balanced_ambient60_frontier_20261018.md),
make two simultaneous changes.

* Remove the half-payment between a whole finite component of rank2,3,4
  or5 and a whole D16 marked at its fork center or a short fork leaf.
  The whole A1/D16 halves remain unchanged.
* At an ideal factor of order17, restore weight1/8 to the four affine-D16
  pairs `(2,15),(2,16),(10,13),(10,14)`. Lower the distinct whole-D5
  long-endpoint/whole-D12 short-leaf coupon from1 to7/8.

Every triangle still receives at least2 and every quadrilateral at least1.
The actual order17 family totals become A85, B111, C108, D229/2.
The whole D5/D12 component-pair cost falls from23/2 to45/4.
No finite same-component weight changes and no new coupon is added.
All statements concern finite-sided, finite-volume real hyperbolic Coxeter
polytopes under the same actual-corner hypotheses as the completed rule.

The complete formal rank60 finite maximum remains815/2, against curvature
405. However, the two smaller-component obstructions disappear:

|Formal forest|Previous cost|New cost|
|---|---:|---:|
|3D16+2A5+A2|1629/4|1413/4|
|4A1+4D5+3D12|406|403|

Only six representative alternatives now exceed405. Each has cost815/2
and consists of two components chosen from A17,D16,D17, together with
D26,D27 orD28 to make rank60. These have actual typed attainers. Thus
this theorem removes the smaller-component obstructions; it does not
exclude dimension60 or establish a new unrestricted bound59.

The [producer](../scripts/fresh_small_d16_coupon_and_ideal17_rebalance_20261018.py)
and [certificate](fresh_small_d16_coupon_and_ideal17_rebalance_20261018.json)
give the exact changed rule, inherited proof bindings, four reconstructed
radical seals, all216 current literal corners, complete new finite
certificates and all2,113 ideal row updates.

## Why the discarded D16 coupons are unnecessary

The broad small-component allowance entered through a counting theorem.
The [first broad budget](fresh_broad_finite_half_coupons_20261017.md)
enlarged the earlier A1/D16, A2/D14, A1/D15 and A2/D13 coupons to all
small components of rank at most4 against D12 through D16. The later
[rank5 budget](fresh_finite8_a17_d11_coupon_budget_20261017.md) enlarged
this further. Both separate the budget enlargement from the geometric
question of whether a whole-component payment is needed. In particular,
the enlargement alone is not a geometric requirement to retain every
added small/D16 pair.

The necessary proof here is an audit of the actual outputs of the
complete face-payment argument. Absence from a finite sample would not
suffice. The following is the complete division of its coupon-producing
interfaces.

|Interface|Cross-component outputs relevant here|
|---|---|
|Universal cycle/feature-free tree short-path replacements|Sealed A1/A13, A2/A13, A1/A15 or A1/A17 coupons; cycle12 also uses A1/A11 units and A2/A11 halves; other outputs have both marks in one component or a retained short path|
|Universal nine-pair affine-D13 transfer, factor order14|Whole A1,A2,A4 orA5 against D11 orD13 fork marks|
|Universal nine-pair affine-D14 transfer, factor order15|Whole small rank at most3 against D12 orD14 fork marks|
|Literal order16 and order17 original replacements|Whole A1/D15, A2/D13, A1/A15; whole A1/D16 and A2/D14, respectively|
|Later selected ideal transfers at orders14–17|Complete literal records, or their separately sealed A4/D13, D5/D12, A4/D12, A5/D11 and E6/D11 roles|
|Order11 and order12 universal spectator transfers|Small whole A components against B9,B10,B11,D10 orD11; the remaining outputs are short paths, whole finite components or retained small affine factors|
|Order13 center and center/leaf spectator transfers|Whole A2/A3 against D10, whole A1/D12-center, whole finite D12/D13, or independent retained ideal/short-path outputs|
|Later finite rank15–17 transfers|Whole same-component alternatives, or the explicit sparse A/D roles listed below; no new small/D16 cross coupon|
|Large finite-D protected-leaf transfer|Whole small rank at most5 against D10 short leaves|

For the first row, the [finite8 cycle theorem](fresh_finite8_cycle_replacement_20261017.md),
[ideal5 local theorem](fresh_finite8_ideal5_local_replacement_20261017.md)
and [cycle14 theorem](fresh_cycle14_ideal5_transfer_20261017.md) classify
all supporting contacts. Their equality coupons are the stated sealed
A1/A2 versus path roles. The [cycle12 theorem](fresh_cycle12_ideal5_transfer_20261017.md)
adds the stated A1/A11 unit and A2/A11 halves; the
[cycle13 theorem](fresh_cycle13_ideal5_transfer_20261017.md) adds no coupon.
A finite D16 may occur as a single whole marked component;
that is a same-component payment, which is untouched. The shortened
finite8-to-short7 and finite7-to-short6 routes inherit these independent
outputs. The current whole-A8 alternative is also a same-component unit.

The two universal fork-transfer rows follow directly from the full
[order14](fresh_affine_d13_long_pair_transfer_20261017.md) and
[order15](fresh_affine_d14_long_pair_transfer_20261017.md) statements.
Their finite equalities include arbitrary supporting spectators and
prove that the named factors remain whole. The large ranks are bounded
by14, not by an unspecified D11-through-D16 budget allowance.

The original [order17 literal replacement](fresh_ideal7_finite_cross_payments_20261017.md)
needs only whole A1/D16 and A2/D14; its positive affine kernel and
full-core radical seal the entire actual supporting diagram, or supply
an independent length-two unit. The corresponding
[order16 replacement](fresh_order16_finite_cross_payments_20261017.md)
has the lower ranks in the table. These are the literal transfers used
by the [complete finite8 theorem](fresh_finite8_ideal5_combined_payments_20261017.md),
whose full low/disjoint/overlap domain is retained here. Low-core,
protected/adaptive and ordinary selected routes give same-component
payments. Cross-component replacements enter through the explicitly
listed universal and sealed literal cases.

The later order16/order17 partial-weight transfers have explicit whole
D5/D12, A4/D13, A4/D12 orA5/D11 outputs. The current refinements add
E6/D11, D5/D11, D5/D12, A4/D11 and A5/D11 roles. The order13 sparse
terminal transfer adds A4/B9. The current universal B/C11 and B/C12
proofs use only the small B/D factors listed above. Their energy-spectrum
boundary cases can give larger single whole components, but never turn
them into a small/D16 coupon. Those single-component outputs retain
the same weights as before this change.

The subsequent finite exchanges use existing small/D13 or small/D14
halves, A2/A13, A1/A13, A2/A14, A4/A12, D5/A11, and D5/D12 roles, or
restored same-component A17/D17 weights. The latest D15 theorem uses
whole A2/D13-center or same-component D16 halves. Its D16 pairs
`(1,12),(1,13)` are unrelated to the discarded cross-component coupons.
The complete current proof bindings retain these whole-component
statements; none is inferred from a truncated component.

As a second, exact check of this inventory, the producer reconstructs
every finite and ideal corner in all216 current literal records from
its Gram matrix and supporting whites. There are exactly ten corners
with marks in separate components, one of type D16. In all ten the
other component is A1. Six are in the original order17 cross-payment
package and four in the old order17 partial-weight package. No rank2–5
coupon is debited. Arbitrary extra whites cannot enlarge these A1
components: their inherited radical/equality arguments seal the whole
core before claiming the coupon. The other later explicitly typed
coupons above are likewise whole-component conclusions.

This is why the removal is universal. Replay the existing full-domain
proof with the same actual-corner priorities and the stated output
inventory; every payment it uses survives. It is not an argument that
a quadrilateral containing a deleted coupon must itself be one of the
216 saved cores. Such a face is covered by the original exhaustive
proof, whose actual paying outputs have just been retained.

## Four ideal eighths and two finite seven-eighths

Use the complete frozen [D114 theorem](fresh_order17_d5_d12_unit_coupon_zero_transfer_20261018.md).
Its 140 eligible native patterns and488 successful induced maps locate
all affected images in four whole-factor contexts:1978,2371,2535,2539,
with32 affected maps. Each has a complete19-root core, a whole17-root
affine-D16 component and a distinct whole finite D5/D12 coupon corner.

The positive affine kernel forces every additional supporting white
orthogonal to the affine component. Pairing with the actual full-core
radical then either gives the independent opposite-side path of length2
or forces all remaining core contacts to vanish. In the latter branch
all displayed finite factors remain whole; an indirect supporting tail
cannot join them. The producer recomputes all four old seal certificates
and compares them exactly with the frozen proofs. It does not rerun the
historical complete embedding census.

The selected ideal pair now contributes1/8 and the distinct D5/D12
corner contributes7/8. Their sum is1 without spending any same-component
finite payment. Intrinsically, the D5 mark is its unique long-arm
endpoint, and the D12 mark is either short fork leaf. There is one
eligible D5 node and two eligible D12 nodes, so the finite component-pair
change is exactly

    1 * 2 * (7/8 − 1) = −1/4.

These two marked roles `(0,10),(0,11)` are disjoint from the reverted
center role `(0,9)` and from the six finite-D17 compensation roles
`{3,4} × {9,10,11}`. The latter still receive full units. Hence the
current D5/D12 total becomes23/2−1/4=45/4. The four ideal pairs gain
4/8=1/2, taking D114 to229/2. All other ideal-family costs stay fixed.

The complete selected/nonselected scope of D114 is inherited: the four
changed pairs have distance11, so basic, short8 and elementary triangle
weights are untouched. Earlier small-ideal actual-corner priorities,
current short8→short7→short6 replacements and the whole-A8 branch remain
independent of these roles. In the selected order17 branch, the four
seals now use1/8+7/8; all other images retain their proof. Together with
the preceding coupon-output inventory this gives a noncircular current
whole-face composition. Every other finite reduction and the current
rank16/rank17 balances keep exactly their old weights.

The new216-corner audit begins with the final407.5 literal totals.
Exactly four ideal credits1/8 and four finite coupon debits1/8 occur;
they cancel in their four records, whose totals stay1. The D16 coupon
removal has zero occurrences. All216 totals remain at least1. Triangles
use retained same-component/ideal elementary payments, so retain2.

## Complete changed finite count

Keep the actual finite base costs from the407.5 rule, including
A15=101, B15=97, D15=201/2. The rank15 representative allowance remains
101; actual D15 does not dominate A15. Every other actual base is
unchanged.

The producer checks all33,489 ordered type pairs for simultaneous
representative domination after both coupon decreases. Each actual
small rank2–5/D16 total loses3r/2; its representative loses the same
amount. Other types mapped to D16 have no such cross coupon. For
D5/D12 the actual and representative allowances both lose1/4. Whole
A1/D16 still costs3/2. No extra cost is concealed in a representative.

The complete fixed-component domain still contains1,044 multiplicity
vectors for ranks9–17. The old seven-coordinate interaction key is
extended by the number of D16 representatives; there are now1,042 keys.
This is necessary because rank2–5 residual costs now distinguish D16
from the other formerly grouped large components.

In eighth units the residual item change is

    −12*r*c16,       2 <= r <= 5,
    −2*c12,          additionally when r=5.

All residual–residual and fixed–fixed interactions remain zero. Every
residual item is checked against its explicit base plus all fixed cross
costs. Generation computes new item-first Bellman tables. Fresh replay
checks every last-item upper inequality and an attaining predecessor;
it does not optimize again. There are1,502,564 upper inequalities,
62,520 attaining predecessors and54,288 item identities. Every fixed
vector reconstructs an attaining representative forest. The unchanged
actual forest2D17+D26 attains815/2, matching the complete upper.

For the first targeted obstruction, each of three D16 factors loses
3*(2*5+2)/2=18 cross payment, so1629/4 becomes1413/4. For the second,
twelve D5/D12 component pairs each lose1/4, taking406 to403. These are
actual typed cost calculations, not representative forests mistaken
for geometric vertices. Neither formal forest is asserted realizable.

## Complete positive ideal comparison

Independently enumerate the terminal rank59 partitions with at most
one part1–8 and all other parts at least9, and compare them bijectively
with the frozen2,113 rows. Their exact convolutions are inherited;
every new scalar row is computed exactly. The only normalized cost
increase is

    count(rank16)/34.

All2,113 margins remain strictly positive. The normalized minimum is
31/6000 at(9,9,9,9,9,14), and the absolute minimum is54 at(59).
The exact allowable total increase of the order17 envelope before a
row reaches zero is897/1360, attained at(9,9,9,16,16); the chosen1/2
is strictly smaller. This threshold is a scalar observation, not a
license to change any additional geometric payment.

The inherited low9 concentration changes only low-rank factors, whose
cost remains r/2. It decreases both normalized curvature and V. Since
every terminal margin is positive, the absolute bound54 extends to
all actual cusps. The finite bound is still405+5/2, so the same global
comparison gives

    (5/2)*v_f >= 4500 + 54*v_i.

The completed unrestricted ceiling60 and all-ideal ceiling59 remain
unchanged. The new work supplies a smaller obstruction list, not a
proof of unrestricted59, a construction, external review or formalization.

## Verification scope

Final source-bound generation42974/98e044 and separate fresh
replay17587/badde0 both terminated exit0/PASS for the four seals,
all216 literal records, all changed finite certificates and all2,113
ideal row updates. Root completed a full producer/report mathematical
composition read, including the broad-coupon and original finite8
provenance and D114 reports. That read checked the exposed D16 count,
simultaneous domination, complete fixed/residual domain, disjoint
fractional roles and ideal row update, without executing arithmetic or
repeating the historical census. A second peer completed a full
producer/report mathematical read and the inherited broad/finite8,
cycle12/cycle14, universal-fork and original order16/order17 provenance
checks. It checked the complete actual-output inventory and arbitrary-
white sealing, four disjoint fractional seals, literal deltas, typed
domination, the new interaction key, Bellman proof and ideal updates.
Its cycle12/cycle14 inventory precision corrections are incorporated
above. Neither peer reran arithmetic, maps, geometry, DP or convolutions.
No mathematical gap was identified within this stated scope.

All176 direct bindings,16 local report links and Python syntax pass.
Source, certificate and report are final/frozen. No historical geometric
census or convolution is repeated by this producer; their complete
frozen certificates remain explicit inputs. These are internal agent
mathematical reviews, not external refereeing or proof-assistant checks.
