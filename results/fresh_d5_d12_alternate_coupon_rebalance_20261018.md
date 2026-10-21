# Alternate whole corners for the D5/D12 coupon

Starting from the frozen [rank13 contextual rule](fresh_rank13_contextual_coupon_rebalance_20261018.md)
and the now proved [order14 D-near/fork rule](fresh_order14_d_near_fork_half_transfer_20261018.md),
make the following simultaneous changes. A component in this statement is
a whole connected component of the actual finite corner diagram.

* Lower the payment between the unique long-arm endpoint of D5 and
  either short fork leaf of D12 from7/8 to1/2.
* Add3/8 between the unique long-arm endpoint of D5 and either endpoint
  of a whole A12.
* Add3/8 between that D5 endpoint and a whole A1.

The added amounts are increments, not replacements for other payments.
All old broad halves, finite same-component weights, ideal weights and
other sparse coupons remain unchanged. Every triangle still receives
at least2 and every quadrilateral at least1, under the inherited
finite-volume actual-corner hypotheses. The four affected sealed faces
each receive1/8+1/2+3/8=1 at three distinct corners.

The complete formal rank59 finite maximum remains **795/2=397.5**, above
the required783/2=391.5. The number of representative large-component
alternatives exceeding391.5 falls from14 to12. This is a local payment
and counting improvement; the completed unrestricted dimension bound
remains59. No new chamber or new dimension bound is claimed.

|Formal rank59 forest|Previous actual cost|New actual cost|
|---|---:|---:|
|3A1+4D5+3D12|791/2|391|
|8A1+3D5+3D12|1557/4|783/2|

The second forest shows why the new A1/D5 cost must be counted. The
statement concerns formal forest costs, not realizability of these
forests by hyperbolic Coxeter polytopes.

The [producer](../scripts/fresh_d5_d12_alternate_coupon_rebalance_20261018.py)
and [exact certificate](fresh_d5_d12_alternate_coupon_rebalance_20261018.json)
contain the reconstructed full-core proofs, marked alternatives, all216
current literal interfaces and the complete finite upper/attainment
certificate. No historical overlap census or ideal convolution is rerun.

## Four distinct whole-corner alternatives

The [original order17 transfer](fresh_order17_d5_d12_unit_coupon_zero_transfer_20261018.md)
has a complete140-pattern selected domain,488 successful native maps,
and32 affected maps in contexts1978,2371,2535,2539. Each affected map
already contains the entire17-root affine factor. Its current ideal
weight is1/8 under the [sparse-D16 and order17 rebalance](fresh_small_d16_coupon_and_ideal17_rebalance_20261018.md).
The latter used a distinct D5-long-end/D12-short-leaf payment7/8.

The new producer reconstructs the four full19-root Gram certificates
using the inherited native routine, compares every reconstructed seal
with the frozen rule, and independently extracts every marked finite
component at all four corners. It checks the full signature(17,1,1),
positive-definite supporting whites and all four strict edge principals.
Here are the three paying corners; numbers are literal core labels.

|Context|Ideal corner, weight1/8|D5/D12 corner, weight1/2|New alternate corner, increment3/8|
|---|---|---|---|
|1978|(1,15)|(1,13)|(3,13): A12 endpoint3 / D5 long endpoint13|
|2371|(2,15)|(2,13)|(3,13): A1 node3 / D5 long endpoint13|
|2535|(5,11)|(3,11)|(3,13): D5 long endpoint3 / A12 endpoint13|
|2539|(5,12)|(3,12)|(3,13): D5 long endpoint3 / A1 node13|

Each row uses three different corners. The fourth corner contributes
zero to this lower bound. No whole D17 or D16 same-component payment
at that fourth corner is used. In particular, no previously deleted
long finite pair is silently restored.

The component roles are intrinsic. D5 has fork-arm lengths1,1,2,
so its unique long-arm endpoint is canonical node0. The two endpoints
of A12 have canonical images0 and11. The two D12 short leaves have
images10 and11. The producer checks all canonical images, rather than
choosing one graph isomorphism. This makes every added and reduced
marked role invariant under the corresponding component automorphisms.

## Arbitrary supporting whites

The entire affine component consists of the core roots other than3
and13. Its positive kernel and nonpositive actual wall products force
any additional supporting white root to be orthogonal to all seventeen
affine roots. The full19-root core has a negative direction and a
one-dimensional radical. Its radical is an actual zero vector relation:
the orthogonal complement of a subspace containing a negative direction
is positive definite. The radical coefficients at3 and13 have opposite
nonzero signs.

Pairing an additional white with this zero relation forces its two
remaining contacts either both to vanish or both to be nonzero. In the
second case the path3–white–13 joins the adjacent black sides of the
actual corner(3,13) and has an independently retained length2 payment.
On the remaining branch every additional supporting
white is orthogonal to the full core. A longer bridge is also impossible,
because its first root touching the core would violate this conclusion.

Consequently all the displayed finite components, including the new
A12 and A1 alternatives, remain whole in the actual corner after
adjoining an arbitrary orthogonal positive-definite spectator forest.
The old seal applies to every corner, not only the previously paying
D5/D12 corner. This supplies the required arbitrary-support argument
for the new three-corner total. It assumes neither a bound on spectator
rank nor a simple ideal link.

## Complete payment-use inventory and current composition

The decreased roles return to their original broad half. Therefore any
earlier proof which used only that half remains valid. The complete
output audit in the frozen [sparse-D16 rule](fresh_small_d16_coupon_and_ideal17_rebalance_20261018.md)
and its [rank13 refinement](fresh_rank13_contextual_coupon_rebalance_20261018.md)
is retained here. The following separates all possible uses by the
actual marked role, rather than by a broad component-cost allowance.

|Geometric interface|Payment retained here|
|---|---|
|Universal short cycle/tree, low/disjoint and protected supporting routes|Their original broad halves, specified A-path coupons, retained short paths or whole same-component alternatives|
|Original order17 partial-weight transfer|Any D5 mark against a D12 center/short leaf only requires the old half|
|Four order17 long-end/short-leaf upgrades,1978/2371/2535/2539|The three distinct whole corners just proved, or their independent length2 alternative|
|Order17 center transfer and its later reversal|D5 long endpoint0 / D12 center9 remains1/2|
|Finite-D17 D5/D12 unit transfer|All six roles {3,4}×{9,10,11} remain full units|
|Other later ideal/finite sparse transfers|Their separately specified A4/D13, A4/D12, A5/D11, E6/D11, D5/D11, D5/A11, small B/D or A-path roles remain unchanged|
|Current rank13 contextual proof|Its A2/A13-middle, A2/D13-center, A1/A13-middle and context-dependent A2/A2 payments remain unchanged|
|Current order14 D-near/fork theorem|A4/D10-short-leaf, A1/D13-neighbor9, whole finite D14, retained D11 half or independent short/BC11 outputs; no D5/D12 output|

The original [partial transfer](fresh_d17_partial_weight_transfer_20261017.md)
explicitly permits an arbitrary marked node of a whole D5 against a D12
short leaf at weight1/2. The only later upgrade of the long-end/short-leaf
roles was the four-seal theorem above. The [finite-D17 unit theorem](fresh_finite_d17_d5_d12_unit_transfer_20261018.md)
instead uses the two D5 short leaves, canonical3/4, and is disjoint.
The long-end/center role0/9 is also disjoint and unchanged. Thus the
current component-pair cost is the broad15/2 plus the six half-unit
increments3, namely21/2; there is no residual charge from either
long-end/short-leaf upgrade.

Use the same actual-corner priorities as the current completed proof.
The independent small-ideal transfers and short8-to-short7-to-short6
routes retain all their protected outputs, including the whole-A8
alternative and the current depth6 finite7 exceptions. The current
B/C/D14 basic/near-fork refinements are resolved by their own independent
outputs. On the remaining selected order17 branch, the complete
140-pattern native map domain has exactly the four affected contexts
above. Their full-core proof now gives1/8+1/2+3/8. Unchanged selected
images retain their prior proof; low, disjoint, protected, immediate
and later routes retain the outputs listed in the inventory.

This is a replacement of the actual paying outputs of the exhaustive
face proof. It is not an assertion that every quadrilateral with a
D5/D12 corner belongs to one of four saved cores. No new arbitrary
supporting-white restriction is inferred from absence in a saved census.
Triangles retain their elementary/same-component payments and total2;
all old broad half payments survive, and the new coupons are positive.

As an exact simultaneous interface check, the producer starts with the
current rank13 literal lower bounds and aligns all216 records and all
marked corners with the frozen order14 D-near/fork interface. The latter
has no literal debit or credit. The new rule has four D5/D12 debits,
two A12 credits and four A1 credits. Four paired records stay at1.
The two additional A1 credits occur in the order16 partial-transfer
records[2372,0] and[2540,0], whose totals rise to11/8. All216 totals
remain at least1. This corroborates the complete output proof; it does
not replace it.

## Exact finite accounting and the complete rank59 certificate

In units of1/16 the actual component-pair changes are

* D5/D12:−12, because two marked pairs each lose3/8;
* D5/A12:+12, because two endpoints each gain3/8;
* A1/D5:+6, because one marked pair gains3/8.

The resulting actual component costs are respectively21/2,11/4,3/8.
The newly proved A1/D13-neighbor half is included: the actual total
A1/D13 cost is2. It was the optional complete-count branch in the
frozen rank13 certificate and now has its independent geometric proof.

The same representative type map is used as in that certificate. The
only changed representative costs are A5/D12:−12 and A1/A5:+6.
The D5/A12 increase is absorbed by the reduced A5/D12 allowance.
The source verifies simultaneous base and cross-cost domination for
all180 actual finite types through rank59, hence32,400 ordered type
pairs, including every added and reduced coupon. The nonlinear context
coupon is separately dominated by the same representative trigger:
any actual whole A10 supplies a fixed rank10 component, and A1/A2
counts are preserved by the representative map.

Fix every count of ranks9 through17 exactly as in the inherited
961-vector domain. Also fix every possible count k5 of A5 in the
representative forest, from0 to the remaining rank divided by5. There
are2,574 such combined vectors and2,572 distinct interaction keys.
The residual item ranks are1,2,3,4,6,7,8 and17 through59: rank5 is
excluded because its entire multiplicity is already fixed. The
fixed-rank17 A17 and residual-rank17 D17 roles remain distinct, as in
the inherited typed count.

All costs involving a fixed A5 are included explicitly. In particular,
each residual A1 gains6*k5 units. The source checks every residual
item cost against the complete fixed-component cross sum, for128,700
item identities. Residual components have no mutual pair cost after
rank5 has been fixed, apart from the already isolated contextual
A2/A2 interaction.

For each key, the producer makes two item-first unbounded tables: one
with all residual ranks and one excluding rank2 as well. The replay
does not optimize. It checks every possible last-item upper inequality
and one attaining predecessor for every rank0 through59. The5,144
tables give6,718,064 exact Bellman upper inequalities and303,496
nonzero-rank attaining predecessors.

The untriggered option uses the ordinary table. If a fixed rank10
component is present, the triggered option reserves two A1 components,
fixes the total number j of A2 components, and uses the no-rank2 table
for the remaining rank. Its contextual increment is4*j*(j−1) units.
The no-rank2 table may add further A1 components, so the reservation
does not constrain the actual A1 count to two. Every triggered forest
has exactly one such total j and is covered; every untriggered forest
is covered by the ordinary table. All6,455 options are checked.
The reconstructed maximizing forest is evaluated with the full
nonlinear rule, including its actual trigger, so an ordinary-table
attainer that happens to trigger the bonus cannot be overlooked.

Each combined fixed-vector row has an exact upper value and an
attaining representative forest. Taking the maximum over k5 gives the
complete961 large-vector alternatives. Their maximum is6360/16=795/2.
The actual typed forests D16+D17+D26, A17+D16+D26 and2D16+D27 all
attain it. Twelve alternatives remain above6264/16=783/2. There is
therefore no dimension58 conclusion from this certificate.

## Verification and limits

Generation session15966 terminated exit0/PASS (terminal chunkb28479;
saved log4503c1). The initial runs stopped at an overly narrow expected
count for positive A1/D5 credits; the final guard records all four.
No failure of a payment inequality or finite certificate occurred.
Fresh saved replay session10149 terminated exit0/PASS (terminal
chunk383bd1; saved loge69b1b), with exact saved-object equality and
all upper/attainment checks. After the descriptive path-terminology
precision, final source-bound generation4937/3a01df and fresh saved
replay81754/2b1408 both terminated exit0/PASS (logs2872db andc64781).

The arithmetic peer completed a full producer/report mathematical read,
including the inherited full-seal and actual/upper/context/item
interfaces. The review passed the three distinct whole-corner payments,
retention of broad halves and finite-D17 short-leaf units, complete
current coupon inventory and literal composition, exhaustive rank5
fixing, the ordinary/no-rank2 contextual split, simultaneous type
domination, all-last-item upper/attainment design and actual397.5
witnesses. It did not execute geometry, census, finite tables or scalar
arithmetic.

The ideal/finite peer also completed a full producer/report mathematical
read. It passed all four arbitrary-support seals and marked roles,
disjoint short-leaf/center obligations, current D14 and contextual
interfaces, exhaustive residual decomposition, type domination and
attainment. Its requested prose precision now explicitly identifies
3–white–13 as the short path at the adjacent black sides of actual
corner(3,13). It did not execute seals, maps, literal or Bellman checks.

Root completed a full current producer/report mathematical composition
read, including that terminology precision. The three distinct payments,
whole-corner sealing, original-half output inventory, actual and
representative increments, complete rank5 fixation and ordinary/no-rank2
context split passed. This review did not execute arithmetic or rerun
the complete inherited census. All three reviews are internal agent
reviews. Final metadata validation checks317 direct bindings,10 local
report links and source syntax. Source, certificate and report are final
and frozen.

The calculation reconstructs four exact native full-core seals and
all216 literal interfaces. Complete historical map enumeration,
all-label supporting classifications, the original exhaustive
low/disjoint/high-overlap reduction and the current triangle and
short-path theorems remain explicitly inherited frozen inputs.
No numerical normal search, geometric census, ideal scalar replay,
external mathematical review or formal proof-assistant verification
is claimed. The attained forest costs are formal counting statements,
not finite-volume or compact hyperbolic examples.
