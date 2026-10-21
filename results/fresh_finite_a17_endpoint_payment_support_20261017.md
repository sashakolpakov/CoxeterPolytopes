# The whole finite A17 endpoint unit is redundant

This package audits the complete support of one proposed change to the frozen
dimension69 face-payment rule: at a finite corner whose two marks are the
opposite endpoints of a **whole finite A17 component**, remove its unit.
Canonical A17 is the path0--1--...--16, so the removed pair is(0,16).
The same-component cost is136 before and135 after this change. The pair is
fixed by path reflection. All102 short8 and actual-leaf9 pairs are retained,
and every cross-component coupon is unchanged, including A17 positions5/11.

The complete local provenance audit below finds no lost triangle or
quadrilateral payment: triangles still receive at least2 and quadrilaterals
at least1. This is a local refinement conditional on the frozen complete
payment theorem and its subsequent transfers. It asserts no ambient
dimension improvement. In particular, it does not recompute the optimum
over finite forests, and the remaining ideal scalar obstructions have not
been discharged by this change.

The [checker](../scripts/fresh_finite_a17_endpoint_payment_support_20261017.py)
and [certificate](fresh_finite_a17_endpoint_payment_support_20261017.json)
bind the [complete finite8/ideal5 theorem](fresh_finite8_ideal5_combined_payments_20261017.md),
its complete all-label overlap and disjoint domains, and every later local
transfer used in the [dimension69 composition](fresh_ideal5_weighted_dimension69_20261017.md).
The distinction between actual whole components and displayed subdiagrams
is retained throughout.

## The first reductions do not need the removed pair

The triangle theorem uses finite marked distance at most8, or retained
ideal elementary pairs and whole affine factors of order at most9. A pair
at distance16 in a finite A17 cannot supply any of these payments.

For a quadrilateral with an adjacent marked path of length at most8, a
finite corner pays directly. For an ideal corner, the universal cycle and
tree replacements use finite short8, distance9 at an actual whole-component
leaf, or the existing coupons between separate components. The finite8
cycle theorem covers the long cycle arcs and radius6 tree branch; the
ideal5 local replacement covers arc6 and radius5 trees. Their finite
single-profile conclusions are stated without blanket small-component
units. The separate cycle13 and cycle14 proofs likewise use finite paths
of length at most7/8 or retained cross-component coupons.

The sealed cycle18 equality is an especially relevant interface: its two
half-payments are between whole A1 and whole A17 at positions5/11. They
do not have two marks in A17, and they do not depend on a same-component
endpoint unit. The checker records both literal corners and confirms their
halves remain unchanged. The cycle16 and cycle14 coupon components are
A15 and A13 respectively. No triangle, short-path or universal cycle
branch loses a payment.

We may therefore assume no marked path of length at most8. The low-core
theorem's66 surviving overlaps have a whole marked-component bound12.
The all-label high-core disjoint theorem excludes that entire branch.
These results are independent of the removed A17 pair.

## All2,614 high-overlap contexts

The frozen82,197-identification census and its no-extra-cross-edge proof
leave exactly2,614 supported contexts. The new checker routes **every**
one as follows.

|Route|Contexts|Why the endpoint unit is unnecessary|
|---|---:|---|
|Whole marked component at most16|2,391|Cannot be whole A17|
|Protected finite leaf9 or basic ideal feature payment|74|The distance9 unit is retained|
|Adaptive finite payment|9|Actual leaf9 or whole component at most16|
|Bound17, no induced finite-A17 embedding|125|The selected connected component cannot lie in A17|
|Bound17, induced finite-A17 maps|15|All48 maps avoid the endpoint pair|

For the140 bound17 cases, use the saved actual core Gram and its selected
connected marked component. In a positively gauged Gram, a native simple
bond is characterized by `4*g_ij^2 = g_ii*g_jj`; this avoids confusing a
rational congruence with norm-two normalization. An induced connected
subdiagram of A17 must be a simple path. A connected path embeds into
A17 precisely as a consecutive interval, in either orientation. The
checker enumerates all such maps and checks every edge and nonedge.
This over-approximates actual completions: excluding the endpoint image
from this larger set requires no additional extension assumptions.

Exactly15 selected components are paths admitting maps, with48 maps and
10 distinct canonical marked pairs:

`(2,11),(2,12),(2,13),(2,14),(3,12),(3,13),(3,14),(4,13),(4,14),(5,14)`.

None is(0,16). The other125 selected components fail the connected simple-path
test and therefore cannot lie in finite A17. This also has a direct
graph interpretation: any connected subdiagram of A17 containing both
endpoints already contains all17 nodes. A proper extension cannot turn
the marks of one of the listed path patterns into both endpoints.

The9 adaptive contexts were not folded into this embedding count. Their
independently proved alternatives are4 protected leaves,3 adaptive choices
of a protected leaf, and2 leaf-or-tail cases with whole bound16. Their
distance and whole-order interfaces are checked separately. The same is
true of the74 protected contexts. Thus the48 maps are not presented as
the entire quadrilateral domain.

## The old literal D16/D17 replacements

The complete small-factor map theorem matches18 literal old replacement
contexts,9 each for affine orders16 and17. Their proof says either an extra
supporting root gives a retained length-two marked path, or every extra
supporting white is orthogonal to the full core. In the latter case the
displayed finite components are whole actual components.

The new checker reconstructs every finite corner from the saved literal
Gram and black marks, verifies its components and positive definiteness,
and subtracts the proposed endpoint unit when appropriate. Seventeen
cases lose no payment. The only affected case is old order17 case662.
It originally had three finite units. Its marks(3,13) are indeed the
endpoints of a whole finite A17, so that unit is removed. The other two
finite marked components are both D17, each retaining its whole-component
unit. The new finite total is therefore2, still above the required1.

This check is necessary: absence from the preferred overlap maps alone
would not have covered this nonpreferred finite corner. The sealing and
length-two alternatives also explain why arbitrary additional supporting
whites cannot change this conclusion.

## Universal D13/D14 transfers and their immediate alternatives

The two universal opposite-fork transfers each have72 sealed equality
configurations. The checker reconstructs all144 complete Grams and all
corners, verifies the saved signatures and components, and recomputes the
retained unit-plus-half total. Every total is at least1.

The D14 equality cores have18 roots and16-root corners, so no A17 component
can occur. The D13 equality cores have19 roots. Exactly8 of their finite
corners have whole A17; the marked positions are(0,13) in the deterministic
path orientation, equivalent to(3,16) after reflection. None has marks
at both endpoints. All8 units survive. The finite nonextreme zero cases
have whole D12 components and also survive.

The universal proofs' immediate ideal alternatives have orders12 and13,
not finite A17. Their refinements in the dimension71 composition are the
universal affine-D12 transfer, the B12 singleton strict gap and the cycle13
transfer. These are included below. Consequently this argument covers
the immediate and universal uses as well as selected-core occurrences.

## Later local transfers through dimension69

The new checker visits the explicitly typed finite-corner interfaces of
all subsequent transfers, rather than inferring safety from their names.

|Transfer|Finite units/coupons that can supply its conclusion|
|---|---|
|Order17 D singleton|Separate A4/D13 halves; length-two alternative|
|Order17 D partial weights|Whole D17 units or separate A4/D13 and A5/D12 halves|
|Order16 D partial weights|Whole D16 units or separate A4/D12 and A5/D11 halves|
|Order15 B three pairs|Whole D14 or whole A12/A13 finite alternatives|
|Order13 B singleton and promoted long pair|Whole B13 or separate A1/D12, A2/B10 coupons|
|Universal affine D12 opposite leaves|Whole D13/D14/A15, full A1/D12, or separate A2/D12 halves|
|Order14/order15 B singleton zeros|Strict projection exclusions; no finite alternative needed|

The D16 partial-weight proof's proper-completion exclusions remain valid
before any finite payment is used. In its10 sealed cases the finite corners
have16 roots. The weighted D17 source supplies12 sealed19-root cores;
its only same-component units have typeD17. The A4/D13 new singleton
coupon remains between different components. The order13 B long-pair
promotion is read in its actual dimension71 composition: two contexts
give whole B13 and two give the full whole-A1/D12-short-leaf coupon.

For the universal D12 transfer, arbitrary one-sided spectators leave a
whole A1/D12 or D13 corner; two-sided energy-one cases give a short path
or two sealed A1 spectators. Its literal same-component units are D13,
D14 or A15. The whole-component arguments and arbitrary-support scope
are inherited from that universal theorem, not extrapolated from the
finite table. None of these routes needs the finite A17 endpoint unit.

## Status and verification scope

This establishes the local refinement by exhaustion of the predecessor's
coverage routes and all later transfer interfaces. It preserves the
triangle minimum2, quadrilateral minimum1 and every ideal pair weight.
Only the one whole finite A17 endpoint pair is changed. It does not claim
that any other unsupported long pair can be deleted.

The checker independently enumerates the140 A17-eligible overlap
interfaces and48 path maps, reconstructs18 literal old replacements and
144 universal D equality configurations, checks all relevant successor
finite-unit types and cycle18 halves, and verifies all direct frozen
bindings. It inherits the full all-label census, the arbitrary spectator
completion/sealing arguments, and the published finite/affine diagram
classifications. It does not rerun those censuses or assert external
refereeing or proof-assistant formalization.

Final generation27989 and separate fresh replay86207 both terminated
exit0/PASS. All78 direct bindings, syntax and4 local report links pass.
Root completed a full source and mathematical report read. That review
passed complete2,614-context routing; the native-simple path test in
positive gauges and all48 induced maps; the sole nonpreferred legacy662
debit retaining two D17 units; all144 sealed D equality configurations
including eight nonendpoint A17 units; cycle18 separate-component halves;
and the successor typed-unit interfaces. The full geometric census and
arbitrary-support sealing proofs remain inherited inputs and were not
rerun by that review. Source, certificate and report are final/frozen.
No global dimension or exact finite-forest maximum is asserted here.

    python3 scripts/fresh_finite_a17_endpoint_payment_support_20261017.py
    python3 scripts/fresh_finite_a17_endpoint_payment_support_20261017.py --verify
