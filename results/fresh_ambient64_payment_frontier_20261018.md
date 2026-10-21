# Exact ambient64 payment frontier

Under the completed [dimension64 payment rule](fresh_deep_finite8_dimension64_20261018.md),
the exact finite-vertex maximum at ambient64 is **472**, against curvature464.
The complete current ideal comparison has **1,513 terminal profiles**, of
which24 have negative margins. Its minimum is−62665/419628 at affine ranks
(10,10,11,16,16). Thus the current rule does not exclude dimension64.

The [source](../scripts/fresh_ambient64_payment_frontier_20261018.py) and
[certificate](fresh_ambient64_payment_frontier_20261018.json) also give
conditional targets. They prove scalar and finite-forest statements only.
No proposed pair deletion, replacement coupon, changed local payment
theorem or dimension63 bound is asserted here.

## Exact current finite maximum

All components of rank at most17 retain their current weights. The
[frozen finite budget](fresh_d16_d17_d14_coupon_forest_budget_20261017.md)
has exact rank64 maximum472 on that restricted domain, attained by four
whole D16 components, each of cost118.

For a forest with a component of rank r>=18, the proof of the
[depth17 budget](fresh_depth17_rank65_conditional_budget_20261018.md)
applies unchanged: its type is A/B/D, it participates in no current
cross-component coupon, and D dominates its modified same-component cost.
Remove that component and bound the remaining modified forest by the old
unrestricted maximum at rank64−r. The47 exact upper bounds have maximum469,
attained at r=64 by whole D64. The next bounds are467 at r=48 and933/2 at
r=47. Hence this branch is strictly below the all-small maximum472.

This proves the sharp current maximum without assuming that the remainder
has small components. A separate conditional optimizer below also returns
472 in its baseline calibration.

In particular, if coupons stay unchanged, the witness4D16 requires reducing
the D16 same-component cost by at least2 to reach464. That change alone is
insufficient: whole D64 still costs469. These are formal finite vertex
forests; the calculation does not assert that any such vertex occurs in
an actual dimension64 polytope.

## Complete current ideal domain

For affine rank r the current normalized factor envelope is r/2 for
r<=10; then16/3,74/13,6,94/15,107/16,7 for r=11,...,16; and163/23 for
r>=17. These are the unchanged envelopes of the
[ambient65 ideal closure](fresh_ambient65_ideal_margin_closure_20261017.md).

The inherited [convex-order theorem](fresh_distance14_unrestricted_dimension156_20261017.md)
allows concentration of two low ranks until one vanishes or the other
reaches10, and of two ranks at least18 until one reaches17. Each move
decreases normalized curvature and preserves cost, because the relevant
cost ranges are linear and constant. Thus a complete terminal domain at
link rank63 has arbitrary multiplicities of10,...,17, at most one rank1,...,9,
and at most one rank>=18. The source enumerates all1,513 such partitions.

For each partition r, it constructs

    H(t)=product_i(1+t+...+t^r_i), V=H(1),
    G=sum_j [t^j]H(t)*min(j,63−j),
    margin=(5G−63V)/(2V)−sum_i c(r_i).

Every polynomial, G, V, exact cost, normalized margin and absolute margin
is saved. There are24 negative margins and no zeros. The main deficits are:

|Affine ranks|Normalized margin|
|---|---:|
|(10,10,11,16,16)|−62665/419628|
|(10,10,10,10,10,13)|−29765/204974|
|(1,10,10,10,16,16)|−50224/384659|
|(10,10,10,10,11,12)|−12997/103818|
|(10,10,12,15,16)|−48131/427856|
|(1,10,10,10,10,10,12)|−15835/190333|
|(10,10,10,11,11,11)|−5365/69696|
|(9,11,11,16,16)|−5/2601|
|(3,10,10,10,10,10,10)|−1927/644204|

All24 negative rows are retained in the certificate, not just this table.
The last profile has absolute deficit21197 and contains no affine rank
above10. Therefore improvements confined to factor orders>=12 cannot
settle the ideal comparison. A reduction at order11, a lower-order change,
an exclusion of that actual cusp shape, or a new global compensation
argument is necessary for this profile. The formal product profile is not
claimed geometrically realizable.

## Conditional ideal targets and the necessary expanded domain

Let a be a reduction in the maximum **total pair cost at factor order11**
(affine rank10), and b a reduction at factor order17 (affine rank16).
Their contributions to a profile's normalized improvement are

    a*count(rank10)/11 + b*count(rank16)/17.

Reducing the rank10 cost destroys the old low-rank linear range through10.
It would be invalid to test only the original1,513 profiles. Instead retain
arbitrary multiplicities of9,...,17, at most one1,...,8 and at most one>=18.
The low range through9 and the high constant range remain unchanged.
The resulting **2,419 profiles** are a complete reduction domain for every
scenario here. They include the1,513 original rows with identical scalar
data; the current baseline has26 negatives on this larger domain.

|Conditional total pair-cost changes|Negative profiles|Minimum margin|
|---|---:|---:|
|Order11:55→54|1|−5/2601|
|Order11:55→54.5; order17:119→118.5|0|167/419628|
|Order11:55→54.5; order17:119→118|0|2385/146432|
|Order11:55→54; order12:64→63|0|177/4624|
|Order11:55→54.5; order12:64→63|2|−9235/427856|

With a=1/2, the exact necessary and sufficient b threshold for nonnegative
margins on this domain is24517/49368, attained at(10,10,11,16,16).
The choice b=1/2 exceeds this threshold and gives the displayed strictly
positive minimum. With b=1, the corresponding a threshold is12997/37752,
attained at(10,10,10,10,11,12). If instead order12 improves by1, the a
threshold is48131/77792 at(10,10,12,15,16). No profile independent of the
varied parameter remains negative in any of these three threshold checks.

These are **all-family envelope targets**, not changes to just one affine
type. In particular affine A10 is an11-cycle: all55 pairs are already at
distance at most5. Lowering the order11 maximum requires a new transfer
for an elementary short5 pair in that cycle, in addition to any needed
B/C/D changes. Deleting unused long pairs only in tree families cannot
achieve this target. The hypothetical order17 reduction must lower its
actual worst family, presently affine D16 of total cost119; the other
families have smaller totals.

Any additional finite coupons needed for these ideal transfers must be
counted separately. The scalar improvements here do not grant a free
finite replacement.

## Conditional finite targets with all existing coupons retained

For clarity, the diagnostic separately computes the complete rank64 finite
maximum under the following proposed same-component changes. Nothing in
this section proves their local geometric availability.

|Deep distance8 rule|D16 cost|D17 cost|Exact rank64 maximum|Attaining forest|
|---|---:|---:|---:|---|
|Current depth17 half|118|125|472|4D16|
|Depth17 half|116|125|469|D64|
|Depth17 zero|116|125|465|D17+D47|
|Depth17 zero|116|124|464|D17+D47|
|Depth11 zero|116|127|464|4D16|

The last row is the parent's proposed two-unit exchange: remove two D16
units and restore two D17 units, combined with a prospective depth11 zero
rule. Its finite arithmetic is sufficient even though its D17 cost rises.
Its local proofs are not inputs to this diagnostic and remain separate.
For every row, A16=116, B16=104, A17=114 and B17=105 stay unchanged.

The all-rank counts of affected distance8 pairs in finite A/B/D are
max(0,r−2t−8), max(0,r−2t−9), max(0,r−2t−10), respectively, for feature
depth t. They follow from the same path-start intervals as in the frozen
depth17 budget. A zero weight subtracts two per pair from doubled old cost;
a half subtracts one. No rank<=17 component is affected at either depth.

For each scenario the source independently regenerates all753 multisets
of components D10,...,D16,A17 with total rank<=64, and checks equality with
the frozen domain restricted to64. Their unchanged coupon interactions
give348 five-coordinate environments. Residual types are A1,...,A9, with
E6 replacing A6 when it dominates, and D17,...,D64. A17 remains a separate
large type because of its distinct coupon. Finite type domination and all
38,025 type/partner coupon inequalities are checked in each scenario;
the uniform I2 representative covers all dihedral labels.

An item-first integer Bellman calculation is followed by all597,168
last-item upper inequalities and22,272 attaining-predecessor checks per
scenario. The literal component/coupon formula is checked42,921 times;
all753 global alternatives have a reconstructed typed attaining forest.
These checks account for arbitrary numbers of components and every
unchanged coupon, including E6/D11 and A1/A2–D14 roles. They are not a
sample of promising forests. The saved tables and the fresh replay are
part of the certificate.

If a future compatible local theorem gives finite maximum<=464 and the
positive ideal target above, the inherited global additive constant at
ambient64 is5120. Zero finite slack would be sufficient: the ideal surplus
and positive global constant would give a contradiction. This diagnostic
does not establish those local hypotheses or a dimension63 theorem.

## Verification scope

Final generation16c58d and separate fresh saved replay244610 both
terminated exit0/PASS. A separate full source/report mathematical read
passed the current strip bound, corrected2,419-profile concentration
domain, all-family target distinction, complete753-vector/348-environment
finite decomposition, type domination, explicit coupon formula and
Bellman upper/attainment proof. The reviewer compared the coupon formula
with its frozen provider definitions; it did not independently rerun
the scalar or dynamic-programming arithmetic. No gap was identified in
this diagnostic scope. The current finite proof binds the frozen sharp
all-small table and the proved deep17 large-type/coupon interfaces. The
ideal and conditional finite arithmetic here imports no project arithmetic
module; it uses exact Fraction/integer operations and explicit formulas.
The inherited convex-order theorem, finite type classification, current
geometric payment theorem and coupon interpretations remain bound proof
inputs. No historical geometric census is rerun. Direct source bindings,
Python syntax and local report links pass. Source, data and report are
final/frozen. Reviews are internal; external refereeing and proof-assistant
formalization are absent.

    python3 scripts/fresh_ambient64_payment_frontier_20261018.py
    python3 scripts/fresh_ambient64_payment_frontier_20261018.py --verify
