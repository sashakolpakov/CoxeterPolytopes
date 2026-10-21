# Measured payment targets from the ambient58 obstruction witnesses

The most promising coupon target in this saved-witness diagnostic is a
common decrease of1/8 **per small marked node** in the ordinary payments
between whole components of ranks3–5 and whole D15. It clears all33
recorded witnesses it affects, including every saved maximizer in seven
of the23 obstruction classes. The exact required decrease for those33
witnesses is19/160<1/8. Sixteen classes still have a recorded counterexample.
This is a target for a new geometric proof, not a proved payment reduction.

The pure-base obstruction A16+3D14 is independent of those coupons.
If d16 and d14 denote proposed decreases in the corresponding component
base costs, it imposes

    d16+3*d14 >=19/4.

No single one of the49 individual or grouped cost features considered
here clears all recorded witnesses. The current completed unrestricted
bound remains58, and no dimension57 statement follows.

The [source](../scripts/fresh_ambient58_obstruction_payment_sensitivities_20261018.py)
and [data](fresh_ambient58_obstruction_payment_sensitivities_20261018.json) read only the frozen
[D12 type-split count](fresh_d12_type_split_ambient58_forest_frontier_20261018.md).
That certificate lists26 above-curvature representative classes; exact
actual-type refinement leaves23. The present source extracts all114
stored actual attaining forests with charge above377, including64
class-maximizing witnesses. It recomputes their full current charge,
rank and coefficient multiplicities. It does not run a new dynamic
program, optimizer, geometric scout or cusp convolution.

## What the sensitivities mean

A component-base decrease d contributes d times the number of those
components. A decrease q in an ordinary **total component-pair cost**
contributes q times the number of those unordered component pairs.
A per-small-mark decrease delta against D15 contributes

    delta * (3*c3+4*c4+5*c5) * cD15,

where cr counts all whole components of rank r, across their actual
finite types. The source checks that their present D15 component-pair
cost is3r/2 for every actual type of rank3,4 or5. That aggregate identity
does not prove which individual geometric coupon role can be decreased.

For example, lowering one intrinsic D15 fork-center role from1/2 to3/8
for every small rank3–5 mark would have the measured delta=1/8 effect,
**if** a complete proof preserves every use of that role. Reducing the
two short-leaf roles by1/16 each would have the same aggregate cost
but is a different geometric proposal. Neither proposal is asserted
valid here; the center and leaf obligations must be audited separately.

Every recorded witness currently has zero contextual charge. This fact
is checked exactly, but does not mean the contextual rules can be
removed from a future complete count: other forests can trigger them.
Any new compensation must be included before claiming an improved
maximum.

The individual features and grouped features are separate, sometimes
overlapping parameterizations. Their reported savings are not a menu
of independent coupons that may be added without checking overlap.
A component-pair cost unit is not automatically one marked-pair weight;
component rank, diagram symmetry and orbit size determine that conversion.

## Concrete necessary constraints

The following inequalities are necessary to remove their displayed
witnesses, assuming only the decreases named in each row. They are not
sufficient for the complete forest problem.

|Witness|Current excess over377|Necessary savings|
|---|---:|---|
|A16+3D14|19/4|dA16+3dD14≥19/4|
|A3+2A5+3D15|9/2|3dD15+3qA3,D15+6qA5,D15≥9/2|
|2A15+2D14|4|2dA15+2dD14≥4|
|A1+9A2+3A13|61/16|3dA13+3qA1,A13+27qA2,A13≥61/16|
|2A4+A5+3D15|7/2|3dD15+6qA4,D15+3qA5,D15≥7/2|
|2A16+D26|7/2|2dA16+dD26≥7/2|
|A16+D17+D25|25/8|dA16+dD17+dD25≥25/8|
|7A1+3D5+3D12|23/8|3dD12+21qA1,D5+21qA1,D12+9qD5,D12≥23/8|

The D15 coupon candidate clears the leading coupon-heavy witness once
delta≥3/26, because its coefficient is39. Across all33 affected saved
witnesses, the sharper requirement is19/160, attained by
5A1+2D5+2D14+D15: its excess is19/16 and coefficient10. Thus delta=1/8
clears all those recorded witnesses.

A literal A5-only coupon reduction would miss equal witnesses using B5
or D5. Likewise an A15 base decrease alone can leave tied D15 forests.
The grouped parameters explicitly expose these type ties.

The next independent coupon family is A2 against A13 or D13. Decreasing
both total component-pair costs by61/432 clears all ten affected saved
witnesses, including every maximizer in two classes. Reducing only A2/A13
leaves alternatives using D13. Its geometry and any increased witness
payments remain unproved.

## Ranking the bounded probes

Each row below probes the specified cost decrease and records how many
of the23 classes still contain one of the114 known counterexamples.
It does **not** report the number of classes whose reoptimized maximum
would still exceed377. Unrecorded forests can replace the old maximizers.

|Cost feature|Probe decrease|Classes with a surviving saved witness|
|---|---:|---:|
|proportional_coupon:rank3to5 / D15|1/8|16|
|base:D14|1|18|
|coupon_group:rank5 / D15|1|17|
|base_group:A15_D15|1|20|
|coupon_group:A2 / A13_D13|1|21|
|base:A16|1|23|
|coupon:A5 / D15|1|23|

The data retain the full49-feature ranking, exact rational thresholds,
coefficient vectors and every saved witness. The ranking measures only
these witnesses; it does not estimate geometric proof difficulty.

## Verification and scope

Final generation9614a5 and separate fresh saved replay36231b both
terminated exit0/PASS on23 classes,114 witnesses and49 features. The source independently
reconstructs all multiplicities and current full charges, verifies every
class-maximizer count against the frozen certificate, and recomputes
all rational probe/threshold statements. It makes no new global-optimum,
local-theorem, universal-coupon or geometric-realizability assertion.
A complete changed-rule count is still required after any candidate
geometry and all compensating costs are proved.

The parent completed a full saved-accounting source/report read and
passed the exact rank/cost/context checks, actual-type ties, per-component
versus per-mark units, multiplicities, grouped overlap,114/64/23 witness
scope and necessary-threshold limitations. That review executed no DP,
geometry or other arithmetic. Source, data and report are final/frozen.

    python3 scripts/fresh_ambient58_obstruction_payment_sensitivities_20261018.py
    python3 scripts/fresh_ambient58_obstruction_payment_sensitivities_20261018.py --verify
