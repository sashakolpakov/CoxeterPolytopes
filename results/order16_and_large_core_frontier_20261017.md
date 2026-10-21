# Order16 replacements, dimension86, and larger arithmetic cores

The unrestricted dimension ceiling is now **86**, improving the completed
91 argument. This is a complete internally checked computer-assisted
proof for all finite-sided, finite-volume real hyperbolic Coxeter polytopes.
It assumes neither arithmeticity nor a restriction on cusp factors.
Primary generation, fresh replay, independent exact reconstruction of
all3,355 final profiles and separate internal mathematical reviews pass.
External refereeing and proof-assistant formalization remain unestablished.

Compact arithmetic remains at most21, arbitrary compact at most29,
and the compact-facet ambient ceiling remains79. The earlier29 bound
under the minimum affine-factor-rank9 hypothesis is unchanged. No new
finite-volume construction has been found.

## How the91 obstruction was removed

The preceding rule paid every pair in an ideal factor of order16,
costing120. The new proof replaces certain ideal payments with finite
half-payments while preserving the total finite budget9n.

The [complete embedding census](fresh_order16_preferred_embeddings_20261017.md)
retains all260 preferred-component contexts whose extension bound can
reach16. Its complete induced labelled maps yield pair costs112,107,103,120
for affine A,B,C,D on16 nodes.

The [universal cycle replacement](fresh_cycle16_half_coupon_transfer_20261017.md)
proves that an antipodal pair in the16-cycle either has another short
paid corner or lies in the unique determinant-equality configuration.
In that configuration two finite A1/A15 corners have their A15 marks
at the unique middle node. Arbitrary supporting spectators either give
a length-two payment or are orthogonal to the entire core, so those
whole finite components persist. Half a unit at each finite corner
replaces the missing ideal payment. The cycle needs only112 pairs.

For D, the [combined local theorem](fresh_order16_finite_cross_payments_20261017.md)
removes nine pairs and leaves111. They occur only in nine literal
shared-edge L25 stars. The affine kernel forces additional supporting
whites to be orthogonal to the entire affine factor. A lone contact
with either remaining black root would have projected square at least4,
contradicting its norm2; contacting both gives a paid length-two path.
Thus otherwise the entire core is isolated from extra support. Five
cases already have a finite unit payment; four receive two A1/D15 or
A2/D13 half-payments.

The whole finite-component accounting includes these new payments and
the preceding A1/D16 and A2/D14 payments. A rank/slack inequality proves
total cost at most9n for every n<=108. The ideal rank envelope becomes

    b(r)=r/2 for r<=14;  b(15)=7;  b(16)=127/17;
    b(r)=217/29 for r>=17.

Its large-rank bound follows from all-rank affine path/fork formulas,
with the maximum at order29. The replacement arguments apply to actual
face stars in arbitrary ambient polytopes, not merely to saved small
examples.

## The dimension comparison

The [dimension86 proof](fresh_order16_piecewise_dimension86_20261017.md)
uses a new piecewise concentration argument. Reduced shapes have any
number of ranks14,16,17 and at most one each of a smaller rank, rank15,
and a rank above17. In particular, replacing two15 parts by14 and16
increases the budget while lowering curvature, which is the useful
direction for the contradiction.

All3,355 reduced profiles for ambient dimensions87 through91 have
positive curvature minus cost. The minimum is exactly146417/911250
at ranks(2,14,14,14,14,14,14), ambient87. The completed91 theorem supplies
the tail. The [independent implementation](fresh_order16_piecewise_dimension86_20261017_review.md)
uses different enumeration, sliding-window convolution and a cumulative
identity without importing project arithmetic; every value agrees.

At ambient86, the formal profile(1,14,14,14,14,14,14) has margin
-112409/1518750. The present comparison cannot lower the ceiling85.
This formal profile is not an example. Further local progress must
address order15 factors or reduce the distance-seven payments themselves;
global face-deficit allocation remains another possible route.

## Larger arithmetic cores remove the saved witness

The [new raw178/296/297 constraints](fresh_raw178_296_297_actual_corners_20261017.md)
cover all47 positive explicit/unique classes of those families in the
latest arithmetic witness:376 actual corners,7,411 literal orientations
and236 finite vertex forests. Core orders9/10 and spectator ranks6/5
are handled explicitly. Full contact vectors and future-height signs
are used throughout; all contact roots lie in Q(sqrt(5)).

Of761,872 literal pairs,760,379 are excluded by exact geometric,
arithmetic or unique-face tests. Unknown pairs remain allowed. The
saved incompatibility cliques and conservative compatibility colorings
give469 necessary capacity columns,91 of which violate the old witness.
Keeping all4,647 face masses fixed and raising65 vertex allowances
reduces its exact ratio from12.4575435307 to

    2877693523508930/247081190427751 = 11.646752707185115 < 12.

Within-family constraints alone already reduce it below12; mixed-family
constraints improve it further. This removes one fixed witness. It is
not an optimized feasible weighting rule, complete implicit-family
coverage or a new arithmetic dimension theorem. No new optimizer was
run. Any next enlarged solve must preserve all inherited obligations
and these469 columns, and any changed face masses need fresh implicit
coverage.

The signed columns in the inherited model are nonnegative variables
with some negative face coefficients. They are not free variables.
Their dual conditions are inequalities, so increasing a vertex allowance
preserves them. The new certificate explicitly binds the constructor
and solver/verifier documenting that convention.

## New construction geometry, still no closed chamber

The [A53 deformation](fresh_dotted_a53_disk_20261017.md) genuinely changes
the previously frozen prefixes by freeing A's old53 product. Nineteen
sampled values leave A native-admissible, but only the old control value
allows C under the saved AC4/5 recipes. This is a bounded failure,
not an all-parameter obstruction.

The [old10-anchored successor](fresh_c10_endpoint_disk_20261017.md)
instead prescribes C's old10 contact, then solves for A. It finds27
new native-admissible A/C pairs with A53 dotted and254 coherent five-cap
prefixes outside the old20-prefix exact obstruction. Every saved prefix
passes complete numerical cone and separate positive-definite support
checks, but all have24 finite rays,8 exterior rays and12 open ends.
All5,730 tested next-Z branches fail native contact gates.

These are new numerical domains, not exact constructions or all-label
exclusions. The certified control remains63 walls,199 finite vertices
and10 open target-face ends, with infinite volume. Another construction
step should change the mandatory XYRZ incidence or its defining finite
contacts, rather than repeat these blocked prefixes.

## Verification and archive

The [archive auditor](../scripts/audit_order16_and_large_core_frontier_20261017.py)
and [saved result](order16_and_large_core_frontier_20261017_audit.json)
extend the [preceding checkpoint](ideal7_cross_and_raw24_model_frontier_20261017.md)
without changing its frozen files. Eight new packages are included.

Independent small checks validate the260-domain filter and736 emitted
maps, reconstruct16 cycle rectangle values and the D-core inverse
contact block, and check all27 finite D corners, their whole components
and payments. They compare the two3,355-profile certificates and
validate all268 capacity cliques,201 colorings, actual saved occurrence
multiplicities and the exact arithmetic-witness repair. Construction
checks remain consistency checks of saved numerical records.

The archive does not claim to duplicate the complete embedding search,
native pair-projection census, physical construction arithmetic or
inherited geometric classification. Their generation/replay checks and
mathematical/source review scopes are stated in their reports. The
arithmetic source also received a later separate narrow review of its
variable core orders, whole-component maps, signed-column convention
and capacity construction; that review did not rerun the pair census.
All reviews described here are internal.
