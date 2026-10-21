# Sparse finite coupons lower A16/D16/D17 costs to114/114/125

Relative to the completed dimension63 payment rule, finite same-component
costs can be changed to **A16=114, B16=104, D16=114 and D17=125**, with
the six whole-component coupons below. Every triangle still receives at
least2 and every quadrilateral at least1. Ideal weights are unchanged.
A separate complete finite-forest count is required; this local theorem
does not assert an ambient dimension bound.

The [source](../scripts/fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018.py)
and [certificate](fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018.json)
bind the [completed63 rule](fresh_order11_order17_dimension63_20261018.md),
the [preceding finite exchange](fresh_finite_d16_d17_two_unit_exchange_20261018.md),
and its complete selected and literal interfaces. Every frozen predecessor
is preserved.

## The new intrinsic finite rule

Use path numbering0 through n−1 for A_n. For D_n use the spine0 through
n−3, with its two short leaves n−2,n−1 at n−3. Delete these unit pairs:

|Whole component|Deleted pairs|Old cost|New cost|
|---|---|---:|---:|
|A16|(0,12),(3,15)|116|114|
|D16|(0,10),(0,13)|116|114|
|D17|(1,12),(1,13)|127|125|

These sets are invariant under the complete labelled automorphism groups.
Every marked distance is at least10. All finite short8 pairs and protected
whole-leaf9 pairs remain units. Existing halves, other deletions, deep11
distance8 zeros, and all old coupons retain their values.

Add the following payments between distinct **whole** finite components.
An endpoint or second node is counted at either end of an A path; the
D5 long endpoint is the leaf on its length2 arm. The two E6 long endpoints
are the leaves on its two length2 arms.

|First marked component|Second marked component|Weight|Total per component pair|
|---|---|---:|---:|
|A4 endpoint|A12 endpoint|1|4|
|A4 endpoint|A13 second node|1|4|
|D5 long endpoint|A12 second node|1|2|
|either A2 node|A13 endpoint|1/2|2|
|either A2 node|A14 second node|1/2|2|
|E6 long endpoint|A11 second node|1|4|

These roles are disjoint from the old coupons on the same type pairs.
For example the old A2/A13 cycle coupon marks the middle of A13, whereas
the new half marks an endpoint. No cost is assigned merely because a
larger component contains one of these diagrams as a subdiagram.

## Complete selected support and proper D16 completions

The complete finite induced-embedding routine is rerun on every eligible
selected component in the frozen2,614-context high-overlap domain.
The two A16 deletions occur only in contexts1651,2026. The D16 deletions
occur in501,504,1031,1282,1534,1537,1614,1621,1688,2547. No selected D17
image uses either removed restored pair; those units were introduced
solely for the previous exchange's four contexts1681,1685,2500,2532,
which are handled explicitly below. Thus map absence is not used to
discard their actual predecessor obligation.

Six of the D16 contexts are proper completions with a nonsingular
Lorentzian core. All roots are normalized to square2. The inherited
contact theorem partitions possible contacts into two ends. A root
touching both ends gives an actual short8 path; within each end the
inverse block is entrywise positive. Hence unspecified additional
contacts can only increase a projection estimate. A first target root
adjacent to the old marked component is genuinely new: an old other
corner component is orthogonal to that component, and an omitted black
root does not occur at this corner.

For501 and1534, a missing short fork leaf must have a unit contact at
an old node with inverse diagonal6. Its projection has square at least6,
contradicting its norm2.

For504 and1537, two distinct missing short fork leaves are orthogonal
in the required D16 diagram. Both contact the same old node with inverse
diagonal3/2. Each residual square is at most1/2, but their residual inner
product has magnitude at least3/2. Indeed both contact vectors lie in the
same positive inverse end block, while their actual inner product is0.
The residual two-by-two determinant is at most−2. The orthogonal complement
of the Lorentz core is positive definite, so this is impossible.

For1031 and1614, the first missing root is the fork center. Its forced
unit contact has projection square4/3. A nonsimple contact would already
give at least8/3. Any additional contact in the same end raises the
projection beyond2; the other end gives a short8 path. Thus this center
has exactly the one core contact and residual square2/3.

Its two prospective short leaves cannot be old other-corner roots,
since the center is orthogonal to all those roots. A new leaf has no
contact with the old marked component, as prescribed by the actual D16
diagram. For every remaining permitted core position, the exact test
gives one of three obstructions: projection greater than2; projection2
and a strictly positive product with the forced center; or a short path
of length at most8 through the required center–leaf bond. Same-end
multiple contacts also have projection greater than2. Therefore both
new leaves are orthogonal to the entire core. Their actual mutual product
is0 and each has product−1 with the center. The residual principal is

    [ 2/3  -1  -1 ]
    [  -1   2   0 ]
    [  -1   0   2 ],

with signature(2,1,0), contradicting positive definiteness of the residual
space. These arguments allow arbitrary extra supporting roots and all
native contact magnitudes; unspecified products are never set to zero
without the stated necessity.

## Ten sealed radical alternatives

All remaining selected cases, and all four old exchange cases, have
the same19-root L25/L25 core identified at node9. Its signature is(17,1,1),
with radical

    (−6,−3,−4,−2,−5,−4,−3,−2,−1,0,
      6, 3, 4, 2, 5, 4, 3, 2, 1).

The radical represents the zero vector in an actual Lorentz realization:
it is null and orthogonal to a span containing a timelike vector. Every
additional supporting white must annihilate it. All permitted contact
positions have nonzero radical coefficients; every opposite-sign pair
gives a checked short8 path. Therefore, unless such a path pays, every
additional white is orthogonal to the whole core. Its displayed finite
corner components are whole and cannot acquire an indirect tail or bridge.

The source reconstructs all four finite corners, their full components,
and the marked roles, obtaining these alternatives:

|Contexts|New finite payment replacing the deleted unit|
|---|---|
|1651,2026|A4-endpoint/A12-endpoint unit|
|1282,1621|A2/A13-endpoint half plus distinct A2/A14-second-node half|
|1688,2547|E6-long-endpoint/A11-second-node unit|
|1681,2500|A4-endpoint/A13-second-node unit|
|1685,2532|D5-long-endpoint/A12-second-node unit|

These ten alternatives use no deleted same-component unit. The latter
four discharge exactly the obligations previously paid by restoring
D17(1,12),(1,13); the restoration is not simply undone without replacement.

## Literal and universal compatibility

All216 literal configurations of the preceding finite theorem are
reconstructed and debited simultaneously. Two require the new A2/A13
endpoint half: cases394 and648 of the order16 finite-cross package lose
their D16(0,13) unit but retain an old A2/D13 short-leaf half at another
corner. The new half restores total1. The old affine-kernel/inverse
argument in that package seals every finite component at these corners,
so the new coupon is applied to whole factors. Two other configurations
gain the A4/A12 unit. The remaining212 totals are unchanged, and the
minimum over all216 remains1. Literal compatibility is not a claim that
these216 configurations enumerate every quadrilateral.

The complete low-core bound is12; the all-label high-disjoint domain is
excluded. The protected/adaptive branches use basic ideal or whole-leaf9
payments. Their two bound16 trees have distance9 and miss every deletion.
Universal cycle/tree alternatives retain finite short paths, protected
whole-leaf9 pairs, small components, or separate coupons. The later A16
long-pair uses are among the literal configurations just checked. The
later baseline63 order11 outputs are whole finite components of order
at most14, A9/A18 endpoint8, D18 leaf9, E8 path6, or typed coupons. The
order17 E6 transfer uses no finite same-component D16/D17 unit.

If a required short8 payment is a deep finite8 zero, the independent
depth11 corridor supplies its constructive short7 alternative. The
completed63 small-tree, cycle-tail and protected-feature interfaces pay
that alternative without using the present distance10-or-longer deletions.
Every finite short7 unit and every ideal weight is unchanged. Thus this
argument does not invoke the unchanged-corner branch of the new theorem
to justify a payment used by its changed-corner branch.

Triangles retain the completed finite7/ideal5 coverage. The new finite
coupons only add payments. The selected-component classification,
proper-completion exclusions, radical alternatives, complete literal
audit and universal interfaces therefore prove simultaneous face payment
under the baseline63 ideal rule. Future ideal refinements must provide
their own compatible changed-ideal-first composition; none is assumed here.

## Verification and counting scope

Generation42537 and separate fresh saved replay28025 both terminated
exit0/PASS. The complete selected-map domains contain394 eligible A16
patterns/238 maps,394 D16 patterns/578 maps, and140 D17 patterns/168 maps.
The new checks reconstruct eight proper-map exclusions in six contexts,
ten radical transfers, complete canonical roles and all216 literal totals.
All143 direct bindings, Python syntax and four local links pass.

Root and a separate reviewer completed full producer/report mathematical
reads. Root also read the inherited complete selected/literal routines and
the whole-component sealing proof for old cases394/648; the separate
reviewer read the inherited endpoint-constraint helper. Both checked the
first-entry argument, all extra/native contact alternatives, two-leaf
Cauchy and fork residual exclusions, radical whole-factor control, exact
new coupon roles, restored-unit obligations, simultaneous literal debits
and baseline63 short8/deep11 interfaces. Neither reviewer re-enumerated
selected maps, inverses, historical classifications or literal arithmetic.
Historical all-label classifications remain inherited theorem inputs.
Source, data and report are final/frozen. All reviews are internal; no
external refereeing or proof-assistant verification is claimed.

The new actual coupons must be counted, rather than declared free. Some
can be dominated by existing representative-type budgets; the E6/A11 unit
may require a larger conservative E6/D11 representative allowance. Such
an allowance is an upper bound, not an extra actual payment rule. A
separate complete forest certificate must include all interactions before
these local costs are used in a dimension comparison. No global optimizer,
new ambient bound, or example is supplied by this packet.
