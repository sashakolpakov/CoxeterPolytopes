# Twelve half-weights on the order17 affine D factor

**Completed local argument; verification is recorded below.** Starting from the frozen complete finite8/ideal5 face-payment theorem and its [D17 singleton successor](fresh_d17_singleton_pair_coupon_20261017.md), reduce twelve further ideal pair weights from1 to1/2. They form four automorphism orbits in affine D16, whose diagram has17 nodes. Its total pair weight falls from126 to120. No new finite coupon is needed.

This is a weighted rule. The twelve pairs are still paid by1/2; deleting them would not be justified by the present argument. Every affected selected ideal half has a different, actual finite corner carrying an existing half. Their sum remains1. All elementary triangle payments and all marked paths of length at most8 retain weight1. The [source](../scripts/fresh_d17_partial_weight_transfer_20261017.py) and [certificate](fresh_d17_partial_weight_transfer_20261017.json) save the complete weight vector and the literal supporting cases.

The finite-vertex budget is unchanged. At a cusp with link vertex count V, this factor's contribution becomes `120*V/17`, by the existing factor-fiber multiplicity. This local result does not assert an ambient dimension bound or any reduction of the whole order13 cycle rule.

## The complete affected support

The canonical affine-D16 diagram has spine0–1–…–12, leaves13,14 at0 and leaves15,16 at12. The previously removed singleton `(1,11)` remains at weight0. The new half-weight orbits are:

|Representative|Distance|All pairs|Selected contexts|
|---|---:|---|---|
|0,10|10|(0,10),(2,12)|1314,2529|
|0,11|11|(0,11),(1,12)|1311,2494|
|1,15|12|(1,15),(1,16),(11,13),(11,14)|1974,2367,2504,2511|
|2,15|11|(2,15),(2,16),(10,13),(10,14)|1978,2371,2535,2539|

The source obtains the context support from every emitted marked image in the frozen complete1,153-pattern mapping certificate, and verifies equality with these lists. Each orbit is invariant under every diagram automorphism. All other weights remain as in the singleton successor. In particular the last four-pair nonbasic orbit is retained at weight1.

Every listed context has the same literal19-root core: two L25 diagrams sharing their far node9. The overlap certificate supplies absence of additional cross-core edges whenever no retained short payment has already resolved the face. The proof below does not assume that an abstract finite pattern is the entire ideal factor; in these twelve cases the literal preferred component is already the full affine D16 diagram.

## Actual-kernel sealing applies uniformly

Number the first L25 by0 through9 and the second by10 through18 with shared9. The affine component F omits roots3 and13. The full19-root Gram has signature `(17,1,1)`, and its radical is an actual zero relation because its span contains a timelike vector. The radical coefficients at3 and13 are nonzero and have opposite signs.

For any additional actual supporting white root e, the positive affine kernel of F and positive semidefiniteness of the corresponding actual corner force e to be orthogonal to every F root. Pairing e with the full-core radical then forces its contacts with3 and13 either both to vanish or both to be nonzero. In the latter case the adjacent marks3,13 have a path3–e–13 of length2 and receive a retained unit.

Otherwise every additional supporting white is orthogonal to the entire19-root core. Each remaining corner is consequently a direct orthogonal sum of the displayed finite core corner and the positive-definite remaining white diagram. All displayed components remain whole components. This proves both the finiteness of these corners and the absence of indirect tails or bridges that could invalidate a finite coupon. No simple-cusp hypothesis is used.

The checker verifies this setup separately for all twelve marked contexts: the full signature and radical; positive affine kernel; actual white positive definiteness; all four strict edge principals; exact corner connected components and marked roles. These are48 literal corner decompositions. They share a geometry proof, not an assumption that the marked roles are interchangeable.

## The existing finite halves

In every affected context one adjacent finite corner has an existing half-payment between a whole small component of rank at most5 and a fork center or short fork leaf of a whole D11–D16 component. Representatives of the four groups are:

|Representative context|Whole small component|Whole large component and mark|Finite half|
|---|---|---|---:|
|1314|D5, arbitrary marked node|D12, fork center|1/2|
|1311|A4, endpoint|D13, fork center|1/2|
|1974 or2367|A4, endpoint|D13, short fork leaf|1/2|
|1978 or2371|D5, arbitrary marked node|D12, short fork leaf|1/2|

The reversed contexts give the same component roles. These payments are already included in the [whole-small-component finite budget](fresh_finite8_a17_d11_coupon_budget_20261017.md). The D5 mark need not be a leaf: the inherited rule permits any mark in the whole small component. No D12 short-leaf upgrade is used here; its old half suffices, even when the small component is D5 rather than A1.

The selected ideal corner now receives1/2, and this distinct finite corner receives1/2. Hence the quadrilateral still receives at least1. Four of the twelve contexts have an additional whole finite D17 unit; retaining that extra payment is harmless. The other eight have exactly the required two halves under the recorded rules.

## Transfer to the complete face-payment theorem

The original [complete finite8 reduction](fresh_finite8_ideal5_combined_payments_20261017.md) first resolves short marked paths, then classifies the remaining low/disjoint/overlapping minimal-core domains. The newly reduced pairs all have distance greater than9. Thus every short8-first alternative remains a unit; the elementary ideal triangle rule also retains all its units. The source checks those inclusions explicitly.

In the remaining complete overlap domain, any selected marked image whose weight changed is one of the twelve contexts above, by the complete frozen map support. If an extra white yields the length-two path, the quadrilateral receives a unit. Otherwise the preferred ideal corner contributes its half and the sealed finite corner contributes the existing half. Selected pairs outside these four orbits retain their previous weight and proof. The prior zero-weight singleton keeps its own separately validated two-finite-half replacement.

This argument preserves the inequalities “triangle total at least2” and “quadrilateral total at least1.” It does not reinterpret a half-weight as a binary paid unit, and it adds no cost to any finite vertex. The ideal-factor sum is linear in the weights, so the twelve half-unit reductions subtract exactly6 from126. A new scalar comparison must use this weight sum, not the number of pairs with nonzero weight. In particular the old concentration move eliminating repeated rank15 parts can change direction when the rank16 coefficient falls, and is not imported here.

## Verification scope

Generation25567 and a separate fresh saved replay16515 both terminated exit0/PASS for the complete saved-map support, all four symmetry orbits, all twelve sealed contexts and48 literal corner checks, and the exact120 total. Syntax, direct bindings and report links passed. Every input status and recorded direct binding is verified.

Two separate agents completed full source and written-mathematics reviews. Both passed the weighted composition: selected ideal halves combine with distinct sealed finite halves, short8 and triangle units remain unchanged, the previous zero singleton retains its own proof, and the actual factor-fiber cost is linear in the weights. The root review also checked that the inherited immediate whole-affine alternatives have order at most13, so no unlisted immediate D17 obligation is being changed. Neither review repeated the48-corner arithmetic or the1,153-pattern map census. The complete overlap and induced-map searches remain explicit frozen inputs; they are not repeated here. Source, data and report are final and frozen. This report makes no new scalar, dimension, external-review or formalization claim.
