# Deleting the affine D12 center half with full finite coupons

Relative to the completed dimension62 rule, the unique fork-center pair
in a whole affine D12 factor of order13 may be reduced from1/2 to0,
provided precisely these whole finite coupons are upgraded:

|Separate whole components and marked roles|Old weight|New weight|Old pair total|New pair total|Increment|
|---|---:|---:|---:|---:|---:|
|A2 either endpoint / D10 fork center|1/2|1|1|2|1|
|A3 unique middle / D10 fork center|1/2|1|1/2|1|1/2|

Thus the ideal D12 total falls from147/2 to73. The four order13 family
totals become65,73,73,73, with all-family normalized envelope73/13.
For a finite forest with multiplicities a2,a3,d10, the added charge is

    a2*d10 +(a3*d10)/2.

This local theorem does not bound that enlarged finite charge globally,
lower an ambient dimension, or assert that the ideal improvement alone
resolves the current frontier.

The [source](../scripts/fresh_order13_d_center_pair_zero_transfer_20261018.py)
and [certificate](fresh_order13_d_center_pair_zero_transfer_20261018.json)
compose the frozen
[universal D13 half theorem](fresh_order13_d_center_pair_half_transfer_20261018.md)
with the [completed dimension62 rule](fresh_sparse_coupon_order11_13_dimension62_20261018.md).
The arbitrary-label and arbitrary-spectator proof is inherited in full,
not replaced by the finite list of literal corner roles.

## Every branch has a distinct full unit

Use the preceding theorem's actual quadrilateral star W+a+b+c+d, with
the selected ideal factor F containing a,b. In its canonical numbering,
F has spine0 through8 and two short leaves at each end; the changed
marks are a=0,b=8. All roots have physical norm2. The actual varying
factor is whole, its other supporting components form the positive
definite spectator Z, and both opposite-root heights are positive.
All retained weights below are evaluated in the simultaneous new rule.

If an adjacent black product is negative, or c,d have a common
spectator contact, the preceding proof gives a **distinct** actual
corner a path of length1 or2. It still carries a unit. No original
ideal payment is used.

For the immediate end4 or next3 profiles, the old positive-kernel
argument seals a whole order11 affine B10 or D10 factor at another
corner. Its marks are respectively center/terminal distance8 or
center/leaf distance7. These are retained units under the completed62
rule, not any of the order11 terminal/fork-leaf halves or the deleted
C11 endpoint pair. The source compares all four reflected immediate
records with their inherited unit payments.

After these cases the exact R=8 determinant and native profile bounds
leave only the following possibilities. These are precisely the full
analytic branch partition from the half theorem; no additional rank
or label cutoff is introduced.

* Two nonextreme ends are impossible: their cost sum is at most6,
  giving a strictly positive determinant gap at least1/4.
* One extreme and one nonextreme end force the extreme side's entire
  spectator contact vector to vanish. The weighted spectator energy
  is at most1, while any nonzero extreme coordinate would cost at
  least2 by Gram-metric Cauchy. In the adjacent corner the extreme
  root and its old fork leaf form a whole A2, separate from the whole
  D10 containing the far mark at its fork center. The upgraded coupon
  now gives a **distinct unit**, replacing the old coupon half plus
  original half. No condition on the other spectator vector is needed.
* If both ends are extreme, the unscaled spectator energy equals1.
  If one support is empty, the same whole A2/D10 corner gives a unit.
  If both are nonempty, the bound inherited energy-one classification
  gives a mixed whole A3 or D4–D7, yielding a distinct path4 unit, or
  one contact on each of two separate whole A1 components. In the
  latter case the selected adjacent corner has a whole A3 with its
  new mark at the unique middle, separate from whole D10 marked at
  its fork center. The upgraded A3-middle coupon gives a distinct
  full unit.

Whole-component scope is unchanged by upgrading weights. In the
one-sided branch the chosen opposite root has no Z contact at all;
the other opposite root is absent from that corner. In the two-A1
branch the equality classification forces all other spectator contacts
to vanish and identifies those A1s as entire components. Consequently
neither small A2/A3 nor far D10 can acquire an undisclosed supporting
tail or bridge. Only the middle of A3 is charged; its endpoints never
need an upgrade.

The four immediate affine records, eight finite coupon records and
five mixed energy-one records are checked as17 literal interfaces.
Each now has a distinct payment1 and the original ideal contribution0.
Those records corroborate the new roles; completeness comes from the
preceding universal proof and the analytic branch partition just given.
In particular this theorem is not the invalid inference that an old
half-transfer automatically permits deleting its original half.

## Simultaneous coverage and finite cost

For every actual quadrilateral, first test whether any actual corner
contains the newly zeroed D13 center pair. If one does, the proof above
supplies a distinct full unit under the new rule. If none does, every
actual ideal weight agrees with the completed62 rule, and finite
weights only increased. The preceding whole-face payment at least1
therefore transfers directly. This also covers older proofs which may
have selected a D13 corner: the case distinction is about all actual
corners, not the provenance of an earlier selected witness.

The new outputs are short paths of length at most4, the retained
order11 center units, and the exact whole A2/A3-middle/D10-center
coupons. They avoid sparse A16/D16/D17 long deletions, deep11 finite8
zeros, and the already deleted C11/B12 pairs. No changed branch invokes
the unchanged-corner case recursively. The order13 distance8 center
pair is absent from the inherited finite7/ideal5 triangle proof, so
triangle payment2 remains unchanged. Every other pair weight and every
other finite coupon is retained.

The finite increments are literal counts of invariant marked roles:
A2 has two interchangeable endpoints and D10 one fork center, giving
two upgraded halves; A3 has one middle and D10 one fork center, giving
one upgraded half. They are distinct from the old D10 short-leaf
coupons, which keep their old values. The
[preceding finite budget](fresh_sparse_coupon_rank63_forest_budget_20261018.md)
counted center totals1 and1/2 only. It cannot be reused as an enlarged
budget without adding the displayed increments and recomputing a
complete forest bound. That is a separate task.

## Verification scope

Final generation61a495 and a separate fresh replay892cd1 terminate
exit0/PASS. Root completed a full source/report mathematical read:
every upgraded whole-component unit genuinely dispenses with the
original ideal half; the17 literal interfaces corroborate the complete
inherited analytic branch partition; changed-actual-D13-first priority
is noncircular; and the exact increments1 and1/2 concern center roles
distinct from the old short-leaf coupons. That review did not rerun
matrices, arithmetic or classification. The source imports only
standard-library arithmetic and checks frozen
bindings, the sole ideal weight change, all four family totals, exact
marked-role multiplicities, every inherited immediate/coupon/mixed
interface, and the new finite increment. It does not rerun the preceding
profile matrices, spectator classification, historical censuses or any
finite-forest optimization. Source, data and report are final/frozen.
No global dimension claim, external review
or proof-assistant formalization is made.

    python3 scripts/fresh_order13_d_center_pair_zero_transfer_20261018.py
    python3 scripts/fresh_order13_d_center_pair_zero_transfer_20261018.py --verify
