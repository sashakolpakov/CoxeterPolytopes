# Combined rank15/order15 frontier: finite413 and49 negative cusp profiles

The new finite-rank15 and affine-order15 local theorems compose. With all
compensating finite payments included, the complete ambient60 finite
maximum remains **413**, against curvature405. The number of negative
terminal cusp margins falls from77 to **49**, among all2,113 profiles;
none has zero margin. These improvements do not yet exclude dimension60.
The completed unrestricted dimension bound remains60, and no new
finite-volume example is established.

The [source](../scripts/fresh_rank15_order15_ambient60_frontier_20261018.py)
and [certificate](fresh_rank15_order15_ambient60_frontier_20261018.json)
bind the frozen [finite exchange](fresh_rank15_sparse_payment_exchange_20261018.md),
[order15 B/D theorem](fresh_order15_bd_terminal_coupon_transfer_20261018.md),
[depth10 local theorem](fresh_depth10_finite7_zero_transfer_20261018.md)
and [preceding complete frontier](fresh_ambient60_payment_frontier_20261018.md).
They leave the [completed dimension60 proof](fresh_order11_order16_dimension60_20261018.md)
unchanged.

## One simultaneous local rule

The finite exchange gives A15=101, B15=97, D15=102 and A17=116. The last
cost includes the restored whole A17 pair(1,15) at weight1. All other
same-component finite costs follow the current rule, including the
depth10 distance7 deletions. The affine-order15 costs become
A75/B92/C91/D92, with maximum92 and normalized envelope92/15.

Exactly three actual finite component-pair costs increase:

|Whole component pair and changed marks|Added cost|New total actual cost|
|---|---:|---:|
|A2, either node / D13, short fork leaf|2|5|
|A4, endpoint / D11, short fork leaf|2|8|
|A5, second-from-end node / D11, short fork leaf|2|19/2|

Each increase upgrades four half payments to units. Every other mark
keeps its previous weight, including fork-center marks. These costs are
counted together with the A17 increase.

For a changed order15 D corner, the universal local proof supplies a
distinct whole A4/D11 or A5/D11 unit, or its independent short-path
escape. The changed B corner has an impossible forced completion or
the same independent short-path alternative. Neither proof consumes
a same-component finite long-pair unit. Thus the new rank15 finite
decreases cannot invalidate the new ideal transfers.

Conversely, the finite exchange uses unchanged whole A16 halves, the
restored whole A17 unit and the upgraded whole A2/D13 unit, together
with the independently established short8 route. Its depth8-to-short7
and depth10-to-short6 replacements retain the necessary protected
outputs. The new ideal deletions have distance10 and do not remove
those short-path outputs. Current earlier ideal-corner priorities and
the complete selected/nonselected and later literal interfaces retain
their established scope. Hence the two local theorems apply under a
single rule, with every triangle paid by at least2 and quadrilateral
by at least1. This composition uses their universal proofs, not merely
their finite lists of literal configurations.

## Complete finite count with new residual tables

The preceding representative rule already allowed5 for A2/D13 in
order to dominate A2/A13. The new actual A2/D13 cost fits that allowance.
The two D11 additions genuinely raise the upper rule: they increase
the rank4 and rank5 residual item values by4 times the number of fixed
D11 factors in doubled units. The source checks both ends of every
typed component-pair domination simultaneously.

D15 at102 dominates the new A15/B15 values101/97. A17 remains a
separate fixed species at116. The complete fixed species are
B9,D10,...,D16,A17. There are1,044 multiplicity vectors and792 residual
interaction keys. The residual ranks are1..8 and17..60. All fixed-fixed
and residual-residual upper cross costs are zero, and every residual
item is checked directly against its base plus all fixed interactions.

Generation constructs fresh item-first unbounded-knapsack tables for
all792 keys. Replay does not repeat that optimization. It checks every
possible last-item upper inequality and an attaining predecessor at
each rank, proving the tables by induction. This covers1,142,064 upper
inequalities,47,520 attaining predecessors and54,288 item identities.
Every fixed forest is then combined with its entire residual optimum.

The complete maximum is413. The formal forest2D17+D26 attains it and
contains neither a newly charged component pair nor a deep distance7
pair. The actual witness costs recorded for other reconstructed forests
subtract the depth10 deletions; the representative upper is allowed to
ignore those nonnegative savings. Thus attainment is for the current
actual payment rule, not just for a conservative upper formula. No
geometric realization of these formal forests is asserted.

There are15 representative alternatives above405. In particular, four
rank15 factors now have costs404 for all-A,388 for all-B and408 for
all-D. The finite obstruction has narrowed, but the worst deficit
413−405=8 remains.

## Complete ideal comparison and remaining targets

The low9-only concentration theorem still applies: costs through rank9
remain r/2, and no high-rank factor is concentrated. Independently
enumerating every terminal partition of rank59 gives2,113 profiles,
with at most one part1 through8 and all other parts at least9.

For each profile, generation and replay reconstruct the exact integer
convolution

    H(t) = product_i (1+t+...+t^(r_i)).

The coefficient sum, symmetry, weighted curvature moment and independent
absolute-deviation/CDF identities are checked against the prior exact
row. The only normalized cost change is

    saving = 2 * number of factors of rank14 / 15.

This includes the actual factor-fiber multiplicity. All other envelopes
and family tail formulas are unchanged. Exact recomputation resolves28
previously negative profiles and leaves49 negative, with none zero.
The minimum normalized margin is still−1823/11560 at(9,9,9,16,16).
The minimum terminal absolute margin is now−123905 at(9,9,9,9,11,12).
Negative terminal absolute margins are not universal absolute bounds
for unreduced cusps.

The remaining part of the preceding sufficient envelope scenario is:

|Affine order|Current maximum|Conditional target|Required saving|
|---|---:|---:|---:|
|12|63|125/2|1/2|
|13|73|145/2|1/2|
|14|83|331/4|1/4|
|16|105|415/4|5/4|
|17|116|229/2|3/2|

Order15 already meets its former target92. Applying these further
conditional savings to all2,113 profiles gives normalized minimum
73/104000 and terminal absolute minimum54. These are exact sufficient
scalar targets; the corresponding local theorems and any added finite
costs remain unproved here. They are not an optimality claim.

The global two-extrema correction remains4500 and is counted once.
The finite deficit and49 negative ideal profiles prevent the usual
pointwise curvature contradiction, so no bound59 is claimed.

## Verification scope

Generation session92271, final chunk71fb95, terminated exit0/PASS.
Fresh replay session4000, final chunk589cdf, terminated exit0/PASS.
An independent internal agent completed a full mathematical source/report
composition review, including the inherited local interfaces, typed costs,
complete finite domain and ideal comparison: PASS. That review did not
execute a second geometric census, Bellman calculation or convolution.
The binder validates the frozen local inputs and exact typed increments,
all finite domination and Bellman certificates, all ideal convolutions
and scalar comparisons, and the single global constant. It does not
rerun a geometric census or prove the local theorems anew. The packet is
now frozen. No external refereeing or formal verification is claimed.

    python3 scripts/fresh_rank15_order15_ambient60_frontier_20261018.py
    python3 scripts/fresh_rank15_order15_ambient60_frontier_20261018.py --verify
