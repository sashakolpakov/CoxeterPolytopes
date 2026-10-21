# The order13 improvement leaves seven ambient65 cusp obstructions

The completed [order13 B theorem](fresh_order13_affine_b_fork_projection_transfer_20261017.md)
changes the A/B/C/D factor costs to65/73/73/74 without changing any finite
weight. Its normalized rank12 envelope is therefore74/13. On a complete
domain of1,650 cusp profiles at ambient dimension65, this reduces the
number of negative margins from16 to7. The smallest new margin is
−35835/454597 at(10,10,12,16,16).

The [source](../scripts/fresh_order13_b73_scalar_frontier_20261017.py)
and [certificate](fresh_order13_b73_scalar_frontier_20261017.json)
record every profile and exact rational comparison. This is a scalar
diagnostic, not a new dimension theorem. The finite maximum remains487
against curvature480, independently of the remaining ideal deficits.

## A complete domain that survives further cost changes

Use the inherited convex-order reduction: spreading two rank summands
decreases normalized curvature. Two ranks below10 can be concentrated
until one vanishes or the other reaches10, preserving their linear
cost. Two ranks at least18 can be concentrated until one reaches17,
preserving their constant cost. Every move increases the bounded sum
of squared ranks, so it terminates.

The resulting complete domain has arbitrary repetitions of every rank
10 through17, at most one rank1 through9 and at most one rank at least18.
It has1,650 profiles of total rank64. It deliberately retains repeated
rank13 parts: the previous duplicate13 reduction becomes invalid in
one of the hypothetical improvements below. No such reduction is used
to justify this domain or its results.

For H(t)=product(1+t+...+t^r), put V=H(1) and
G=sum_j[t^j]H(t)*min(j,64−j). Each recorded margin is

    (5G−64V)/(2V) − sum_i c(r_i).

The source constructs every integer polynomial directly and uses exact
rational costs, with the proved new c(12)=74/13. The seven negative
profiles are:

|Ranks|Margin|
|---|---:|
|(1,10,10,11,16,16)|−3563/209814|
|(2,10,10,10,16,16)|−1108/67881|
|(10,10,10,10,10,14)|−14986/219615|
|(10,10,12,16,16)|−35835/454597|
|(10,10,14,14,16)|−1376/30855|
|(10,11,11,16,16)|−3869/114444|
|(16,16,16,16)|−96/4913|

These are comparison obstructions, not geometrically realized cusps.
Six have no rank12 part, so further order13 improvements alone cannot
resolve the ideal comparison.

## Exact simultaneous targets

Starting with the proved order13 maximum74, suppose one could also
prove the maximum order15 cost at most94 instead of96, and the maximum
order17 cost at most119 instead of120. Then all1,650 margins would be
positive, with minimum61/1936 at(9,10,10,10,10,15).

These are maxima over every applicable affine family. Improving only
the B family at order15 is insufficient when C or D still has a larger
cost. With the order17 reduction and only a one-unit order15 reduction,
the single profile(10,10,10,10,10,14) remains negative at−23/14641.

The hypothetical changes leave the low linear and high constant costs
used in this domain proof unchanged. Their geometric validity and any
additional finite coupons remain obligations. Even if these ideal
targets were proved without new finite expense, the finite deficit7
at ambient65 would still prevent a new global bound by this comparison.

## Verification scope

Generationfb4876 and a separate fresh saved replay026822 terminated
exit0/PASS on all1,650 profiles and all conditional scenario rows.
The source imports no project arithmetic. Its actual local input is
the frozen order13 B theorem; hypothetical order15/17 drops are clearly
separated. A separate full source/report read passed the complete domain,
direction of concentration, all-family targets and finite-deficit scope.

An independent standard-library reconstruction also passed all1,650
actual rows and all three hypothetical scenario tables. It enumerates
optional low/high parts before multiplicities of10 through17, uses
sliding-window convolution, and independently computes G through an
inclusion-exclusion cumulative-composition identity. It imports no
producer arithmetic. Every domain element, polynomial, V, G and margin
agrees, including the old16/current7 negatives. That separate run
terminatedfcd4f6 with exit0/PASS.

All11 direct source bindings, three local report links and Python syntax
pass. Source, certificate and report are final/frozen. No new ambient
bound, external review, formalization or example is asserted.
