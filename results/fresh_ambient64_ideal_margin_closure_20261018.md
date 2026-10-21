# Positive ideal-vertex margins at ambient dimension64

The new order11 and order17 factor envelopes give **positive margins for
every possible ideal-vertex rank profile at ambient dimension64**. The
minimum normalized margin is2385/146432; the minimum absolute surplus is
240/23. The complete finite coupon budget and simultaneous global payment
composition are separate obligations. This scalar packet alone asserts
no new dimension bound.

The [source](../scripts/fresh_ambient64_ideal_margin_closure_20261018.py)
and [certificate](fresh_ambient64_ideal_margin_closure_20261018.json)
independently reconstruct the complete corrected profile domain. They use
rank multiplicities and sliding-window polynomial multiplication, rather
than the earlier diagnostic's sorted-rank recursion and direct convolution.
No project arithmetic helper is imported.

## Exact factor envelopes

The [affine11 cycle transfer](fresh_cycle11_half_pair_transfer_20261018.md)
and [coupled order11 tree transfers](fresh_order11_bcd_half_pair_transfer_20261018.md)
give total costs A=99/2, B=54, C=109/2, D=53. Hence the maximum divided by
factor order is109/22 at affine rank10. The tree transfers add a finite
half coupon between a whole A1 and the doubled terminal of a separate
whole B10.

The [order17 refinement](fresh_order17_e6_half_coupon_transfer_20261018.md)
gives A85/B111/C108/D118, with envelope118/17 at affine rank16. Its finite
E6/D11 coupon increases from a quarter to a half on the specified two-by-two
marked roles. These coupon costs must be included on the finite side;
they do not enter an ideal factor's pair total.

All other ideal factors retain the [previous completed envelopes](fresh_ambient65_ideal_margin_closure_20261017.md).
For affine rank r the complete normalized upper envelope is

|Rank r|c(r)|
|---|---:|
|1 through9|r/2|
|10|109/22|
|11|16/3|
|12|74/13|
|13|6|
|14|94/15|
|15|107/16|
|16|118/17|
|17 and larger|163/23|

Exceptional affine factors have rank at most8 and retain their earlier
envelopes. The factor-fiber multiplicity remains V/(r+1) for each marked
pair, so the cost at an actual cusp is at most V times the sum of c(r).

## Complete rank-profile reduction

Use the inherited [all-rank convex-order concentration theorem](fresh_distance14_unrestricted_dimension156_20261017.md).
For two positive ranks below9, transfer rank from the smaller to the
larger until the smaller vanishes or the larger reaches9. For two ranks
at least18, transfer until the smaller reaches17. Every transfer decreases
normalized curvature and preserves the sum of costs: the respective cost
ranges are linear and constant. It strictly increases the bounded sum
of squared ranks, so the procedure terminates.

The total affine rank is63. A terminal profile therefore has arbitrary
multiplicities of9 through17, at most one rank1 through8, and at most one
rank at least18. There is no additional restriction on repeated middle
ranks. The source enumerates each of the2,419 profiles exactly once by
choosing the nine middle-rank multiplicities and then the two possible
extreme factors.

The older1,513-profile domain, based on a linear envelope through rank10,
would not be complete after changing the rank10 cost. The larger domain
is essential here.

## Exact curvature calculation and absolute surplus

For a profile(r_i), set

    H(t)=product_i(1+t+...+t^r_i),  V=H(1),
    G=sum_j [t^j]H(t)*min(j,63-j),
    margin=(5G-63V)/(2V)-sum_i c(r_i).

Every polynomial, V, G, curvature, charge and margin is computed exactly.
Every row agrees with the corresponding scalar data in the
[earlier conditional frontier](fresh_ambient64_payment_frontier_20261018.md).
The new normalized margins agree with its simultaneous cost drops of
1/2 at factor order11 and1 at factor order17. This comparison uses all
2,419 rows, not only the extremal examples.

The unique minimum normalized margin is2385/146432, at ranks
(10,11,12,15,15). The unique minimum absolute margin V times the normalized
margin is240/23, at the single-factor profile(63). All rows are positive.

The absolute lower bound also holds before concentration. For a transfer
delta from a<=b to the larger factor,

    (a+1)(b+1)-(a-delta+1)(b+delta+1)
      =delta*(b-a+delta)>0.

Thus V decreases at every transfer, including disappearance of a factor,
whose contribution then becomes1. Normalized margin also decreases and
its terminal value is positive. Their product consequently decreases.
The terminal absolute minimum240/23 therefore bounds every actual cusp,
without a simple-cusp, arithmetic or number-of-factors hypothesis.

The [shape-curvature extremum estimate](fresh_cusp_link_shape_curvature_20261017.md),
with the quantitative strengthening recorded in the previous ideal-margin
proof, supplies constant5*floor(64²/4)=5120. The finite side must have
maximum at most464 to combine this positive constant and ideal surplus
into a dimension64 contradiction. This packet does not establish that
finite hypothesis or the full simultaneous face-payment theorem.

## Verification scope

Final source-bound generation621472 and separate fresh replay902f1b both
terminated exit0/PASS on all2,419 profiles, both unique extrema and every
field comparison with the independently implemented earlier conditional
calculation. The certificate binds the completed local proof packets.
Its rank-multiplicity recursion and sliding-window polynomial arithmetic
are separate implementations of those scalar checks.

A full separate source/report mathematical read, including the inherited
convex-order and shape-curvature arguments, passed the complete threshold9/17
domain, preserved costs, V-product monotonicity, positive absolute surplus,
factor envelopes and constant5120. That reviewer did not rerun the2,419
rows. Source, data and report are final/frozen. All review is internal;
external refereeing and proof-assistant verification remain absent.
The geometric factor-fiber theorem and all-rank concentration theorem
remain explicitly inherited inputs; scalar arithmetic does not prove
those geometric statements.
