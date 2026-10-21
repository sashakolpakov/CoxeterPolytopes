# Exact positive ideal margins at ambient dimension63

For the specified order11/12/13 total-cost envelopes54,63,147/2, every
ideal rank profile at ambient dimension63 has positive margin. The unique
normalized minimum is3505/201344 and the unique absolute minimum is139/46.
This is an exact scalar conclusion. Geometric validity of the new local
weights, their finite coupon costs and a simultaneous global theorem are
separate obligations; this packet asserts no new dimension bound.

The [source](../scripts/fresh_ambient63_ideal_margin_closure_20261018.py)
and [certificate](fresh_ambient63_ideal_margin_closure_20261018.json)
reconstruct every profile by a multiplicity recursion and every polynomial
by sliding-window sums. They import no project arithmetic module and agree
with all scalar rows of the independent
[ambient63 frontier](fresh_ambient63_payment_frontier_20261018.md), whose
implementation uses sorted-rank recursion and direct convolution.

## Specified envelopes and complete domain

Let r denote affine rank, so the factor has r+1 nodes. The normalized
envelope used in the calculation is

|r|c(r)|
|---|---:|
|1 through9|r/2|
|10|54/11|
|11|21/4|
|12|147/26|
|13|6|
|14|94/15|
|15|107/16|
|16|118/17|
|17 and larger|163/23|

Only the three middle values differ from the
[completed dimension64 ideal calculation](fresh_ambient64_ideal_margin_closure_20261018.md).
They require total pair-cost decreases1/2,1,1/2 at orders11,12,13
respectively. Finite coupons do not enter these ideal costs, but must be
included separately at finite vertices.

The [all-rank concentration theorem](fresh_distance14_unrestricted_dimension156_20261017.md)
transfers rank from a smaller to a larger factor while decreasing normalized
curvature. For two positive ranks below9, stop when one vanishes or the
larger reaches9. For two ranks at least18, stop when the smaller reaches17.
These transfers preserve the sum of c(r), since the respective ranges are
linear and constant. They strictly increase the bounded sum of squared
ranks, so they terminate.

Every terminal rank62 profile consequently has arbitrary multiplicities
of9 through17, at most one rank1 through8, and at most one rank at least18.
The source chooses the nine middle multiplicities, then the remaining
low/high factors, generating each of the2,204 profiles once. The ranks
always sum to62. Changing only the three listed middle costs does not
alter either concentration boundary or the completeness of this domain.

## Exact calculation and absolute bound

For a profile(r_i), define

    H(t)=product_i(1+t+...+t^r_i), V=H(1),
    G=sum_j [t^j]H(t)*min(j,62−j),
    margin=(5G−62V)/(2V)−sum_i c(r_i).

The factor-fiber formula gives at most V times the cost sum at an actual
cusp. All coefficients, V, G, curvature, cost and margins here are exact.
Every polynomial and scalar field matches the frozen frontier, and each
new margin matches its three-order conditional scenario. The gain from
the old normalized margin is exactly

    count(rank10)/22 + count(rank11)/12 + count(rank12)/26.

All2,204 margins are positive. The unique normalized minimum3505/201344
occurs at(10,10,12,15,15). The unique absolute minimum139/46 occurs at(62).

For any concentration transfer delta from a<=b,

    (a+1)(b+1)−(a−delta+1)(b+delta+1)
      =delta*(b−a+delta)>0.

Thus V decreases together with normalized margin. Since the terminal
margin is positive, their product decreases as well. The terminal absolute
minimum therefore gives surplus at least139/46 for every actual cusp
profile, including those with many small or large factors. No simple-cusp
or arithmetic assumption is needed for this reduction.

The inherited [shape-curvature estimate](fresh_cusp_link_shape_curvature_20261017.md)
supplies the additive constant5*floor(63²/4)=4960. To use it in a new
dimension theorem requires a complete compatible payment rule with all
three specified ideal envelopes and finite maximum at most899/2. This
scalar calculation does not supply those geometric or finite hypotheses.

## Verification scope

Generationd7483c and separate fresh replay14112e both terminated exit0/PASS
on every profile, all independently implemented field comparisons and both
unique extrema. A separate full mathematical/source read passed the complete
multiplicity domain, sliding-window implementation, concentration boundaries,
cost preservation, V monotonicity, positive absolute surplus, constant4960
and explicit conditional geometric scope. The reviewer did not rerun the
2,204-row arithmetic.
The actual factor-fiber, concentration and shape-curvature arguments are
inherited written mathematical inputs, not consequences of this numerical
table. No external refereeing or proof-assistant verification is claimed.
Source, certificate and report are final/frozen after these checks.
