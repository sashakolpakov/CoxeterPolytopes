# What an affine order12 cost of65 would suffice to prove

If every affine factor of order12 can be paid with cost at most65, with
the other current ideal costs unchanged, all ideal comparisons are
strictly positive in ambient dimensions68 and69. This is a conditional
calculation. The unrestricted ceiling remains69; no ceiling68 is asserted.

The [source](../scripts/fresh_rank11_cost65_scalar_frontier_20261017.py)
and [certificate](fresh_rank11_cost65_scalar_frontier_20261017.json)
use the independently reconstructed scalar identities from the
[completed69 review](fresh_ideal5_weighted_dimension69_20261017_review.md).
Only the rank11 cost changes, from11/2 to65/12. The complete A/B/C local
transfers and any resulting finite-coupon budget are still required.

## The concentration domain must change

The old reduction concentrated low ranks up to11 using the linear cost
r/2. With rank11 now cheaper by1/12, that cost-preservation argument no
longer applies. Instead concentrate ranks at most10 only up to10. Their
costs stay linear, and the inherited convex-order theorem decreases
normalized curvature under concentration. This leaves arbitrarily many
rank10 parts and at most one lower part.

Rank11 parts are retained without an elimination rule. Two rank12 parts
can still become11 and13: the cost increase is now

    65/12 + 43/7 - 2*(75/13) = 23/1092 > 0.

The unchanged constant-cost high-rank rule leaves arbitrary copies of17
and at most one rank at least18. Hence a sufficient terminal domain has
arbitrary copies of10,11,13,14,15,16,17; at most one12; at most one part
between1 and9; and at most one part at least18. Every nontrivial move
strictly increases the bounded integer sum of squared ranks, proving
termination. The curvature decreases and cost is preserved or increased,
so positivity on this complete terminal domain implies positivity before
the reductions.

The checker enumerates this new domain explicitly. For every row it
computes the link polynomial, vertex count and curvature by the inherited
separately implemented convolution, cumulative-composition and absolute-
deviation identities, then replaces the cost exactly. No old terminal
profile list is reused.

|Ambient dimension|Terminal profiles|Minimum ideal margin|Minimizing ranks|
|---:|---:|---|---|
|68|1,906|411499/3234816|(11,11,11,11,11,12)|
|69|2,060|3641/13608|(11,11,11,11,11,13)|

## Finite and local obligations remain distinct

The [revised exact finite rule](fresh_a17_endpoint_d11_upgrade_forest_budget_20261017.md)
removes the A17 endpoint unit and includes the full A1/D11-short-leaf
upgrade. It has maximum544 in both ambient ranks68 and69. The finite
curvature coefficients are527 and544 respectively. Thus the finite
comparison still fails by17 at68, whereas it has zero slack at69.
At69 the strict global constant340 and the positive ideal margins would
give a contradiction if all local coverage and final finite costs were
established. The completed ceiling69 would supply the higher-dimensional
tail, yielding ceiling68 under those outstanding hypotheses.

The newly completed D11 transfer provides affine-D order12 cost62.
It does not reduce the maximum over A/B/C/D by itself. The other local
branches can introduce additional finite coupons, so the existing544
finite maximum must not be silently reused for a larger rule. The next
work is to close those local branches and verify the actual combined
finite rule.

Generation and fresh separate-process saved replay both terminated
exit0/PASS. A separate full source and written-proof review passed the
changed concentration domain, duplicate12 cost gain, termination,
optional-first enumeration, exact cost replacement and finite/ideal
separation. That review did not repeat the1,906/2,060 scalar calculations.
Source, certificate and report are final/frozen. No new ideal envelope
or ambient theorem is promoted by this conditional package.
