# Exact scalar frontier after the order14 B reduction

The new [order14 B transfer](fresh_order14_affine_b_fork_projection_transfer_20261017.md)
changes the A/B/C/D costs from70/86/84/82 to70/84/84/82. Thus the normalized
affine **rank13** coefficient falls from43/7 to6. This resolves the former
worst formal cusp profile in ambient66, but ten negative ideal margins
remain. There is no new ambient dimension bound.

The [source](../scripts/fresh_order14_b84_scalar_frontier_20261017.py)
and [certificate](fresh_order14_b84_scalar_frontier_20261017.json) evaluate
the complete old, corrected new and common expanded concentration domains
with exact rational arithmetic. All other ideal coefficients are those
of the completed ceiling66.

## Correct rank convention and required improvement

Five rank13 factors have total link rank65. Each factor has14 nodes.
Their old limiting type was B with cost86; the order13 costs65/75/73/74
instead belong to affine rank12 and are not the coefficient in this
profile.

The old margin of(13,13,13,13,13) is-1565/10976. A total pair-cost decrease
delta in each worst order14 factor adds5*delta/14. Resolving this single
profile requires delta>313/784. The proved two-unit reduction gives new
margin6275/10976. B and C now tie at cost84; A70 and D82 remain below them.

## The concentration domain must change

The inherited convex-order theorem decreases normalized curvature under
concentration(a,b) to(a-1,b+1), for1<=a<=b. A valid margin reduction also
requires that total cost be preserved or increased.

Previously two rank12 factors could become11 and13, with cost increase
23/1092. With the new coefficient this increment is-19/156. The old move
is therefore invalid. Indeed the move already fails for an order14 cost
drop greater than23/78, which is smaller than the313/784 needed to resolve
the former worst profile. Reusing the old terminal list would miss a
real proof obligation.

The new valid middle move is13+13 to12+14, whose cost increase is11/65.
Low parts still concentrate only up to10, preserving their linear cost;
two high parts at least18 still concentrate until one is17, preserving
their constant cost. Every move strictly increases the bounded integer
sum of squared ranks and therefore terminates.

The old terminal domain retains arbitrary copies of10,11,13,14,15,16,17
and at most one12. The corrected new domain retains arbitrary copies of
10,11,12,14,15,16,17 and at most one13. Both have at most one low part1--9
and at most one high part at least18. These are complete reductions for
their respective costs. A third, common expanded domain permits arbitrary
copies of every rank10--17, applying only the low/high moves. It provides
a common exact comparison without relying on either middle move.

## Complete old and new comparisons at ambient66

|Domain and cost rule|Profiles|Negative margins|Minimum margin|Minimizing ranks|
|---|---:|---:|---|---|
|Old terminal domain, old costs|1,625|23|-1565/10976|(13,13,13,13,13)|
|New terminal domain, new costs|1,659|10|-31895/248832|(10,11,11,11,11,11)|
|Common expanded domain, old costs|1,796|28|-1565/10976|(13,13,13,13,13)|
|Common expanded domain, new costs|1,796|10|-31895/248832|(10,11,11,11,11,11)|

There are no zero margins. All ten remaining negative profiles have no
rank13 factor, so their numerical margins were already negative before
the improvement. Eight were explicitly present in the old terminal list:

|Ranks|Unchanged negative margin|
|---|---|
|(10,11,11,11,11,11)|-31895/248832|
|(11,11,11,16,16)|-1055/10404|
|(10,10,11,11,11,12)|-3535/41184|
|(10,11,12,16,16)|-2330/41327|
|(9,11,11,11,11,12)|-12799/539136|
|(1,10,10,11,11,11,11)|-449/50688|
|(1,10,11,11,16,16)|-103/28611|
|(11,11,12,15,16)|-863/509184|

Two have repeated rank12 and are newly admitted to the terminal list:

|Ranks|Unchanged margin|Their old concentration image|Old image margin|
|---|---|---|---|
|(10,10,10,11,12,12)|-205/4719|(10,10,10,11,11,13)|-481/5544|
|(11,12,12,14,16)|-11/1690|(11,11,13,14,16)|-24277/514080|

The checker verifies those old images and the required margin directions.
These are newly necessary terminal checks, not negative margins created
by lowering a cost. For every profile in the common domain the new margin
equals the old margin plus its number of rank13 factors divided by7.

## A quantitative next obligation, without a local proof

The order12 maximum remains65, shared by B and C; A and D cost60 and62.
Changing that maximum would change c(11). The new worst profile alone
would require an order12 cost drop greater than6379/20736. That is not
the strongest remaining requirement.

Across the complete corrected domain, the exact required uniform drop
in the order12 maximum is greater than27960/41327, controlled by
(10,11,12,16,16). Every remaining negative profile contains rank11.
A hypothetical decrease65 to64 makes all1,659 corrected-domain margins
positive, with minimum13367/495924 at that controlling profile. The
concentration proof remains valid because its low and high moves and
the duplicate13 increment are unaffected by c(11).

This is only a conditional ideal comparison. Both B and C would need
geometrically valid reductions to at most64, with full payment provenance
and any new finite coupons counted. This package provides no such second
local theorem. The completed finite rule also still has the separate
D66 cost496 versus curvature495 obstruction. Neither obstruction can be
discarded on the basis of these scalar calculations.

## Verification scope

Final generation042067 and separate fresh replay668fe2 both terminated
exit0/PASS on all three complete domain enumerations, exact profile
curvatures, old/new margins and concentration increments, all ten remaining
negative profiles and their old-domain provenance, and the conditional
next-order12 comparison. A separate implementation check uses the frozen
global66 reviewer's optional-first multiplicity enumeration and its
sliding-window/CDF curvature arithmetic: all1,796 common-domain profiles
and both terminal subdomains agree exactly. This implementation check
was run by the producer; it is not a separate mathematical reviewer.

Root completed a full source/report mathematical read. That review passed
the corrected repeated-rank12 and at-most-one-rank13 concentration domain,
the common1,796-profile comparison, all ten remaining negative profiles
and their provenance, and the conditional order12 threshold27960/41327
and positive maximum64 comparison. Root did not repeat the arithmetic.
All eight direct bindings, Python AST and three local report links pass.
The local B14 report is final/frozen and its final hash is bound here.
Source, certificate and report are final/frozen. The all-rank convex-order
theorem and existing ideal cost envelopes remain frozen inputs.

No hyperbolic realization of a formal profile, new ambient bound, external
refereeing or formal verification is asserted.

    python3 scripts/fresh_order14_b84_scalar_frontier_20261017.py
    python3 scripts/fresh_order14_b84_scalar_frontier_20261017.py --verify
