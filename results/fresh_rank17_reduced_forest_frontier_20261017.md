# What reduced finite rank17 costs would suffice to prove

Reducing whole finite A17/B17/D17 costs to114/105/127, while preserving
all current coupons and face coverage, would make the existing comparison
exclude ambient dimensions67 and68. With the completed ceiling68, that
would give ceiling66. This package proves the exact finite counting under
explicit assumed rank17 costs. Geometric coverage and the final global
composition remain separate obligations; it makes no new dimension claim.

The [source](../scripts/fresh_rank17_reduced_forest_frontier_20261017.py)
and [certificate](fresh_rank17_reduced_forest_frontier_20261017.json)
evaluate two abstract rules. All costs other than the specified rank17
bases, and all finite cross-component coupons, are unchanged from the
[completed68 proof](fresh_ideal5_weighted_dimension68_20261017.md).

|Rule|A17 base|D17 base|Required B17 upper bound|
|---|---:|---:|---:|
|Mild reduction|131|131|131|
|Full proposed reduction|114|127|127|

The proposed geometric B17 cost105 would satisfy the second upper bound.
Neither rule is inferred merely from an absence in a partial selected-map
list. The geometry must independently establish full face coverage.

## Complete abstract finite-forest count

Keep the inherited finite-type domination in ranks other than17. At
rank17, keep A17 as a separate coupon-bearing type and D17 as the residual
type. B17 has no such coupons and is dominated by D17 under the explicit
base-cost hypothesis. Do not identify A17 with D17 merely because they
have the same rank.

The finite forest is partitioned into small residual components of ranks
1 through9, residual D components of ranks17 through75, and the large
coupon-bearing types B10,D11,...,D16,A17. Enumerate every feasible vector
of counts of the latter types:1,608 vectors through rank75. Their coupon
environments give675 distinct additive Bellman tables, each with76 states
including rank0. Every state maximizes over every residual component rank
that fits; the saved last-rank predecessor supplies an attaining forest.
Combining every large vector with every residual rank proves the upper
bound and attainment at each total rank.

Writing a for total rank in components of rank at most5 and t for the
number of D11 through D16 components, the doubled cross-cost is

    a*(3*t + 2*b_A17)
    + 2*a1*(b_D11+b_D12) + 2*a2*b_B10 + 2*a4*b_D13.

The full D11 and D12 upgrades replace their old halves rather than count
them twice. The two A11 cycle coupons are also included in the bound:
the frozen A11-to-D11 substitution still preserves rank/base and dominates
every interaction, since no rank11 cost or coupon was changed. Every
saved maximizing forest omits A11, so the maxima remain attained.

|Ambient rank|Mild maximum|Mild finite slack|Full maximum|Full finite slack|
|---:|---:|---:|---:|---:|
|66|502|−7|496|−1|
|67|513|−3/2|504|15/2|
|68|524|3|512|15|
|69|528|16|520|24|

The full-rule maxima in these four ranks are attained by the single
finite components D66,D67,D68,D69. Their unchanged cost is8n−32. Hence
removing the rank17 obstruction exposes a different finite obstruction
at ambient66: D66 costs496 against curvature495.

The mild rule would suffice only to exclude68, giving ceiling67 once
coverage were proved. It fails the finite comparison at67 by3/2. The
stronger rule fits both67 and68; these statements must not be conflated.

## The unchanged ideal comparison

The source freshly enumerates the completed68 theorem's concentration
domain and imports its exact polynomial arithmetic. Arbitrarily repeated
ranks are10,11,13,14,15,16,17; at most one part lies in1 through9, at most
one is12 and at most one is at least18.

At ambient67 all1,761 terminal profiles have positive margin. The minimum
is851/165888, attained by six rank11 factors. At68 all1,906 are positive,
with minimum411499/3234816. At69 all2,060 remain positive. The exact
finite slacks above and the positive global constants330/335 would
therefore exclude67/68 under the full rule. Ceiling68 supplies n>=69.

Ambient66 still has a negative ideal margin−1565/10976 at five rank13
factors, in addition to its finite deficit1. No ceiling65 follows from
these costs. The ideal comparisons are inherited-formula calculations,
not a new independent implementation of that arithmetic.

## Verification scope

Generation22631 and fresh separate-process replay38315 terminated
exit0/PASS after explicitly distinguishing the mild and full rules.
A [separate standard-library certificate verifier](fresh_rank17_reduced_forest_frontier_20261017_review.md)
independently reconstructs both1,608-vector domains and675 coupon tables
per variant. It checks every Bellman upper alternative and attaining
predecessor, every global upper bound and every typed attaining witness.
Its report records final bound generation and replay. A full mathematical
source/report read passed the rank17 A/D separation, coupon counts,
unchanged A11 domination and the distinct mild/full conclusions. The
ideal profiles are not claimed as independently reconstructed by that
finite-budget verifier. Source, data and report are final/frozen.
The finite-type classification,
unchanged type domination, A11 coupon replacement and convex-order lemma
remain frozen proof inputs. No external refereeing or formalization is
claimed, and this abstract budget alone is not an ambient dimension proof.
