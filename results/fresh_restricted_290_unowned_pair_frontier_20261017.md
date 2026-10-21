# Restricted witness: complete optimistic unowned-pair frontier

The remaining genuinely unowned pairs cannot remove the current named 248-Q restricted witness. Even if **all 570,398 such pairs** are declared incompatible, its exact repaired ratio is

\[
\frac{18670750401222348}{1549726013368729}
=12.047775051950415\ldots>12.
\]

This is a counterfactual ceiling for this fixed face-mass witness and its current conservative graph domain. It proves no new physical incompatibility, exports no actual capacity column, and leaves both actual witnesses unchanged. The actual full-mass ratio remains 11.999994075588122 and the actual restricted ratio remains 12.055550455551511. Neither value establishes target feasibility or a dimension theorem.

The executable source is [the restricted frontier script](../scripts/fresh_restricted_290_unowned_pair_frontier_20261017.py); its exact data are [the saved certificate](fresh_restricted_290_unowned_pair_frontier_20261017.json). The baseline is [the actual eight-pair 2599 successor](fresh_2599_raw110_raw24_pair_geometry_20261017.md).

## Complete current domain and ownership

All 571 forests of the 290-class graph are reconstructed from the frozen ancestry. The chronological 118-pair and eight-pair successors are applied after the raw125/raw205 domain. Each newly acquired pair is removed from the genuinely unowned set regardless of whether the test excluded it. In particular, the two surviving owned unknowns at forest1763 remain compatible and cannot be deleted in this diagnostic. All older full-domain owners and explicit tested pairs retain their outcomes.

The remaining domain contains 570,398 pairs in 854 orbit-pair domains. Only 1,828 pairs, in 31 domains, have positive restricted demand at both endpoints. The nine recent forest patches account for exactly126 acquired pairs, including the two owned unknowns. Current native stars, physical pair certificates and actual-corner multiplicities are inherited; this program does not recompute physical geometry.

## Why zero-demand projection is exact

For each forest, start from its current attaining restricted cover and charge. Set charge coordinates with zero restricted demand to zero. Nonnegativity means this cannot violate a charge inequality; the program checks every inequality on the **complete original maximal count-vector domain**, and the objective is unchanged.

Delete zero-demand stars from each literal clique in the attaining cover, retaining its coefficient. The resulting cliques remain compatible, cover every positive demand and have the same total mass. Matching charge and cover objectives establish the same optimum. Consequently deleting an edge with a zero-demand endpoint cannot increase the required restricted covering mass. These projected cover vectors need not be maximal: they supply a valid upper certificate, while the lower certificate was checked on the complete original domain.

This projection prevents irrelevant zero-demand stars from creating artificial splitting costs. An original legacy row omits a duplicate `demands` field; its demand is reconstructed from current face masses and actual occurrence coefficients and its exact certificate is still checked. No missing-field default is used as mathematical evidence.

## Strongest deletion and exact witness ceiling

Deleting more compatibility edges can only increase the covering optimum. Therefore the graph with every genuinely unowned edge deleted is contained in every graph obtainable by testing only a subset of those pairs. An explicit feasible cover for the strongest deletion remains feasible for each weaker deletion.

For all 571 forests the program supplies an exact baseline lower certificate and a valid hypothetical cover upper bound. Where needed, it enumerates the complete changed graph and obtains matching nonnegative charges and covers. Seven grand scenarios have complete exact cone certificates. In all forests, the quantity `max(old_mu, required_cover_mass)` is determined exactly, including cases in which the cover mass itself is only enclosed but remains below the old multiplier. There are no unresolved multiplier intervals.

Only six forests increase the restricted multiplier under the strongest deletion:

| Forest | Exact multiplier increase |
|---:|---:|
| 788 | 809431/2000000000 |
| 944 | 631907/4000000000 |
| 2365 | 29743/2500000000 |
| 2521 | 6061/1250000000 |
| 2622 | 488801/10000000000 |
| 2902 | 37741/2000000000 |

Their sum is 12943777/20000000000. The normalizer increases from 1548726494908789/1544400000000000 to 1549726013368729/1544400000000000. The other565 multipliers need no increase. Since the resulting ratio still exceeds12, **every weaker deletion of currently unowned pairs also leaves an above12 witness**. No physical pair-testing shortlist can overcome this obstruction within the stated domain, independently of the one-domain-per-forest shortlist search.

Fixed face masses and nondecreasing vertex multipliers preserve all439,410 inherited model columns, including zero-cost and positive-cost columns with signed face entries. This hypothetical witness additionally satisfies every nonnegative orbit-charge inequality of the strongest conservative graph cones. The result does not exclude new constraints from additional actual classes, higher-order incompatibilities, stronger tests deliberately revisiting owned unknowns, cross-incidence relations, or other missing geometric structure.

## Verification and review scope

Generation session90245 completed with exit0/PASS in60.39 seconds. Separate fresh saved verification session58895 completed with exit0/PASS in40.17 seconds. It reconstructs ownership and graph domains, checks literal covers and all saved rational certificates, and never invokes an optimizer. No global model solve was run.

A separate [standard-library `Fraction` saved-data checker](../scripts/check_restricted_290_unowned_saved_20261017.py), importing no producer or project arithmetic, passed from its final repository location: all14 direct bindings, 571 demand reconstructions, 590 charge/cover certificates, 1,838 charge inequalities, 1,615 mixture terms and593 scenario-bound records, and the exact normalizer and above12 ratio. This is an accounting check, not another graph or physical reconstruction.

The root agent read the complete source and final report, including chronological ownership, complete-domain validation of the zeroed charge, projected literal covers, exact `max(mu,t)` accounting, strongest-deletion containment and the inherited scenario interface. The construction agent independently read `current()` and `positive_cover()` and their surrounding calls and guards, with the same conclusion. These are source/mathematical reviews, not independent physical-pair, graph-enumeration or LP replays. The construction review does not claim a fresh audit of the inherited scenario/shortlist implementation.

All predecessor packages remain frozen and unchanged. The JSON binds the exact producer and helper sources and the three direct actual-geometry predecessors.
