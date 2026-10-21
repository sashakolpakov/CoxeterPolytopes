# Exact budget 35 on the selected H3-center model

The selected payment model has an exact rational feasible rule with maximum
marked-vertex payment **35**. This is a restricted-model achievement, **not a
compact dimension-21 exclusion**. Newly enumerated triangle and paid-quadrilateral
rows fail this rule, as the [literal diagnostic](fresh_sqrt5_h3_center_selected_diagnostic_20261016.md)
records. The previously frozen compact arithmetic dimension bound of 21 is
unaffected: that certificate uses no discounts.

The [source](../scripts/fresh_sqrt5_marked_h3_center_sparse_recovery_20261016.py)
and [exact data](fresh_sqrt5_marked_h3_center_sparse_recovery_20261016.json)
retain all 32,054 payment variables of the
[H3-center family](../scripts/fresh_sqrt5_marked_h3_center_cross_20261016.py).
The distance-two discount is zero. The adjacent-3 discount is approximately
0.635928530811348; its full rational value is certificate.primal[1]. All bonus
coefficients are nonnegative.

The LP fixes the target at 35 and minimizes the sum of all bonus coefficients.
This objective is bounded below by zero and bounds every bonus on each objective
sublevel; both discount variables have finite bounds. The final exact recovery
uses 7,409 basis equations, 7,009 sparse rational pivots, and a 400-variable FLINT
remainder. It checks every LP inequality exactly. The complete finite-forest DP
then verifies maximum 35 across 14,545 context/mode fronts, and each chosen front
is independently recounted inside the producer by literal unordered root pairs.

The face domain consists of 196 selected triangle families, 541,151 marked
triangle contexts, and 1,979 inherited baseline-zero quadrilateral rows with
227,408 common spectator contexts. These yield 60,763 distinct face profiles.
The selected triangle set is deliberately overcomplete: one later-excluded
conjugate-inadmissible row remains imposed, which does not invalidate a feasible
rule. It is removed when replacing this set with the complete affected-triangle
catalogue in the subsequent repair model.

Generation and a fresh saved-certificate replay pass. This is a primary exact
certificate; it has no separate independent full LP review. The full newly
affected face domain has not been certified, and the rule already has explicit
new-face counterexamples.
