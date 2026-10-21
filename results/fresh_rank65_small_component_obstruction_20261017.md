# The finite rank65 obstruction also lies in small components

Reducing only payments involving whole finite components of rank at
least18 cannot improve the present exact rank65 maximum487. Even
setting every such payment to zero leaves that maximum exactly487,
against the required curvature allowance480.

The [source](../scripts/fresh_rank65_small_component_obstruction_20261017.py)
and [certificate](fresh_rank65_small_component_obstruction_20261017.json)
check the frozen exact budget and the two explicit smaller-component
forests below. This is an obstruction to a payment strategy, not an
assertion that either forest occurs at an actual hyperbolic vertex.

The completed finite rank17 rule pays127 pairs in D17. D16 still has
every pair at unit weight, costing120. No cross-component coupon joins
components of ranks16 and17. Therefore

    D17 + 3 D16: rank65, cost127+3*120=487.

This forest is unchanged by any modification confined to components
of rank at least18. Since deleting nonnegative weights cannot increase
the old exact maximum487, this unchanged attaining forest proves the
optimistic modified maximum is still exactly487. A new optimization
or rank cutoff extrapolation is unnecessary.

A second obstruction is

    A1 + 4 D16: rank65, cost4*120+4*(3/2)=486.

Each D16 has three eligible fork-region marks for the unchanged A1/D16
half-coupon, so the four component interactions add6. This forest exceeds
the finite allowance by6. The first exceeds it by7.

Thus lowering the large D65 cost alone cannot close the finite comparison.
Any successful local-cost approach must lower D16 or D17 same-component
cost: changing small-component coupons alone leaves the first forest
untouched. Coupon changes can address the second forest. Alternatively one needs a proved
restriction on actual vertex diagrams or a justified global compensation
argument. The new finite7 triangle theorem does not supply those missing
quadrilateral and budget obligations.

Final generation445279 and separate saved replayb39174 terminate exit0/PASS. A separate full source/report
read passed the exact monotonicity argument and its scope. An independent
saved-data arithmetic check counted the127 retained D17 pairs,120 D16
pairs and the three D16 fork-region coupon roles, giving both costs and
deficits exactly. That check did not rerun the finite optimization or
geometry. All seven direct bindings, two report links and Python syntax
pass. Source, data and report are final/frozen. No new geometric exclusion,
dimension theorem or finite-volume example is asserted.
