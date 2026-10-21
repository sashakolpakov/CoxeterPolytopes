# Raw357: a qualified inherited-charge barrier

There is an explicit fixed-mass feasible dual above12 for every extension
in the charge-template class defined below:

    4275709870626110 / 355157472905517 = 12.0389128677...

This is a barrier for that stated template, not for arbitrary new geometric
inequalities. It is neither an exported joint charge nor the least repair
of a new model. In particular, merely validating an old-coordinate charge
against a larger geometric domain does not establish the template's
domination hypothesis.

The [producer](../scripts/fresh_raw357_inherited_charge_barrier_20261018.py)
and [saved calculation](fresh_raw357_inherited_charge_barrier_20261018.json)
bind the frozen [full-model certificate](fresh_all_marker_families_full_model_20261018.md)
and the complete [raw357 incidence inventory](fresh_raw357_extra_vertex_incidence_pilot_20261018.md).
All47 source classes and all141 pointed extra-vertex domains are used.
Only nine classes have positive mass in this certificate.

## The precise extension class

Each additional row has one unit of budget at one target vertex. Its
nonnegative raw357 marker coefficients are at most one. After removing
those marker terms, its remaining face coefficient vector must be
coordinatewise dominated by a convex combination of already enforced
unit-budget columns at that same target. All old columns remain unchanged.

This last condition is stronger than geometric validity. A genuinely new
ordinary-star charge or a cross-vertex incidence inequality may fall
outside it. Every proposed future row must separately establish this
condition before this barrier applies to that row. The condition includes
downward scaling/omission of a single existing charge as a special case.

Write the frozen face masses as alpha and vertex allowances as mu. If
the old part of a new row is h, the hypothesis gives h.alpha <= mu_v.
For each target v, let e_v be the sum of alpha over its raw357 extra
occurrences. Each occurrence contributes once, regardless of how many
native orientations encode it. Since each marker coefficient lies in
[0,1], its total contribution is at most e_v. Therefore

    mu'_v = mu_v + e_v

satisfies every row in this template, simultaneously. Existing constraints
remain satisfied: alpha is unchanged, every nonzero vertex-budget
coefficient in the frozen matrix is positive, and zero-budget left sides
are unchanged. The producer also checks all14 saved multi-owner costs
directly with the extended allowances.

## Exact accounting and limits

The total source mass is83284683/10000000000. Each source quad has exactly
three extra actual vertices, giving total added demand249854049/10000000000.
The numerator remains61081569580373/4950000000000; the frozen normalizer is
115500010035889/115500000000000. Adding the full extra demand gives the
displayed ratio, strictly above12. Per-target overlap or smaller allowable
marker coefficients can only reduce the necessary increase and strengthen
this particular lower bound.

The numerical calculation reconstructs source/target multiplicities and
all current masses from the frozen artifacts. It inherits full old-matrix
feasibility and the weighted numerator from their completed exact replay;
it does not repeat that matrix replay. It invokes no optimizer and creates
no new geometric family or capacity claim. The next useful strengthening
must escape at least one of the explicit template hypotheses, or change
the face masses through genuinely new constraints elsewhere.

Generation66e092 and separate saved verification11e9bf terminated
exit0/PASS. Mathematical review is pending; this diagnostic is not yet
frozen.

    python3 scripts/fresh_raw357_inherited_charge_barrier_20261018.py --verify
