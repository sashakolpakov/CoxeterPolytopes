# A raw23-only witness survives the complete current coface matrix

Keeping every positive simplex and edge multiplier, and only the50
unambiguous raw23 Q multipliers, gives the exact feasible dual ratio

`4201253466877460/346500030107667 = 12.124828576701757… >12`.

The [producer](../scripts/fresh_raw23_surviving_coface_witness_20261018.py)
and [certificate](fresh_raw23_surviving_coface_witness_20261018.json)
reconstruct the complete frozen
[coface model](fresh_raw357_coface_full_model_20261018.md) and check all
442,226 column inequalities. This is the same finite matrix, with free
face masses, all signed variables, and all14 multiple-owner costs.
No optimizer is invoked and no constraint or historical compatibility
outcome is removed.

This is a weaker feasible witness than the frozen12.3313814866 witness,
deliberately chosen to isolate the remaining obstruction. Weak duality
shows that this finite matrix cannot attain target12. Neither ratio is
a new optimum or a geometric impossibility theorem. The result supplies
no new hyperbolic dimension bound.

## Exact deletion and matrix verification

The inherited provenance inventory distinguishes explicit, unique saved
context, and merged Q rows. Exactly50 positive Q rows have unambiguous
raw23 provenance, canonical core19 and requirement10. They are retained.
All other2,331 positive Q rows are set to zero, including every merged
row. The1,561 simplex and1,904 edge multipliers and every vertex
multiplier are unchanged.

The retained weighted objective is

|Contribution|Exact value|
|---|---:|
|Simplex rows|28636392980839/2475000000000|
|Raw23 Q rows|6932123/12500000|
|Edge rows|0|
|Total|30008953334839/2475000000000|

The normalizer remains115500010035889/115500000000000.
The deleted weighted demand is1032764639/5000000000.

Deleting positive face multipliers is **not assumed** to preserve a
signed system. The producer reconstructs the original matrix, checks
the inherited primal and dual, computes every candidate column load
and its full vertex cost, and finds zero violations. It separately
records that none of the deleted Q rows has a negative coefficient.
The original442,226-column prepared object matches the frozen one
exactly, including96,192 zero-budget signed records,52 inherited
positive-budget signed variables and14 multiple-owner columns.
The prior exact primal remains a feasible upper certificate and is
verified again; no new primal or optimum is claimed.

## What this rules out as a next tactic

Constraints involving only positive coefficients on the deleted Q
families cannot defeat this particular witness if their other loads
remain nonpositive. More concretely, the separate
[seven-profile coface pilot](fresh_core1_mixed_coface_capacity_pilot_20261018.md)
produces62 represented mixed capacities, including three violated by
the old witness. Every one is satisfied by this raw23-only witness.
All50 surviving raw23 masses are also at most each of their two actual
simplex-base masses, even with capacity1.

The qualified conclusion concerns these one-simplex coface inequalities,
or further nonnegative combinations dominated by the two checked
capacity1 inequalities. It is not a barrier to arbitrary new constraints.
A restriction on a raw23 singleton, a relation involving multiple bases
or vertices, or a new recoverable actual incidence could still exclude
the witness. Merely noting that raw23 prisms form a partial matching
between their two simplex-base types gives only the two existing
one-base bounds and their nonnegative combinations.

## Verification scope

Generation55964/3e77ea and separate full-matrix replay4387/d3b1ed
both terminated exit0/PASS. The saved checker5abc17 also passed.
An earlier computation reached the feasibility checks but failed while
writing metadata because the support inventory has no standalone MD.
That path was corrected to its combined report, checked before building,
and a short stale-path retry was interrupted. Neither failure launched
an optimizer or changed any frozen artifact.

The [standalone saved checker](../scripts/check_raw23_coface_witness_saved_20261018.py)
checks the exact deletion, retained masses, objective, normalizer and
the companion coface arithmetic. It does not reconstruct the old
matrix, prove historical geometry, or independently verify the old
primal. Those matrix obligations belong to generation and full replay.
Root's full mathematical producer/report/checker read passed, including
the mixed-cap producer and inherited prism/base and cross-coface reports.
It checked actual S/q/b recovery, intrinsic roles, provenance and
incidence multiplicity, the target-row distinction, exact lower and
qualified one-base scope. That review repeated no computation, full
matrix reconstruction or historical census. Construction's independent
full producer/report/checker review also passed, with both inherited
incidence/capacity reports and generic full-budget load/verification
interfaces reread. It repeated no execution or historical census.
Its requested saved-domain guards explicitly check raw list lengths,
five-coordinate contacts, and two distinct singleton-base checks per
retained face; the hardened checker957448 passed before freezing.

    python3 scripts/fresh_raw23_surviving_coface_witness_20261018.py --verify
    python3 scripts/check_raw23_coface_witness_saved_20261018.py
