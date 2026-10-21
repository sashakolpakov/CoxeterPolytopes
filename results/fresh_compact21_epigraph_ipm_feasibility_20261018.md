# Complete compact21 epigraph: reduced-column IPM feasibility attempt

18 October 2026. This packet records a completed bounded numerical attempt,
not a dimension improvement or an exact infeasibility theorem.

## Exact formulation and safe reduction

The source uses the independently reviewed complete Bellman formulation in
`compact21_full_forest_bellman_formulation_20261018.md`. It retains all 1979
native faces and 102917 spectator cases, deduplicating only exactly identical
face inequalities (51443 distinct rows). All 68 finite component types and
729 exact hub contexts are represented through the Bellman epigraph. The
proved D5 active/inactive rules, H4 funded transfer, and strengthened A5
capacity one are retained. Literal reconstruction of all 510 saved forest
support rows passes before solving.

Every nonnegative bonus variable that has no negative coefficient in any
native face inequality may be set to zero: it changes no face payment and
can only decrease the left side of every other inequality. The worker checks
all face bonus coefficients are nonpositive and all 167865 nonface bonus
coefficient occurrences are nonnegative. The D5 inactive extreme-E8
coefficient is zero after its two-unit cancellation. This elimination keeps
z, t, and all Bellman potentials. Thus it preserves feasibility and the
optimal value; it is not a heuristic truncation.

The retained matrix has 79274 rows, 32260 columns, and 378326 nonzero
coefficients. Of 93045 original bonus variables, 30599 are retained and
62446 eliminated. The full original global column space has 94706 columns,
including z, t, and 1659 auxiliary rank states. The source constructs an
explicit compressed-to-global column map; no primal certificate was
obtained, so no large redundant map is emitted in this failure packet.

## Measured outcome

The feasibility objective is zero and z is fixed to 3999999/100000
(39.99999). SciPy linprog calls bundled HiGHS with method highs-ipm,
one computation thread, presolve enabled, and crossover disabled. A separate
tiny local test had checked that the installed wrapper accepts the string
option run_crossover=off.

The completed attempt (tool session 47813; final chunk c4ad66) built the
matrix in 4.031 seconds. HiGHS presolve returned status 2, infeasible,
after zero IPM and zero crossover iterations, without a primal vector.
The parent measured 4.283 seconds total and a combined process-tree peak
RSS of 266960896 bytes. Initial resources were 196.90 GiB free disk and
6.49 GiB available memory. The live cap was 1 GiB RSS and 59 seconds;
the computed numerical solve cap was 38.9694 seconds. Both temporary
certificate paths were removed and the owned process group was stopped.

This is numerical evidence that this complete conservative payment model
cannot meet the chosen fixed z. There is no rational Farkas certificate,
no exact lower bound, and no conclusion that alternative payment rules
cannot establish dimension at most 20. No new dimension ceiling is claimed.

An earlier attempt (session 78684) stopped after matrix construction
because its driver omitted the OUTPUT_LIMIT constant. The original source
and compact measured failure ledger are preserved separately under
`*_driver_bug.py` and `*_driver_bug_failure.json`. The corrected source
adds the intended 512 KiB limit. No scientific row or solver option changed.

## Reproduction and next obligation

The frozen corrected source can reproduce the capped attempt with
`python3 -B scripts/fresh_compact21_epigraph_ipm_feasibility_20261018.py`.
A repeated attempt must preserve this evidence and use a new result stem
if any scientific strategy changes. No retry was launched after the
numerical infeasibility result. A useful next step is an exact compact
obstruction or an independently checked feasible vector at a higher z,
so that the responsible geometric support can be identified. Neither is
asserted by this packet.
