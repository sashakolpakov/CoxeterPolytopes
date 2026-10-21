# All five complete marker export packets in one changed full model

This packet incorporates every frozen raw56/57/76/108/112 export in the
current compact arithmetic H13/Q(sqrt5) finite matrix and permits all face
masses to change. Exact primal and dual certificates give the interval

    4275709870626110/346500030107667
        ≤ finite-model optimum ≤1233971157349/100000000000,

or12.339709954130539... through12.33971157349. The exact gap is
56110810193068294783/34650003010766700000000000,
approximately1.61936e−6. In particular, this finite matrix has no feasible
target12 primal. The varying-mass dual restores an obstruction above12;
the earlier fixed-mass11.9775975 witness did not settle this question.
This is not an arithmetic H13 exclusion or a complete geometric model.

The [producer](../scripts/fresh_all_marker_families_full_model_20261018.py)
preserves the [current327 model](fresh_current327_full_model_20261018.md).
It does not repeat the timed-out raw76 solve or reuse its optimization
status. Its new matrix includes these five frozen export sets:

|Packet|Export origins|
|---|---:|
|[Complete raw76 family](fresh_raw76_extra_vertex_complete_family_20261018.md)|242|
|[Complete raw56/57 families](fresh_raw56_raw57_extra_vertex_family_20261018.md)|1212|
|[Joint raw56/57/76 family](fresh_raw56_raw57_raw76_joint_marker_family_20261018.md)|1430|
|[Complete raw108 family](fresh_raw108_complete_marker_family_20261018.md)|415|
|[Complete raw112 family](fresh_raw112_complete_marker_family_20261018.md)|529|

The3828 origins give2354 distinct charges after exact normalization of
face coefficients and vertex cost together. All2354 are absent from the
old matrix. Exactly112 distinct charges violate the old optimized exact
dual. This is a changed optimization problem, not a repeat solve with the
same constraints.

## What is inherited and what is reconstructed

Each frozen family proves its actual extra-vertex incidence, complete native
role transport, conservative pair outcomes and shared marker capacity.
An aggregate with three quadrilaterals contributes once at each proved
extra vertex; it is not charged three times at that vertex. These
geometric theorems and their complete relative-pattern proofs are inherited.
No physical pair tests, new graph census or enlarged-cone enumeration are
performed here.

The producer reconstructs the chronological600 ordinary count-vector
domains from the latest local records and the two subsequent raw76 patch
packets. For every export it checks the actual-origin coefficients of the
ordinary orbits, nonnegative original weights, and the exact saved
count-vector hash. There are7622 ordinary-vector charge inequalities across
all export origins. Compatible ordinary coordinates are zeroed by the
recorded omission union, including all older marker families. No old
ordinary Q–Q outcome is changed, and no retained unknown becomes an
exclusion.

The marker coefficients satisfy5447 exact saved boundary inequalities.
The completeness and geometric validity of those boundaries are inherited
from the frozen family theorems. Nonnegative coefficients and their
ordinary-omission union give one shared vertex allowance. The exported
face coefficients are recomputed from retained ordinary weights and the
once-per-extra-vertex marker terms. These inequalities are valid at
arbitrary masses. Earlier attaining covers or fixed-mass optima are not
reused as covers at new masses.

Positive rational coefficients are normalized jointly with the complete
vertex cost: clear denominators, then divide every face and cost entry by
their common integer gcd. Duplicate origins remain recorded. Matching
uses every face coefficient and every budget coefficient of every old
positive-face column; it cannot match a zero-budget transfer. No such
old match occurs in this packet.

## Complete finite matrix preservation

|Quantity|Count|
|---|---:|
|Inherited variables, preserved literally|439856|
|Appended family variables|2354|
|Total variables|442210|
|Face rows|82806|
|Vertex budgets|6558|
|Active variables|176384|
|Zero-budget signed variables|96176|
|Inherited positive-budget signed variables|52|
|Preserved multi-owner columns|14|
|Face-matrix nonzeros|833679|

Every inherited variable record, signed face coefficient, face requirement,
budget coefficient, zero-transfer record and provenance map is checked
against the reconstructed current327 input. All14 multi-owner columns
retain every owner cost. Their metadata is not substituted for the actual
budget matrix. The inherited positive repair variables remain unchanged.

Before optimization, the entire exact new matrix is checked against the
frozen raw112 repair at its old face masses. Its ratio is
14671331978990040/1224897727063181=11.977597520869008...
This only supplies a feasible starting dual and an accounting check; it
does not bound the varying-mass optimum from above.

## Numerical proposal and exact certification

The installed environment is Python3.14.6, NumPy2.4.4, SciPy1.17.1 and
direct HiGHS1.14.0. The numerical matrix has89364 rows,176385 columns
including the objective allowance, and920459 nonzeros. Preparation took
77.08 seconds in the authoritative attempt.

The first numerical method was direct HiGHS dual simplex, strategy1,
four threads, fixed seed0 and a900-second solver limit. Its initial primal
is the inherited exact primal, extended by zero on every new variable.
That extension is first verified against the entire exact matrix. Floating
conversion is only a numerical warm start. The direct interface saves
the returned primal and row-dual vectors before rational recovery, even
if a time limit leaves the numerical optimization unresolved. It reached
the time limit after301,714 feasibility-phase iterations. Its displayed
objective12.000055256735232 was not a feasible optimum. Exact repair of
that iterate produced bounds weaker than the inherited certificates;
the iterate and its exact recoveries remain saved without alteration.

The [continuation](../scripts/fresh_all_marker_families_dual_continuation_20261018.py)
uses the identical exact matrix and model digest, with a different
numerical formulation. For face requirements r, face matrix A and full
vertex-budget matrix B, it minimizes

    −r·alpha,
    subject to Aᵀalpha−Bᵀmu≤0, sum(mu)≤1, alpha≥0, mu≥0.

This is the normalized dual of the original primal minimization.
Columns omitted from its numerical inequality matrix have no positive
face entry; this is checked explicitly, so their inequalities hold
automatically for nonnegative alpha and mu. The remaining numerical
matrix is the transposed original active matrix:176385 rows and89364
columns, with the same920459 nonzeros. It uses HiGHS IPM with presolve,
four threads, seed0, a900-second limit and crossover disabled. No basis
is supplied, allowing presolve to run. The exact frozen feasible dual
is retained as a candidate regardless of numerical progress.

There are two chronological numerical invocations on one changed model,
with separate immutable returned-vector snapshots. No geometry or pair
census is repeated. Continuation column values feed the original exact
dual recovery; negatives of its row multipliers feed the original exact
primal recovery. The final interval selects the strongest verified dual
and the smallest verified primal maximum among both methods, the frozen
raw112 fixed-mass dual, and the inherited current327 primal.

The continuation finished95 IPM iterations in76.49 seconds. HiGHS returned
`Optimal` with numerical objective−12.339711099623605, after an IPX
imprecision warning and postsolve status change. Both facts are retained
in the [solver transcript](fresh_all_marker_families_dual_continuation_20261018_solver.txt).
The exact interval above, not numerical status or a zero-gap assertion,
is the certificate. The [final saved result](fresh_all_marker_families_dual_continuation_20261018.json)
contains both chronological solver statuses and selects the continuation's
verified primal and dual. The first method's
[snapshot](fresh_all_marker_families_full_model_20261018_numerical.json)
and [transcript](fresh_all_marker_families_full_model_20261018_solver.txt)
remain unchanged.

For primal recovery, the inherited exact engine rounds nonnegative
variables at denominators10^8,10^10 and10^12 and adds the least common
uniform increment needed on the inherited positive repair variables.
It checks every face inequality and recomputes every actual budget,
including all owners. The smallest recovered maximum is an exact feasible
upper bound for this matrix.

For dual recovery, face and budget multipliers are rounded at denominator
10^10. The inherited simplex transfer order repairs every zero-budget
signed inequality; those inequalities are then checked exactly. Remaining
deficits are repaired by increasing a vertex multiplier through the full
budget coefficient list. Every nonzero budget coefficient is positive,
so an increase cannot invalidate another column, even if that column has
signed face coefficients. The weighted numerator and normalization are
recomputed from the complete matrix. This gives an exact feasible lower
bound, independently of the solver's numerical termination status.

A separate continuation `--verify` invocation reconstructs the same model digest and
checks every recovered primal and dual obligation without invoking the
optimizer. Exact feasibility, an exact interval, and numerical optimality
are distinct claims. No exact optimum is asserted.

## Independent checker and support scope

The [standalone checker](../scripts/check_all_marker_families_model_saved_20261018.py)
uses only the standard library and no producer imports. It independently
reconstructs all3828 ordinary/marker coefficient calculations, ordinary
and boundary inequalities, joint primitive normalization and complete
origin accounting. It checks the exact new-column and14 multi-owner
dual loads, nonnegative primal data and interval arithmetic.

It intentionally does not rebuild the old matrix, verify the old signed
constraints, recompute the weighted objective or check every primal face
obligation. Those belong to the separate full-matrix replay. It also
inherits the frozen geometric envelopes, historical owner chronology
and ordinary graph/count-vector completeness. Its scope is not a second
independent reconstruction of all historical mathematics.

The [saved-support inventory](../scripts/fresh_all_marker_families_model_support_20261018.py)
compares a recovered exact dual with current327. It keeps explicit,
unique-context and merged-unexpanded provenance separate, and distinguishes
the327 ordinary classes from the2035 source classes with marker constraints.
It verifies the weighted-demand numerator. Missing-family rankings are
diagnostics, not incidence theorems or guaranteed improvement amounts.

The [saved support](fresh_all_marker_families_model_support_20261018.json)
has5855 positive face multipliers, compared with4624 previously:
2471 become positive,1240 become zero, and3384 retained positive entries
change. There are2432 newly positive entries outside both ordinary327
and marker2035 coverage. Of the positive entries,596 have merged,
unexpanded provenance and are excluded from the unambiguous-family
ranking. There remain1529 unambiguous uncovered classes in214 raw families.

|Unambiguous missing family|Positive classes|Weighted demand in this exact dual|
|---|---:|---:|
|raw77|34|2730923/250000000|
|raw357|9|83284683/10000000000|
|raw14|47|30057453/5000000000|
|raw33|18|2227107/400000000|
|raw135|17|52173831/10000000000|

These are saved multiplier-scale demands, not predicted bound reductions.
The new certificate assigns much less mass to the five marker source
families; only53 of their2035 classes remain positive. This observation
does not prove uniqueness of the optimized masses or a causal necessity.
Raw77 and raw14 already have a
[frozen within-core no-extra-vertex result](fresh_raw14_raw112_raw77_raw108_core_incidence_pilot_20261018.md),
so another incidence mechanism would be needed there. Raw357 is a fresh
catalogue candidate; no new geometric proof or capacity is claimed here.

Even a target12 primal would certify only this finite exported matrix.
The missing ordinary cone facets, missing geometric classes and other
coverage limits of the inherited model would still need their own proofs.
Neither such a primal nor a ratio above12 establishes a polytope,
an arithmetic H13 exclusion or a new dimension theorem by itself.

## Verification

Run27009 completed the first preparation and numerical method, with
terminal chunkf8275b and a solver time-limit status. The equivalent-dual
continuation49513 completed terminal0 at chunk474c79. It is a second
numerical method on the same model, not another changed-constraint
experiment.
An earlier preparation58198 stopped before any optimizer call because
older frozen packets use different saved graph-status field names. The
explicit schema guard was corrected; no frozen packet was changed and
no numerical invocation occurred in that preparation attempt.

The independent saved checker completed terminal0/PASS at05d4fb.
Support generation completed terminal0/PASS at8a5047 and its separate
replay ated4e4a. The separate full-matrix replay60563 completed
terminal0/PASS at9862fa. All numerical and verification processes are
finished.

The root reviewer completed full source/accounting reads of the original
producer, standalone helper and report, then the continuation and final
report. This included the normalized-dual signs, inactive-column argument,
exact recovery interfaces, complete-matrix verification path, preserved
owners and honest numerical-status distinction. The root's final support
review covered the report, not a fresh support-producer read. No execution
or independent matrix reconstruction is attributed to that review.

The order16 reviewer completed full mathematical/accounting reads of the
producer, continuation, helper, support producer and final report. It also
read the inherited canonical normalization, full-budget dual recovery and
primal recovery/verification implementations. Its final delta review
covered both immutable histories, the IPX warning, exact interval and
support-provenance scope. It did not execute optimization, full-matrix,
support arithmetic, native geometry or historical census checks.

Final integrity checks rehash all87 direct source bindings, covering44
distinct bound files, parse the four scripts and check all16 local report
links. The packet is final and frozen; the exact feasible interval is
the result, with no exact-optimum or new geometric-bound claim.

    python3 scripts/fresh_all_marker_families_dual_continuation_20261018.py --verify
    python3 scripts/check_all_marker_families_model_saved_20261018.py
    python3 scripts/fresh_all_marker_families_model_support_20261018.py --verify
