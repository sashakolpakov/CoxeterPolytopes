# Changed global model after the raw76 extensions

This packet concerns the finite exported compact arithmetic H13 model over
Q(sqrt(5)). It inserts34 new proved local charge inequalities into the
frozen current324 global matrix. One changed-model global solve gives the
exact enclosure

    12.34959272328146 <= optimum <= 12.349603919894.

Both endpoints remain above12. The previous exact lower bound was
12.35590648805293, so the new constraints improve the finite-model objective
but do not remove its obstruction at12. No new dimension bound or chamber
construction is claimed.

The [latest local repair](fresh_new5073_raw76_third_class_local_covers_20261018.md)
has ratio11.952110486939478 at the fixed5073 face masses. That particular
least monotone witness repair is below12. The global optimization is a
separate question because its face masses and vertex allowances can change.
Neither a new above-12 dual nor a target12 primal would by itself settle
geometric feasibility outside this finite exported matrix.

## The genuinely changed constraints

The [previous optimized model](fresh_current324_full_model_20261018.md)
contained439822 columns. Since that solve, three completed local packages
produced34 exact violated exports:

|Source|Saved local records available|New exported charges|
|---|---:|---:|
|Current324 cones at the5073 masses|600|25|
|Raw76 classes42980 and43028|6|6|
|Raw76 class42946|3|3|

The [producer](../scripts/fresh_current327_full_model_20261018.py)
checks every one of those34 charge certificates against its complete saved
count-vector domain and its old attaining cover. This verifies530 rational
charge inequalities and244 mixture terms. The current324 domains are
transported chronologically from their frozen source records; the two raw76
packets provide their complete enlarged domains. No native pair, clique
census or local optimization is repeated.

Each nonnegative orbit charge w exports coefficient sum_i(w_i c_if) for
face f and vertex cost1, where c_if counts actual corner occurrences. Its
validity against every complete-domain vector is independent of the masses
at which its attaining cover was found. Those old covers are not reused
as covers at newly optimized masses.

All34 exports strictly violate the previous optimized exact dual. Joint
primitive normalization of both face and vertex-cost coefficients gives34
distinct charges; none matches any column in the entire inherited matrix.
Each is appended, retaining its source record, exact scaling, and old-dual
load and cost.

## Complete matrix preservation

|Quantity|Count|
|---|---:|
|Inherited columns|439822|
|Appended columns|34|
|Total columns|439856|
|Active columns|174030|
|Face requirements|82806|
|Vertex budgets|6558|
|Zero-budget signed columns|96176|
|Positive-cost signed columns|52|
|Preserved two-owner columns|14|
|Face-matrix nonzeros|829635|

The exact matrix digest is
`06c861427de31f38993789b9bc23e7c74c8a04ecbe0b13bc818c5ecde35ac8f3`.
Every inherited variable record, face requirement and coefficient, budget
coefficient, zero-budget record and provenance mapping is preserved
literally. The14 two-owner columns keep both vertex costs. The producer
also verifies the latest fixed5073 repair against the entire439856-column
matrix, including signed face loads and zero-budget constraints.

The resulting numerical proposal matrix has89364 rows,174031 variables
including the objective allowance, and914061 nonzero entries. Its size is
bounded by the previously completed model plus the34 new charges.

## One solve and exact recovery

The invocation guard allows one300-second HiGHS interior-point solve.
Primal recovery uses inherited positive old repair columns and verifies
all face requirements and vertex budgets exactly. Dual recovery applies
the inherited zero-budget simplex correction and raises only nonnegative
vertex allowances using the full budget matrix; every column is then
checked exactly. No fictitious single-owner cost replaces a two-owner cost.

The sole source-bound generation/solve session71699 completed with result
chunk1deb15 and terminal exit0 chunk91a75d. HiGHS terminated optimally after
58 iterations in66.758 seconds, at numerical objective12.349603919309306.
Exact recovery gives

    lower = 102699323852930280/8316008969212517,
    upper = 6174801959947/500000000000.

The displayed lower decimal is rounded outward below the exact fraction.
The exact gap is46555565141152765056599/4158004484606258500000000000,
approximately0.0000111966. The primal has23528 positive rounded variable
numerators and repair delta3/500000000000. The dual has4624 positive face
masses,1045 vertex allowances and1201 exact allowance-repair records.
Thirteen of the34 new columns have positive primal weight. No second
optimizer call is made.

Separate full-matrix replay session43639, terminal chunk249934, passed
exit0. It rebuilt the matrix and checked all primal and dual inequalities
exactly without an optimizer. The standard-library saved checker
session69312 printed PASS in chunk16cdd4 and exited0 in chunk7b524f.

## Changed mass support

The newly optimized certificate has4624 positive face masses, compared
with5073 previously:886 formerly zero masses become positive and1335 become
zero. Of the327 represented actual classes,224 have positive new mass.
Of the886 newly positive faces,864 lie outside those327 classes. These are
saved-support comparisons, not new geometric provenance or coverage claims.

All three raw76 classes receive zero mass in the new optimized certificate:

|Face|Previous5073 mass|New optimized mass|
|---|---:|---:|
|42980|76140277/10000000000|0|
|43028|53883347/10000000000|0|
|42946|31942181/10000000000|0|

This is a property of the saved certificate, not a uniqueness claim about
all optimizers. The complete local graph-cone covers at these4624 masses
have not been recomputed.

## Scope of saved checks

The [standard-library helper](../scripts/check_current327_model_saved_20261018.py)
imports no producer or solver. It independently verifies the34 saved local
charge/cover certificates, primitive scaling and all source identities,
old-dual violations, new dual-column loads and all14 inherited multi-owner
loads, saved primal nonnegativity and new-column budget contributions,
and exact interval arithmetic. It also reports the change in positive face
support. Historical native geometry and complete clique enumeration are
inherited from the bound local packets. Matching against the full old
matrix, old-column dual inequalities, all primal face obligations and the
face-weighted numerator remain obligations of the separate full-matrix
replay.

These selected exports do not exhaust the complete implicit packing cones
for the327 represented classes. If the optimized masses change, complete
local optima and attaining covers at those new masses require a separate
calculation. Missing actual classes, other finite forests and fields, and
global geometric consistency remain outside this finite matrix. The old
raw4 fixed-mass barriers are not asserted at the new optimized masses.

All computational checks are complete. The root completed a full
producer/report/helper mathematical accounting read: mass-independent
charge validity,34 certificate reconstructions and primitive matches,
strict old-dual violations, literal preservation of439822 columns and14
two-owner budgets, the full fixed-repair inequality loop, the invocation
guard, full-matrix verification and the helper's limited scope. That
review did not execute the solver or replay, or re-review the inherited
optimizer implementation. The construction peer independently completed a
full producer/report/helper mathematical accounting read, checking all34
charges,530 inequalities and244 terms, actual coefficient transport,
primitive normalization, inherited-column preservation and matching,
all14 two-owner budgets, the signed/full fixed-repair gate, changed-mass
scope and the helper's limited primal/old-matrix scope. That peer did not
execute arithmetic, optimization, native geometry or clique replay, and
did not newly read the inherited optimizer implementation. Final metadata
checks verify all19 direct bindings, seven local links and both Python ASTs.
The packet is final; frozen predecessor files are unchanged. No external
referee or proof-assistant validation is claimed.

Inputs: [new5073 separation](fresh_current324_new5073_local_covers_20261018.md),
[two raw76 classes](fresh_new5073_raw76_high_mass_local_covers_20261018.md),
[third raw76 class](fresh_new5073_raw76_third_class_local_covers_20261018.md).

    python3 scripts/fresh_current327_full_model_20261018.py --verify
    python3 scripts/check_current327_model_saved_20261018.py
