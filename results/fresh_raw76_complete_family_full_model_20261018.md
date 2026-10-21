# Complete raw76 family in the exported global model

This packet adds all242 frozen raw76 family inequalities to the finite
compact arithmetic H13 model over Q(sqrt(5)). Every one is new to the
current327 matrix. The one authorized numerical optimization reached its
300-second limit without returning a primal or dual proposal. It was not
restarted. The changed matrix nevertheless has an exact fallback enclosure

    102699323852930280/8367588524424317 <= optimum
        <= 6174801959947/500000000000,

approximately12.2734672664 to12.3496039199. The lower endpoint exceeds12,
so this exported matrix still has a certified obstruction to a target12
primal. The optimum within that interval is unresolved. No new optimized
face masses, chamber, dimension bound, or complete enlarged packing cones
are claimed.

The lower certificate is the fixed-face-mass repair already proved in the
[complete raw76 family packet](fresh_raw76_extra_vertex_complete_family_20261018.md),
now checked against every column of the enlarged global matrix. The upper
certificate is the [current327 primal](fresh_current327_full_model_20261018.md),
extended by assigning zero to all new variables. Thus these endpoints are
exact feasible certificates, not outputs recovered from the timed-out solve.
In particular, the entire242-column extension is insufficient to eliminate
the target12 obstruction in this exported model, irrespective of the
unresolved location of its optimum inside the certified interval.

## All family charges and all inherited columns

The [producer](../scripts/fresh_raw76_complete_family_full_model_20261018.py)
reconstructs the complete inherited matrix and all242 saved family charges.
The geometric input covers235 raw76 classes,235 target forests and3204
markers. It establishes one extra actual vertex per source quad, complete
finite-frame recovery, all13 relative marker patterns, and marker capacity
min(a−1,h) in a target with a A1 components and h H3 components. This packet
inherits that proof and does not repeat its physical pair tests.

Each charge combines a nonnegative ordinary-orbit charge with marker weight
1/capacity. The producer rechecks1048 inequalities against the complete
saved ordinary count-vector domains and reconstructs face coefficients
from actual corner origins. At target1304, the one ordinary raw178 orbit
compatible with markers is assigned zero in the new charge. Its six
unknown marker interactions remain allowed. All750 unknown marker-pair
edges also remain allowed. No historical Q–Q outcome is altered.

Joint integer normalization includes both face coefficients and the vertex
cost. All242 charges are distinct and none matches an inherited column;
10 strictly violate the old exact current327 dual. Their validity holds at
arbitrary face masses. Old attaining mixtures are not used as optimal
covers at different masses.

|Quantity|Count|
|---|---:|
|Inherited columns, preserved literally|439856|
|New family columns|242|
|Total columns|440098|
|Active columns|174272|
|Face requirements|82806|
|Vertex budgets|6558|
|Zero-budget signed columns|96176|
|Positive-cost signed columns|52|
|Inherited two-owner columns|14|
|Face-matrix nonzeros|830177|

The exact matrix digest is
`5487730df8c3a7358b636f741ac5db3992d20be09941502f13cc9a89c89053b4`.
Every inherited variable record, face requirement, signed coefficient,
full budget coefficient, zero-budget record and provenance mapping remains
literal. Both costs of each two-owner column are retained.

## The single proposal and exact fallback

The prepared numerical matrix has89364 rows,174273 variables including
the objective allowance, and914845 nonzero entries. Preparation took
71.26 seconds. The sole solve ran in session46090, printed its limit result
in chunkb69cdf and exited0 in chunk843d6e. HiGHS returned status1,
"Time limit reached", after300.084690125 seconds, with objectiveNone and
no solution or marginal vector. The saved reported iteration count is121207;
no convergence claim is inferred from it. The invocation counter remains1.

The [fallback certifier](../scripts/certify_raw76_complete_family_fallback_20261018.py)
makes no optimizer call. It retains the original failed proposal JSON
byte-for-byte, checks the unchanged matrix digest and source bindings, and
checks two certificates against the complete exact matrix:

* The inherited primal has all new variables zero. Every old face
  requirement and vertex budget is unchanged, so the old primal remains
  feasible. The full verifier independently checks those obligations.
* The dual keeps all4624 old positive face masses. Its vertex allowances
  are the frozen family repair. Fixed face masses leave every zero-budget
  signed load unchanged. Every nonzero budget coefficient is positive, so
  raising vertex allowances preserves all other old inequalities, including
  signed face loads and the14 two-owner costs. The certifier additionally
  checks every old and new column directly.

The failed proposal and its solver metadata remain recorded in
[the proposal snapshot](fresh_raw76_complete_family_full_model_20261018_proposal.json).
The [certificate JSON](fresh_raw76_complete_family_full_model_20261018.json)
explicitly sets optimization_resolved=false and changed_face_masses=false.
The exact interval is broader than a completed optimizer recovery; it is
not described as a newly optimized value.
Its exact width is
318540094180158291466831199/4183794262212158500000000000.

## Saved support and the next structural target

The [saved-only support inventory](../scripts/fresh_raw76_complete_family_support_20261018.py)
classifies the surviving fallback face masses using the frozen provenance
map and canonical-core catalogue. It distinguishes ordinary-corner coverage
of327 classes from the extra-vertex marker constraints for235 raw76 classes.
Merged provenance contexts are retained as merged rather than assigned an
unjustified unique raw family. Rankings are exact weighted demand, not a
predicted improvement or a newly proved geometric obstruction.

There are still4624 positive face masses:1949 belong to quad classes,
including421 merged, unexpanded provenance classes. Of the ordinary327
classes,224 have positive mass. Sixteen marker-family faces have positive
mass, totaling57192559/5000000000. There are1288 unambiguous positive quad
classes in213 raw families outside both represented sets. These are
existing support gaps, not a new optimized migration.

The largest missing-family demand totals are:

|Raw family|Missing positive classes|Exact weighted demand|
|---|---:|---:|
|14|62|72193209/5000000000|
|112|17|122289/10000000|
|57|30|121249627/10000000000|
|77|28|14250007/1250000000|
|108|11|92685609/10000000000|
|56|45|720749/78125000|

The highest single missing demand is face48584, raw112, source context68194,
canonical core92: its mass is10403311/5000000000 and required coefficient3.
The largest family total, raw14, has six quad partitions within each
canonical aggregate; treating it as one ordinary quad would miscount it.

A smaller concrete structural target is the pair of six-wall raw families56
and57. The bound canonical catalogue records one quad partition
{1,3}|{0,2} in each core and five positive-definite two-wall deletions:
{0,1},{0,3},{1,2},{2,3},{2,4}. The last is outside the four quad vertices,
so a fifth-vertex recovery mechanism is a specific next question for these
families. Their missing positive demand totals are0.0092255872 and
0.0121249627. In contrast, raw77's saved catalogue has only the four quad
deletions, so it does not offer this same immediate candidate.

This comparison is a read of saved catalogue entries45,46 and66, not a new
native reconstruction or incidence theorem. It does not establish the
missing-normal recovery, marker multiplicities, compatibility or capacities
needed to export any new column. A next investigation should cover the
whole relevant spectator family rather than only the positive faces.

## Verification scope

The [standalone saved checker](../scripts/check_raw76_complete_family_model_saved_20261018.py)
uses only the Python standard library. It independently checks all242
charge reconstructions,1048 old-vector inequalities, primitive scaling,
source identities, old-dual violations, new dual loads and all14 two-owner
loads. It checks saved primal nonnegativity and new-column budget
contributions, exact interval arithmetic, the failed-proposal snapshot,
and literal fallback certificate provenance. It does not reconstruct the
old matrix or native geometry. Full old-column matching, all primal face
requirements, the face-weighted numerator and all dual-column inequalities
belong to the separate full-matrix replay.

Fallback certification session55334, terminal chunk72a550, passed exit0.
Separate complete-matrix replay session14287, terminal chunkae5cdb, passed
exit0. The standard-library checker session92748 printed PASS in chunkc4628c
and exited0 in chunk27a5fe. Saved support generation session95582,
terminal chunk47d99c, and its separate replay session77045, terminal
chunk65451e, both passed exit0. The producer's failed optimization remains
recorded separately; these later passes do not change its status.

The root completed a full producer, fallback certifier, standard-library
checker, saved-support producer and report mathematical accounting read.
That review checked literal439856-column preservation and the242 new
charges, arbitrary-mass validity and joint normalization, all14 owner
costs, exact fallback verification, the zero-extension primal, unchanged
failed-proposal snapshot and one-solve scope, and the distinction between
merged provenance, ordinary327 coverage and marker235 coverage. It did not
execute computation or independently review the native geometry.

The independent peer completed a full producer, fallback certifier, saved
checker and report mathematical accounting read, including the inherited
canonical, cost_rows and verify_dual interfaces. It checked exact charge
scaling, literal matrix preservation, fallback provenance and the helper's
limited full-matrix scope. It performed no numerical, matrix, geometry or
clique replay. Its final full support-inventory source and report-delta
read also passed: unchanged masses, exact weighted-demand accounting,
explicit/unique/merged provenance, distinct ordinary and marker coverage,
and the catalogue-only next candidate are correctly separated. This final
delta review did not execute the saved arithmetic.

The final integrity audit checks26 direct binding entries covering24
distinct files, all four Python ASTs and eight local report links. The
packet is final, with an unresolved numerical optimum and a completed
exact fallback certificate. No frozen predecessor is changed.

    python3 scripts/fresh_raw76_complete_family_full_model_20261018.py --verify
    python3 scripts/check_raw76_complete_family_model_saved_20261018.py
    python3 scripts/fresh_raw76_complete_family_support_20261018.py --verify
