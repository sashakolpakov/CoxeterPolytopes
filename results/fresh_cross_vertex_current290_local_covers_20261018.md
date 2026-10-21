# Fresh local covers at the new arithmetic dual

In the compact arithmetic H13 model over Q(sqrt5), recomputing every
current actual-corner cone at the new5121 face masses
produces16 new violated charge columns. Exact monotone repair gives

    205569148850452800 / 16640368322881609
      =12.353641750091649098... .

This is still above12. It is a fixed-mass certificate for the enlarged
charge family, not a new optimized model value or a dimension theorem.
All571 current forest cones have exact rational charge/cover certificates
at the new demands. No old5053-mass optimum is reused.

The [producer](../scripts/fresh_cross_vertex_current290_local_covers_20261018.py)
and [certificate](fresh_cross_vertex_current290_local_covers_20261018.json)
bind the [new full model](fresh_cross_vertex_current_full_model_20261018.md)
and the complete current native/count-vector ancestry. The
[support migration](fresh_cross_vertex_support_migration_20261018.md)
states which classes these cones cover and which new positive classes
remain outside them.

## Fixed geometric domain, newly computed marginals

The current290 actual Q classes give571 finite-vertex forests. For each
forest, the complete physical corner-orbit coordinates and the latest
complete maximal count-vector list are inherited from their frozen
necessary compatibility graph. The latest explicit pair outcomes,
including owned unknowns, remain exactly as certified. No pair is deleted
or added in this computation, and no partial clique enumeration is used.

For each orbit j, recompute its new marginal demand

    d_j =sum_f m_(fj)*alpha_f,

using the new full model's face masses and the complete actual occurrence
multiplicities m_(fj). These are actual whole corner-orbit coordinates,
not arbitrary single literals. The source assembles every latest saved
domain, checks all count-vector dimensions, and rechecks each saved
literal representative's orbit counts. Physical compatibility and
completeness of the vector list remain the bound geometric inputs;
the historical wall tests and clique searches are not rerun.

The exact new domain has:

|Quantity|Count|
|---|---:|
|Complete forest cones|571|
|Nonzero new demand vectors|512|
|Zero new demand vectors|59|
|Complete maximal count vectors|2200|
|Positive rational cover terms|1379|
|New violated charge columns|16|

The zero-demand cones have the exact zero certificate. The512 other
cones use freshly proposed local LP solutions, recovered and verified
over rational numbers. These are bounded local solves, not a global
model optimization. A progress file is keyed to all sources and the
new alpha; a saved certificate from another demand cannot enter it.

## Exact dual charges and covering mixtures

Let v range over the complete maximal orbit-count vectors of one forest.
The checker verifies nonnegative orbit charges w_j with

    sum_j v_j*w_j <=1  for every v.

Every actual simultaneous corner family lies in the necessary graph,
so its count vector is dominated by a maximal vector. Nonnegativity
therefore makes this a valid charge of total at most1 at that finite
vertex. Pairwise-compatible cliques may be geometrically unrealizable;
allowing them only weakens the charge constraint.

The same certificate contains nonnegative rational mixture coefficients
lambda_v satisfying

    sum_v lambda_v*v_j >=d_j  for every j.

The source checks every coordinate, every charge inequality and exact
objective equality

    sum_j d_j*w_j =sum_v lambda_v.

Weak duality then proves the exact local optimum for the complete saved
cone. A simultaneous whole-packing automorphism average gives the
uniform literal orbit marginals underlying this certificate, as in the
inherited local-capacity theorem. A mixture is a covering certificate,
not a claim that all its cliques occur geometrically.

This procedure reevaluates the full cone rather than just old optimal
charges. The saved-only support diagnostic found one old-charge
violation; the fresh exact cones find16. The largest required increases
are at forests788,2622,2481 and944:

|Forest|Required increase in its vertex allowance|
|---|---:|
|788|15801791/100000000000|
|2622|758199/5000000000|
|2481|210911/2000000000|
|944|875121/20000000000|

All16 columns are exported with their exact positive face coefficients,
one-unit vertex cost, orbit charge vector and load at the new alpha.
They are new reweightings of proved necessary graph cones, not new
physical incompatibilities.

## Preservation of the full model, including both owners

Keep every new full-model face mass fixed. Replace each local vertex
allowance by the maximum of its original allowance and the exact local
optimum. All other allowances are unchanged. The repaired normalizer is

    16640368322881609 / 16632000000000000,

and the unchanged numerator is679792158897/55000000000.

This preserves all439619 certified exported column inequalities by
monotonicity: their vertex-budget coefficients are nonnegative. It also
preserves the zero-budget signed inequalities because alpha is unchanged.
The source explicitly rechecks every one of209 added columns using its
entire budget row. The14 two-owner columns retain both positive vertex
costs, before and after the repair; neither owner is dropped or absorbed
into a fictitious single-owner capacity.

Consequently the repaired dual is valid simultaneously for the inherited
exported model and every nonnegative charge in the571 represented graph
cones. The rational covering mixtures provide the upper bounds needed
for this entire local charge family; the16 exports merely exhibit where
the previously optimized dual is cut. The repaired ratio is a lower
bound for any further optimization including these cones, not an upper
bound on its optimum.

## Verification and precise remaining scope

Generationf12278 and a separate fresh saved replayc98ffe terminated
exit0/PASS. Replay invokes neither local nor global optimization. It
reloads all complete vector domains, recomputes new demands, checks all
rational charge/cover inequalities and objective equalities, rechecks
literal orbit counts, and verifies the complete-owner cost preservation.
No old-alpha cover or local optimum is reused.

The exact-hub peer completed a full producer/report mathematical read,
including inherited `H.exact_profile` and `W.exact_lp/check_certificate`:
PASS. It checked positive-coordinate LP reconstruction against all full
vectors and demands, exact charge/cover equality, the59 zero domains,
and monotone preservation of every owner cost. The review performed no
arithmetic execution, native/owner/clique recensus or LP replay. Direct
metadata audit08f89b passed33 bindings,4 local links and the source AST.
This packet is final and frozen. These cones cover the current290 actual classes only:227 are
positive at the new alpha and63 are now inactive. The32 newly positive
historical classes outside290 have not been merged into these graphs.
Other missing Q classes, merged provenance rows, simplex features and
joint simplex/Q cones remain separate obligations. No claim of complete
physical realizability, target12 primal feasibility, a global optimized
value or a new arithmetic dimension bound is made.

    python3 scripts/fresh_cross_vertex_current290_local_covers_20261018.py
    python3 scripts/fresh_cross_vertex_current290_local_covers_20261018.py --verify
