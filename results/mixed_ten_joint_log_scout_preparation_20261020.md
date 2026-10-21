# Prepared joint-log relaxation: execution contract

20 October 2026. **Source prepared; not executed.** Only Python AST parsing
was performed. No imported model function, rational coefficient generator,
scipy call, optimizer or proof worker was run.

Worker: `scripts/scout_mixed_ten_joint_log_rows_20261020.py`, SHA256
`fa6731c982be9f4778fb7612c6d4d9d052309af3d4a0f18f6353f6137201d26b`.
Mathematical contract: `results/right_angled_ten_joint_log_transport_20261020.md`,
Its source and the ideal-four theorem are both streamed into the runtime
input manifest; the former was extended with the explicit new cover
dependency after the first prepared version.
The two original weaker full-log constants have been strengthened by
`log(1+x)>=2x/(2+x)` before this final handoff.

## Domain and scope

Finite vertex count is normalized to one. The model keeps the eleven
ideal-degree counts f0,...,f10, their irregular finite-edge incidences,
noncentral finite-direction triple counts, and unavailable R7 four-moves.
It includes both actual logarithmic graph rows, all exact degree-defect
and missing-P3 identities, the refined ordinary union budgets, and the
P3/Q8/R7/big-three-face partition. Big-face finite/ideal incidences and
their A0/A1/A2/A3 capacities remain separate.

It also includes every linear Dehn--Sommerville equation of the truncated
simple ten-polytope, nonnegative h and g coefficients, the inherited facet
minima and facet-plus-cusp minima recorded in the pinned older incidence
packet. After the ideal-four theorem and the independently checked parity
seed, it also uses the finite-vertex minima
`(v5,...,v10)=(6,27,178,1825,44206,8500814)` via the six actual face-cover
rows. Their complete paper proof is separately source-bound in
`results/right_angled_finite_vertex_seed_and_face_cover_20261020.md`.
The variable u represents the positive number 1/f; allowing u=0
is an explicit relaxation. R is unrestricted in sign. No dimension-eleven
strict weighted bound, unproved n10 density, actual realization of a
degree distribution, or common null vector between the two graphs is
assumed. Strict inequalities are relaxed. The source has 67 variables,
24 equality rows and 95 weak inequality rows by direct source counting;
the actual `--plan` result is still pending.

Selected-subset logarithmic coefficients are upper bounded by exact
Fraction arithmetic, using a rational enclosing square-root interval and
48 terms of the atanh series with an absolute tail bound. Missing-move
penalties use upper log bounds. Thus coefficient rounding weakens the
necessary rows. The LP conversion itself is floating point and is **not**
an exact certificate of feasibility or nonexistence.

## Root-only execution

The solver-free inspection argv is

    ["python3", "scripts/scout_mixed_ten_joint_log_rows_20261020.py", "--plan"]

The single numerical scout argv is

    ["python3", "scripts/scout_mixed_ten_joint_log_rows_20261020.py",
     "--write", "results/mixed_ten_joint_log_numerical_scout_20261020.json"]

The worker has **owns_lock=false** and must be launched through the
coordinator's shared-lock/resource wrapper. The usual 40GiB disk/6GiB
available-memory start gates and 8/3 abort gates apply, with one job total,
one solver thread, a one-GiB process-tree cap and a 59-second parent cap.
The worker sets a 45-second alarm; the one HiGHS solve has a 25-second
limit. The small model is expected to use well below one GiB. There is no
loop of solvers, background child, census or historical replay.

The model is built before scipy is imported, and scipy is imported only
after explicit `--write`. An existing output is refused. All source hashes
are streamed before and after execution, the output is bounded to one MiB,
and the ledger is exclusively created. Its status is always
`NUMERICAL_RELAXATION_ONLY`. A future infeasibility result would require a
separate exact rational dual and mathematical composition review before
any dimension claim. A feasible output would only be a numerical witness
for this aggregate relaxation, never a hyperbolic polytope or example.
