# Exact obstruction to the complete compact21 conservative payment model

18 October2026. This is an exact rational weak-duality certificate. It
proves that the stated complete charge model cannot attain z<=40, the
sufficient charge threshold for excluding dimension21. It does not prove
that compact arithmetic Coxeter21-polytopes exist, or that dimension20
cannot be proved by stronger geometric or incidence-dependent rules.

## Certified interval and scope

For the complete conservative model described in
`compact21_full_forest_bellman_formulation_20261018.md`,

  516372005670581/12600001078560 <= optimum
     <=16392763329/400000000.

Numerically this is

  40.981901703900085 <= optimum <=40.9819083225.

The upper is the previously frozen exact complete-model primal. The lower
comes from8034 negative rational multipliers on the79274 full inequalities.
All94706 original global columns are checked, including bonuses eliminated
from the numerical solve and every Bellman potential. The sufficient
threshold is z<=40, not only z<40; the exact lower is strictly greater
than40 and therefore rules out both for this model.

The model includes51443 exact native face classes from1979 faces and
102917 spectator cases, all68 finite component types,729 exact hub tuples,
the reviewed D5 active/inactive rules, the H4 funded transfer, and A5 source
capacity one. Its forest constraints are conservative: every abstract
allowed finite forest is imposed, without asserting that every forest
is geometrically realizable at an arithmetic polytope vertex. The exact
obstruction concerns these inequalities only. There is no new dimension
ceiling, existence theorem, or Lean verification here.

## Exact verification independent of numerical correctness

Write the primal as minimize x_z subject to A*x<=b and x>=0. The artifact
supplies y<=0. The verifier sums every stored rational row and checks

  (A^T*y)_j <=1 if j=z, and <=0 for every other global column j,
  b^T*y =516372005670581/12600001078560.

For any feasible x, y^T*A*x>=y^T*b while y^T*A*x<=x_z.
Thus x_z>=b^T*y. These literal rational comparisons are the proof; they
do not assume any numerical solver result or the success of the recovery
heuristic. Zero unused columns also pass the global comparison.

## Recovery from retained IPM multipliers

The previous min-z run retained raw numerical row multipliers. Twenty-five
slightly positive noise values, the largest2.524e-12, are set to zero.
Negative multipliers are represented by positive masses. On a1e10 grid,
face and Bellman masses are rounded down; final-budget and cap masses
are rounded up. This rounding proposes masses, not a proof by itself.

For an auxiliary state U_c(r), its capacity is the sum of mass times its
row coefficient. A negative capacity is repaired by decreasing masses on
Bellman rows whose destination is that state. Each such row has coefficient
-1 there and has a possible +1 only at the strictly smaller residual rank.
Process states in descending rank. The repair improves the current capacity
and may worsen only lower-rank capacities, which are processed later.
Different capped contexts do not interact through auxiliary states. This
produces nonnegative capacity at every potential without any optimizer.

Next repair a deficient t capacity by adding mass to its sole cap row t<=1.
All other nonface bonus coefficients are nonnegative. If a bonus column's
face demand exceeds its nonface capacity, decrease incident face masses
until the excess is gone. Every face coefficient is nonpositive, so this
can only decrease demands in all other bonus columns and cannot damage
auxiliary or t capacities. Finally divide all masses by at least their
total z load, also at least one. The result satisfies every dual column
inequality, which is then recomputed from scratch with Fractions.

The recovery changed205 Bellman rows, removing497/10000000000 total flow
mass. It added1/5000000000 to the t-cap mass, found1857 initially overloaded
bonus columns, and clipped1464 face masses. The normalizer is
1250000107/1250000000. These small adjustments are diagnostics; the complete
literal weak-duality check is the certified conclusion.

## New obstruction location

The largest final-budget masses occur in ordinary hub contexts. Example
representative exact-primal maximizers are

| Representative forest | Approximate budget mass |
|---|---:|
| A1+A5+A6+D7+I2(5) |0.07039|
| 2A1+2D6+D7 |0.05064|
| A1+3D5+H3+I2(5) |0.05043|
| 2A1+2D7+H3+I2(5) |0.04357|
| A1+A6+D7+H3+2I2(5) |0.04153|

The artifact records the exact masses, hub tuples, reserved components
when relevant, and32 such representatives. These are exact-primal
maximizers of the corresponding final-budget rows. They are not asserted
to be a decomposition of the repaired dual Bellman flow into individual
forests. Thus they identify contexts for further geometry without claiming
an unsupported exact forest-mass interpretation.

## Replay and bounded resources

Generation passed in session25113, final756735. A fresh solver-free replay
passed in session27307, finalaee553. Reproduce that replay with

  python3 -B scripts/fresh_compact21_complete_model_dual_repair_20261018.py --verify

It validates source hashes, rebuilds the full rational row iterator,
repeats the510 literal saved-forest coefficient matches, and checks all
global weak-duality columns. It imports no optimizer or numerical package.
The mathematical row/profile construction is shared with the reviewed full
model; the certificate is not an independent reconstruction of universal
geometric coverage.

The guarded generation took5.888 seconds with182255616 bytes process-tree
peak RSS, below the1GiB live cap. The worker took5.738 seconds and peaked
at156450816 bytes. Initial resources were196.84GiB free disk and8.69GiB
available memory. The parent enforced59 seconds, one thread, bounded output,
source ownership, and signal-safe process-group cleanup. The final artifact
is243817 bytes; both temporary paths were removed. No further optimizer
was launched after the preceding min-z run.
