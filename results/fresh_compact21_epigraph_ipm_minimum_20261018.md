# Complete compact21 payment model: exact feasible upper and numerical optimum

18 October 2026. This packet solves the complete conservative payment model
and checks one feasible rational vector. Its upper exceeds40; it does not
improve the compact arithmetic dimension ceiling.

## Formulation and outcome

The source reuses the complete Bellman formulation reviewed in
`compact21_full_forest_bellman_formulation_20261018.md`. It keeps all1979
native faces /102917 spectator cases (51443 exact duplicate classes),
all68 finite component types,729 exact hub contexts, and1659 Bellman states.
All510 saved forest rows are reconstructed literally before solving.
The proved D5 active/inactive rules, H4 funded transfer, and strengthened
A5 capacity one are included. The two decreasing physical forests are
excluded from ordinary budgets and replaced by the proved allowances;
all remaining forests are covered without a multiplicity cutoff.

The safe column reduction is the same as the preceding fixed-z attempt:
bonuses absent from negative native-face coefficients can be set to zero
because they have only nonnegative coefficients in other rows. All167865
nonface bonus coefficient occurrences are checked. Of93045 bonuses,
30599 remain and62446 are eliminated. The numerical matrix has79274 rows,
32260 columns, and378326 nonzero coefficients. The certificate stores the
entire compressed-to-global column map.

The single approved min-z HiGHS IPM run (session38559; final14d626) ended
optimally at numerical z=40.98190831955684 after87 IPM iterations and
zero crossover iterations. It used one computation thread and no new
library. Its rational proposal has

    z =16392763329/400000000 =40.9819083225,
    t =34404357/46510025.

Every native face inequality is checked over rational arithmetic. All
finite forests are re-evaluated by minimal exact Bellman states, discarding
the numerical potentials. The attaining forest is
A1+A5+A6+D7+I2(5), with ordinary budget status. No global bonus scaling
was needed. These checks establish this rational vector's feasibility
for the complete stated payment model, not exact numerical optimality.
The67122 retained numerical inequality duals are not an exact lower proof.
The30594 positive rational bonuses, original raw numerical bonus values,
and global row/column indices are retained for exact recovery.

## Verification and resources

A fresh solver-free replay passed (session91213; final78e69d). It checks
all source pins, rebuilds the complete native face system, repeats all510
literal forest-row matches, and recomputes the exact full forest maximum.
This replay uses the same mathematical row/profile implementation; an
independent full-forest coefficient reconstruction remains desirable.
The earlier root independent native-face evaluator supports that layer
but has not yet been rerun against this new vector.

The measured matrix build was3.257 seconds; the complete guarded solve
and exact postcheck took14.524 seconds. Combined process-tree peak RSS
was332611584 bytes (worker peak300564480 bytes), below the1GiB live cap.
Initial resources were196.89GiB free disk and6.40GiB available memory.
The numerical solve cap was39.7432 seconds; the parent deadline was59s.
The4.30MB final artifact is bounded below8MiB, stdout was787 bytes,
and both temporary certificate paths were removed. The owned process
group was terminated using the signal-safe cleanup helper.

## Limits

The exact upper is greater than40. There is no new dimension ceiling,
Lean verification, or exact optimality claim. In particular the previous
sub40 support/full-face candidates omitted complete forest constraints;
this packet resolves that omission but identifies a larger model budget.
An exact repaired dual can now locate the geometric obstruction without
another numerical optimization. Alternative incidence-dependent charge
rules are outside the present LP.
