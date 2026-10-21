# Finite model with current raw23 mixed actual-corner capacities

The single enlarged solve gives exact certificates placing the represented finite-model optimum between **12.391408496497405** and **12.391413504106**, both above target 12. This improves the preceding bracket 12.41451686839335–12.414523375499. The [source](../scripts/fresh_raw23_mixed_full_model_20261017.py) preserves the full [previous finite model](fresh_raw178_296_297_full_model_20261017.md) and appends all exported columns from the [raw87/88/89 cohort](fresh_raw87_88_89_actual_corners_20261017.md), [mixed larger-core cohort](fresh_mixed_raw87_large_core_corners_20261017.md), and [current raw23 union](fresh_raw23_mixed_current_corners_20261017.md).

The preceding fixed-witness repair lowered its ratio to 11.56686352540649. That removed one particular dual obstruction and justified testing the genuinely enlarged matrix; it did not imply a feasible target-12 primal. The new exact dual supplies a different obstruction above 12. No dimension theorem follows from this partial model.

## Exact append and preservation

All 436,262 inherited variables, 82,806 face obligations, 6,558 vertex budgets and 96,176 zero-budget signed transfers retain their exact coefficients and metadata. The 52 signed positive-cost incidence columns retain both their negative simplex coefficients and their vertex cost. Every variable is nonnegative: “signed” describes face coefficients, not a free variable.

The three packages export 1,002, 263 and 739 capacity columns respectively. The append checks their common 4,410 face masses, their consecutive initial/final vertex allowances, and their exact repaired ratios. All 2,004 input tuples are distinct. A tuple includes its vertex, positive cost and full face coefficients. Each column's cost and coefficients are scaled together by one positive denominator; no unscaled coefficient is mixed with a scaled budget. The full inherited provenance-to-face map is preserved.

The resulting model has 438,266 variables, of which 172,440 are active, and 823,033 face-matrix entries. It retains the same 82,806 faces and 6,558 budgets. Its exact matrix hash is

    ce5e13fc67cc4c46f68293dcb697f14a56eb00077a5e30f08905c3ad0dd27d47.

Every old row and budget is compared coefficient by coefficient with the inherited matrix. Every face retains strictly positive old pair/port coefficients for exact primal repair. Complete matrix data is reconstructed in memory; only sparse append metadata and certificates are saved in the [JSON](fresh_raw23_mixed_full_model_20261017.json).

## Solve and exact recovery

The inherited HiGHS-IPM engine is called once with a 300-second limit. Its numerical result seeds rational certificates. Primal recovery rounds nonnegative variables and repairs all deficits using strictly positive old columns, then checks every face and budget exactly. Dual recovery raises simplex masses for zero-budget transfers and recomputes every positive-cost load, including all 52 signed-cost columns. The successor separately recounts every new capacity load. The saved replay rebuilds the entire matrix and verifies both rational certificates without another optimization.

The exact finite-model bracket is

    893073974017141280 / 72072030735616573
        <= T_* <=
    6195706752053 / 500000000000.

FINAL/FROZEN: preparation 23420, the sole solve 32174 and separate fresh saved verification 15117 all terminate exit 0/PASS. HiGHS reports optimal status 0 after 85.2466 seconds and 15 iterations, with numerical objective 12.391413503532384. Its matrix has 89,364 rows, 172,441 columns including the objective variable, and 905,855 nonzero entries. Exact recovery returns a primal with 23,393 nonzero rounded variables and a dual with 4,431 positive face masses and 761 positive vertex allowances. Both pass the complete inherited exact verifier and the successor's new-load recount. A separate no-producer small check confirms nonnegative masses/variables, valid unique indices, all 2,004 unscaled new loads, exact ratios and both target flags. All 13 direct source bindings and source syntax pass. The complete matrix is reconstructed for fresh verification, without another optimizer invocation.

Root's complete new append-source read passes preservation and the sequential exact witness chain. Construction's separate full source/delta read passes the common scaling, tuple deduplication, inherited metadata and coefficient preservation, signed-cost treatment, unchanged recovery wiring and added load checks. These are source/composition reviews; neither independently rebuilds the matrix or duplicates the solver/recovery arithmetic. Root also independently read the frozen raw23 producer and relevant pair-field helpers, finding no gap in the native field selection, complete recursive reuse, actual occurrence coefficients and conservative capacities. That review does not duplicate its geometric census or contact computations.

## Coverage

This is the correct target-12 feasibility question for the enlarged finite exported matrix. The new actual-corner constraints are valid with zero coefficients outside their exported classes. All prior exported constraints remain present. However, the model does not contain every implicit local-simplex packing or joint-port inequality, and missing global face, forest and arithmetic-field coverage remains unfinished. Full covers at old face masses cannot be reused at newly optimized masses without reconstructing their marginals or proving an appropriate comparison.

An exact dual above 12 obstructs this finite charge family; it does not assert simultaneous ambient realization of its face frequencies. A target-12 primal would still need the missing coverage before a dimension theorem. The separate support diagnostic classifies the new exact dual and gives exact restrictions without another solve.
