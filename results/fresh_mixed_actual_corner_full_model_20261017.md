# Enlarged finite model with mixed actual-corner capacities

One bounded HiGHS solve, followed by exact rational recovery and a fresh saved replay, gives

\[
\frac{661061434478917472}{52852828351105989}
\le \operatorname{opt}_{\rm finite}\le
\frac{6253797331083}{500000000000},
\]

or **12.50758862113165 ≤ opt ≤ 12.507594662166**. The target is 12, so this enlarged finite model still has an exact obstruction above the target. No new dimension theorem follows. The source is [here](../scripts/fresh_mixed_actual_corner_full_model_20261017.py), and the compact model description and both exact certificates are [here](fresh_mixed_actual_corner_full_model_20261017.json).

This is a new optimized witness. The earlier [mixed-corner repair](fresh_mixed_raw23_raw25_corners_20261017.md) lowered a different, fixed witness from 12.38203 to 11.86346. That crossing justified this new solve; it did not imply feasibility. The optimizer shifts the face masses and obtains the obstruction recorded above.

## Preserved matrix and new columns

The source reconstructs the frozen [actual-corner enlarged model](fresh_actual_corner_enlarged_model_20261017.md), verifies its saved prepared description, and preserves its complete coefficients:

| Object | Count |
|---|---:|
| Face obligations | 82,806 |
| Complete rank-13 vertex budgets | 6,558 |
| Inherited variables | 430,666 |
| Inherited zero-budget signed transfers | 96,176 |
| Inherited signed vertex-cost incidence columns | 52 |
| Exported raw25 columns | 500 |
| Exported mixed raw23/raw25 columns | 1,251 |
| Distinct appended columns | 1,355 |
| Exactly duplicated new tuples merged | 396 |
| Total variables | 432,021 |
| Active variables | 166,195 |
| Face-matrix nonzeros | 809,767 |

Every column from the [raw25 package](fresh_raw25_actual_corner_compatibility_20261017.md) and the mixed package is represented. Deduplication merges only an identical tuple consisting of the vertex index, positive vertex cost and complete positive face-coefficient list. The list of original source columns is retained for each merged column. No proportional, approximately equal or geometrically guessed columns are merged, and no inherited variable is removed.

Each new column is a nonnegative payment variable with positive vertex cost. A common positive denominator multiplies its complete face coefficient list and its cost together. All face indices and vertex indices are checked in the inherited model. The additions do not modify any negative coefficient. In particular, the 52 forced-quadrilateral incidence columns with a negative simplex term and a positive vertex cost remain ordinary vertex-budget columns; they are never treated as zero-budget transfers.

The model digest is

`fb780e5f2aa2e16222b58bbd859001b6a4899056891c2c581df00f2321e88967`.

The saved artifact contains append metadata and source hashes, rather than a duplicate large matrix. A replay rebuilds the matrix and checks the digest and every inherited coefficient before checking the certificates.

## Solver and exact validation

There was exactly one solver invocation: `highs-ipm` with a 300-second limit. It terminated successfully in 54.30 seconds with numerical objective 12.507594661601228. This numerical value is not the certificate.

The exact primal rounds nonnegative variable values rationally and repairs any deficits by increasing only the inherited strictly positive pair/port variables. Every face payment and all 6,558 vertex budgets are then recounted exactly; the largest exact budget is the displayed upper bound. Its `meets_target12` flag is false.

The exact dual rationalizes nonnegative face masses, raises simplex masses to meet the inherited zero-budget Q-positive/S-negative transfers, and recomputes every vertex-cost load. This includes all 52 signed vertex-cost columns and all 1,355 new columns. Every zero-budget load must be nonpositive, and every other load must be at most its cost times its vertex allowance. The objective divided by the sum of the exact vertex allowances is the displayed lower bound. Its `excludes_target12` flag is true. Indexing, uniqueness and positivity of all saved masses are checked by the inherited standalone verifier.

Preparation, solve/recovery and a fresh `--verify` replay all finished with exit 0. The final replay status is `PASS_SAVED_MIXED_ACTUAL_CORNER_FULL_MODEL PASS_EXACT_MIXED_ACTUAL_CORNER_MODEL_CERTIFICATES`. The underlying corner packages have their own complete generation/reconstruction and scoped mathematical reviews; no second optimizer or independent duplicate large LP was run here.

## Changed support and remaining coverage

The new exact dual has 4,661 positive face masses and 780 positive vertex allowances. Its 1,908 quadrilateral rows comprise 433 rows in the inherited explicit-provenance map and 1,475 remaining merged rows. The other 2,753 rows are simplex or maximal-edge obligations. Relative to the preceding optimized dual, 1,300 rows become positive and 1,571 disappear.

The mixed exported cohort still contains 56 positive raw23 rows and 48 positive raw25 rows. Thirteen newly positive pure raw23/raw25 rows are outside that selected cohort: seven raw23 and six raw25, including the returned raw23 rows 26457 and 26465. This is a support diagnostic, not a claim that these omissions explain the obstruction. Other families also carry substantial mass: the 31 pure raw39 rows have total alpha 293940597/10000000000.

Most importantly, the complete local simplex-packing covers and millions of joint-family checks certified at the old fixed face masses do not certify these new face masses. The finite model retains every previously exported necessary column, but it does not contain every member of those implicit families. New full-family separation would require reconstruction at the new masses or a valid marginal comparison; it has not been performed in this package.

The inherited geometric domain is partial: it omits 231 rank-13 finite forests, fixed-discriminant classes and other face families. The new actual-corner payments apply to their stated selected full-star classes and use the fixed ground field K = Q(√5). Neither complete face/field coverage nor simultaneous compact ambient realization is claimed. Even an exact target-12 primal for this finite model would require those remaining global obligations before an ambient dimension conclusion.
