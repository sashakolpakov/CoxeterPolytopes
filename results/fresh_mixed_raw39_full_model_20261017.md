# Finite model with mixed raw39 actual-corner capacities

The one bounded solve recovered exact certificates placing this finite model’s optimum between **12.48512744950078** and **12.485147275642**, still above12. The [saved certificate](fresh_mixed_raw39_full_model_20261017.json) contains both rational certificates. This supersedes the earlier fixed-witness calculation as the optimization result; it does not turn the fixed-witness value11.58583 into a feasible target12 rule.

This fresh successor preserves the complete exported matrix of the [previous mixed actual-corner model](fresh_mixed_actual_corner_full_model_20261017.md). It adds the132 saved exact local-simplex charges,183 selected joint-port charges,129 within-raw39 capacities and1,440 [mixed raw39/raw23/raw25 capacities](fresh_mixed_raw39_actual_corners_20261017.md). All1,884 exported input columns are represented as1,732 distinct additions after merging152 exactly identical new vertex/cost/coefficient tuples. The [source](../scripts/fresh_mixed_raw39_full_model_20261017.py) stores only new-column metadata and exact certificates, rebuilding the inherited sparse matrix in memory.

The reason for this solve is new geometric information: the mixed raw39 family lowered the former fixed-witness ratio from12.3820705641 to11.5858256266. Without the genuinely raw39 mixed columns, that witness still gave12.3793273509. This motivates a different finite model; it is not evidence that a globally valid target12 rule already exists.

## Matrix and recovery

Every inherited432,021 variable,82,806 face obligation,6,558 vertex budget and96,176 zero-budget signed column retains its exact coefficients. The52 signed incidence columns that spend a positive vertex budget retain both their negative simplex entries and their vertex costs; they are not treated as free transfers.

Each new input column has nonnegative face entries and a strictly positive vertex cost. A common positive denominator scales every entry and its cost together, producing an integer matrix. Positive features may be restricted to their explicitly exported face classes by assigning zero elsewhere; this is valid underpayment. The build checks the common input face masses and successive vertex-mass interfaces, preserves all old row metadata and coefficients, and verifies that every face still has a strictly positive old repair coefficient.

The inherited one-shot solver uses HiGHS-IPM with a300-second limit. Its output is only a numerical starting point for exact recovery. Primal repair changes strictly positive old pair/port variables, then checks every saved face and budget rationally. Dual recovery repairs the zero-budget signed constraints by increasing simplex masses and recomputes every positive-cost vertex load, including the52 signed incidence columns and all newly appended charges. The producer explicitly rechecks the new column loads. Fresh saved verification rebuilds the matrix and verifies both certificates without another optimizer invocation.

## Scope and status

Preparation39795 and solve41311 both terminated exit0/PASS. HiGHS-IPM returned optimal numerical status0 in93.205 seconds and34 iterations. It was invoked exactly once with the stated300-second limit. Exact recovery then checked both certificates.

The full exported matrix has433,753 variables, of which167,927 are active,82,806 face obligations and6,558 budgets. There are813,374 face entries and891,683 nonzero entries in the complete solver matrix. All96,176 zero-budget signed columns and52 signed positive-cost columns remain. The model hash is `91edd9e5ec2ce9e039b37ec8826d1ddb10c5ba2b3d1cf1ce90e3d0c0cda1c227`.

The exact finite-model bounds are

\[
\frac{44991485646505450}{3603606437217783}
\le T_*\le
\frac{6242573637821}{500000000000}.
\]

The new dual has5,059 positive face masses and1,068 positive vertex allowances:2,886 simplex/edge rows and2,173 quad rows. Of those Q rows,461 use the explicit provenance map and1,712 remain merged. Exactly117 positive Q classes in the new dual lie in the exported mixed cohort:61 raw23,49 raw25 and7 raw39. The other2,056 positive Q rows do not receive these new mixed-family coefficients. Compared with the preceding fixed support,1,341 rows are newly positive, including640 simplex/edge rows and701 Q rows. These are support diagnostics, not a claim that those classes admit simultaneous ambient realization.

Fresh saved verification20085 terminated exit0 with `PASS_SAVED_MIXED_RAW39_FULL_MODEL`. It reconstructed the matrix and recounted every certificate inequality without invoking an optimizer. Direct source hashes, Python syntax and local links pass. No frozen predecessor was modified. The inherited recovery engine and earlier append conventions retain their previously recorded review scopes; this package’s full matrix and certificate replay was performed by the producer.

This is a finite partial model in dimension13, with target12. It retains all previously exported obligations and adds every stated exported column. It does not include every possible local packing charge, every implicit joint-family inequality, or all missing face and arithmetic-field contexts. The earlier complete local packing covers and selected joint scans apply to their saved face masses; they cannot be reused silently at a newly optimized dual. A new exact dual is an obstruction to this finite charge family, not a simultaneous realization by an ambient Coxeter polytope. A target12 primal would still require complete global coverage before a dimension theorem.
