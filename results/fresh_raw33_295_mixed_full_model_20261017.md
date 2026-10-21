# Finite model with the current union and raw33/295 capacities

The enlarged finite model has exact primal and dual certificates placing its optimum in

    [1338621568130792400 / 108108118477533743,
     6191130293139 / 500000000000]
        = [12.382248317539402..., 12.382260586278].

Both endpoints exceed 12. The new capacities improve the preceding finite model, but do not achieve its target-12 feasibility condition. The earlier 11.779623678658464 calculation repaired one fixed witness; optimizing the enlarged model produces different face masses. No dimension theorem follows from either calculation.

**FINAL/FROZEN:** generation, exact recovery and fresh complete-matrix replay all pass. This successor is outside the preceding 24-package checkpoint.

Final preparation 66507 and the sole solve/exact recovery 50219 terminate exit 0/PASS. Preparation checks each input package's recorded source hashes directly. The [source](../scripts/fresh_raw33_295_mixed_full_model_20261017.py) rebuilds the complete frozen [current-raw23 model](fresh_raw23_mixed_full_model_20261017.md), verifies its prepared identity, and appends every column from the frozen [current220 union](fresh_current220_actual_corner_union_20261017.md) and [raw33/295 mixed package](fresh_raw33_295_mixed_actual_corners_20261017.md). The [model metadata and certificates](fresh_raw33_295_mixed_full_model_20261017.json) record exactly one solver invocation. No large matrix copy is saved.

| Prepared quantity | Count |
|---|---:|
|Inherited variables, retained literally|438266|
|New export tuples / distinct added tuples|795 / 795|
|Total variables|439061|
|Active variables|173235|
|Face obligations|82806|
|Vertex budgets|6558|
|Zero-budget signed columns|96176|
|Positive-vertex-cost signed columns|52|
|Face-matrix nonzeros|826459|

The matrix identity is

    641d6291f9eaef6c29d6c6e8476076e0f03c90cba8a39ead2ec936b953f52573.

The two added inputs contain 443 and 352 columns. Their 795 vertex/cost/face-coefficient tuples are distinct. Positive denominator scaling applies equally to each face coefficient and its vertex cost. Every inherited variable, face demand, face coefficient, vertex-budget coefficient and zero-budget transfer is checked unchanged. The 52 signed-cost columns remain ordinary vertex budgets; signed face entries do not make their nonnegative primal variables free. Exact recovery continues to use the inherited strictly positive old pair/port columns, and the appended verifier explicitly recounts all new dual loads.

The two input capacity packages share exactly their 4,431 baseline face masses. Preparation verifies the full sequential allowance chain from the optimized parent through current220 to raw33/295. Their fixed-witness ratio justified this changed finite target-12 problem; it did not predict the optimum.

HiGHS-IPM reports numerical optimality at 12.382260585718997 after 89.84 seconds and 83 iterations, within its one 300-second limit. The numerical matrix has 89,364 rows, 173,236 columns including the objective variable, and 910,076 nonzeros. Exact rational recovery produces the bracket above. The primal has 23,356 saved nonzero rounded variables and uses repair delta 3/500000000000; the dual has 4,713 positive face masses and 986 vertex allowances. All original and new matrix constraints are verified by the exact checker, including all 795 explicitly recounted new loads. The target flags correctly report primal failure to meet 12 and a dual obstruction above 12.

All earlier exported constraints remain represented. Complete implicit local-packing and joint-family certificates at older masses are not certificates for these new masses. Their features and marginal obligations require reconstruction or an explicit valid comparison before a full-family claim. Coverage of omitted geometric classes and arithmetic fields is also separate. The dual is an exact obstruction to this finite charge matrix, not a globally realizable polytope.

All ten direct prepared bindings, five local links and source syntax pass. Root's full append-source/report read passes preservation, scaling, sequential allowance consistency and recovery scope, without claiming an independent matrix reconstruction. A separate no-producer saved check verifies nonnegative/indexed dual and primal data, both exact ratios and target flags, all 795 new loads, and the maximum of the saved primal vertex budgets. This small check does not independently reconstruct the full inherited matrix. Fresh complete-matrix replay 6358 terminates exit 0/PASS, rebuilding the prepared identity and verifying both exact certificates. There are no live jobs belonging to this primary package and no second optimizer was invoked.

    python3 scripts/fresh_raw33_295_mixed_full_model_20261017.py --prepare
    python3 scripts/fresh_raw33_295_mixed_full_model_20261017.py --solve
    python3 scripts/fresh_raw33_295_mixed_full_model_20261017.py --verify
