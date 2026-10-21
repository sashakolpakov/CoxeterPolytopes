# Finite model enlarged by raw24 mixed actual-corner capacities

The single bounded solve recovers exact certificates placing this finite model's optimum between **12.457543530716451** and **12.457547514692**, still above12. The [saved certificate](fresh_mixed_raw24_full_model_20261017.json) contains both rational certificates. The previous fixed-witness repair11.0634964874 is therefore not a feasible target12 rule: optimization moves to different face masses.

This successor preserves every exported obligation and column of the [frozen mixed-raw39 model](fresh_mixed_raw39_full_model_20261017.md). It adds78 [local-simplex charge columns](fresh_mixed_raw39_full_model_20261017_local_cones.md),68 [selected joint-port columns](fresh_mixed_raw39_full_model_20261017_joint_ports.md),380 [within-raw24 capacities](fresh_raw24_actual_corner_compatibility_20261017.md) and1,793 [mixed raw24/raw23/raw25/raw39 capacities](fresh_mixed_raw24_actual_corners_20261017.md). All2,319 input exports are represented as2,040 distinct additions after merging279 exactly identical new vertex/cost/face-coefficient tuples. The [source](../scripts/fresh_mixed_raw24_full_model_20261017.py) stores only sparse append metadata and exact certificates; the full inherited matrix is reconstructed in memory.

The new mixed geometry—not a repeated solve of an unchanged model—motivates this step. It lowered the preceding fixed witness from12.417934 to11.063496; omitting the genuinely raw24-mixed columns left that witness at12.417927. This supplied new necessary constraints but did not determine the enlarged model's optimum.

## Preservation and exact recovery

All433,753 inherited variables,82,806 face obligations,6,558 vertex budgets and96,176 zero-budget signed transfers retain their exact coefficients. The52 signed incidence columns with positive vertex cost retain their negative simplex entries and positive budget cost. They are not converted into zero-budget transfers.

Each added column has positive face entries and positive vertex cost. The source scales every entry and cost by the same positive common denominator. Restricting a nonnegative payment to its exported face classes, with zero elsewhere, is valid underpayment. The build checks the common input face-mass vector and the successive vertex-mass interfaces, compares every old row and budget coefficient, preserves all metadata and checks that each face still has a strictly positive old repair coefficient.

The inherited one-shot solver uses HiGHS-IPM with a300-second limit. Numerical output only seeds exact recovery. Primal recovery rounds rationally and repairs deficits using old pair/port variables whose coefficients are strictly positive. Every face and budget is then checked exactly. Dual recovery raises simplex masses to meet zero-budget transfers and recomputes every vertex-cost load, including all52 signed incidence columns and every addition. The new verifier also explicitly recounts the appended column loads. Fresh saved verification reconstructs the matrix and checks the certificates without another optimizer invocation.

The full exported matrix has435,793 variables,169,967 active variables,82,806 face obligations,6,558 budgets and817,197 face entries. The solver matrix has89,364 rows,169,968 columns including its objective variable, and897,546 nonzero entries. Its exact model hash is

    65e1790b7a9cae086d2829dc6084faf860cd2aed4924b5cb23a5d319df72d6fb.

The exact finite-model bounds are

    5755387047017860 / 462000155393947
        <= T_* <=
    3114386878673 / 250000000000.

Preparation68369 and solve95324 both terminate exit0/PASS. HiGHS-IPM is invoked exactly once and returns optimal numerical status0 in51.752 seconds and14 iterations. Exact recovery checks both saved certificates. The new dual has4,647 positive face masses and807 positive vertex allowances. Fresh saved verification26453 terminates exit0/PASS, rebuilding the full matrix and rechecking both certificates without another optimizer invocation. Root also independently reads the full new append source and mathematical composition with no identified gap; that review does not duplicate the full matrix or solver replay.

## Coverage limits

This is a finite partial model in dimension13 with target12. It contains every earlier exported constraint and all2,319 stated new exports. It does not contain the entire infinite family of local charges, every implicit joint-port inequality, or complete missing face/forest/field coverage. The preceding full local-packing mixtures and3,318,122 joint checks apply to their saved5,059 face masses. They cannot be reused at this new4,647-face dual without reconstructing the marginals or proving a valid comparison.

An exact dual above12 obstructs the stated finite charge family; it does not describe a simultaneous ambient Coxeter-polytope realization. A target12 primal, if one were found in a future enlarged model, would still require complete global coverage before a dimension theorem. The [separate support diagnostic](fresh_mixed_raw24_full_model_20261017_support.md) identifies where this new witness leaves the previously exported geometric cohorts and gives exact restrictions without another optimizer. No repeated optimization is performed here.
