# Exact partial-model result after cross-coface cap constraints

One bounded solve of the genuinely enlarged rank-13 model gives the exact bracket

\[
\frac{171523868647821690}{13513502427307591}
\le T_*\le
\frac{3173194699319}{250000000000},
\]

or **12.69277669282928 ≤ T* ≤ 12.692778797276**. The target 12 is still excluded within this specified partial model. This is not a dimension theorem or an obstruction to every possible geometric weighting method.

The [checker](../scripts/fresh_sqrt5_cross_coface_cap_full_model_20261017.py) and [saved model and certificates](fresh_sqrt5_cross_coface_cap_full_model_20261017.json) preserve the prior face domain and add new, independently reviewed geometric information. The numerical `highs-ipm` solve terminated successfully in 69.153 seconds and 16 iterations under a 300-second limit. There was exactly one optimizer invocation. Both a rational primal and a rational dual were recovered and checked; fresh saved replay passes without optimization.

## Preserved domain and additions

The successor retains coefficient by coefficient all **82,806 face obligations, 6,558 complete vertex budgets, 428,310 prior variables, and 94,481 prior signed transfers** of the [frozen full flag model](fresh_sqrt5_flag_local_packing_full_model_20261017.md).

It adds:

* **241 cross-coface cap capacities**, each with its actual quad-provenance coefficients and one negative containing-simplex coefficient. The [geometric package](fresh_sqrt5_cross_coface_cap_capacity_20261017.md) and [independent family review](fresh_sqrt5_cross_coface_cap_capacity_20261017_review.md) establish these capacities. They couple different isolated `A1` cofaces of one maximal simplex and therefore go beyond the previously exhausted within-one-coface packing cones.
* **217 saved positive cuts:** 31 from restricted local simplex packing, 24 strengthened three-port cuts, 47 four-port cuts, and 115 five-port cuts. Each rational column is multiplied by the same positive common denominator in its face entries and its vertex cost.

The resulting matrix has **428,768 variables, 162,942 active variables, 94,722 signed columns, and 803,682 face entries**. All inherited row indices, demands, budgets, coefficients, and common actual-provenance maps are checked. New positive charges have zero coefficients outside their exported full-star classes; this is valid underpayment, explicitly weaker than adding their entire geometric families.

## Exact certificates and replay

The rational primal satisfies every saved face inequality and vertex budget. Its maximum budget is the upper endpoint above. The rational dual has **3,447 positive face masses and 770 positive vertex masses**; direct exact checking of every saved positive and signed column gives the lower endpoint.

Signed dual recovery is monotone because all signed positive terms occur on quad rows and their negative term is on a simplex row. Raising simplex masses repairs those loads. Primal recovery raises only inherited positive repair variables, with positive entries explicitly checked in every face row. Inactive columns have no positive face entries and can be fixed to zero. The saved replay reconstructs the entire matrix, compares its digest and source hashes, and verifies both certificates.

The [independent source review](fresh_sqrt5_cross_coface_cap_full_model_20261017_review.md) checks the common provenance map, coefficient preservation, scaling, class restrictions, and recovery logic. It does not duplicate the numerical optimization or full matrix replay.

## What this establishes

The preceding [fixed-witness repair](fresh_sqrt5_cross_coface_cap_repair_20261017.md) lowered that witness's ratio to 10.2470841112 while retaining its complete local and joint checks. That inefficient repair did not establish feasibility. The new exact lower bound shows that a different dual witness survives in this enlarged finite-cut model.

The new optimized face masses have **not** been separated against the entire local simplex packing cone or every joint port inequality. Only the saved finite cuts stated above are part of this optimization. The missing complete geometric face domains also remain open. Thus neither the primal nor dual promotes the ambient compact arithmetic dimension bound. The complete unrestricted dimension-163 proof and the compact arithmetic dimension-21 proof are unaffected. No further optimizer run was launched.
