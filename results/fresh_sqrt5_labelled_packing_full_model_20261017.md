# Preserved rank13 model with labelled base transfers and integer packing

This fresh successor combines the [vertex-labelled shared capacities](fresh_sqrt5_vertex_refined_base_capacity_20261017.md) and [integer simplex-port packing](fresh_sqrt5_integer_simplex_port_packing_20261017.md). Their [exact union diagnostic](fresh_sqrt5_vertex_refined_integer_union_20261017.md) reduced the previous fixed-alpha obstruction to11.95359. That lower bound alone did not prove primal feasibility.

The sole `highs-ipm` solve terminated successfully in47.57seconds, within its300-second limit. Exact rational certificates give the bracket

\[
 12.908239002262698\ldots
 =\frac{2790968068837079040}{216216020508130303}
 \le T_*\le
 \frac{3227060041779}{250000000000}
 =12.908240167116.
\]

Thus target12 is impossible for this restricted-feature model. The earlier11.95359 was one fixed dual witness, not the optimized lower bound. This is a partial weighting-model obstruction, not a dimension bound. Generation and fresh saved verification passed; there was no second solve.

## Preserved and added domain

| Item | Exact count |
| --- | ---: |
| Face inequalities |82,806 |
| Inherited face obligations |82,530 |
| Extra actual-provenance copies |276 |
| Variables |350,744 |
| Active variables in the numerical solve |150,303 |
| Safely zero-only variables |200,441 |
| Full rank13 forest budgets |6,558 |
| Signed zero-budget variables |17,090 |
| Face-matrix nonzero coefficients |710,548 |

The source checks every inherited coefficient in every old row or its replacement copies. Each of the1,582 selected merged quad rows becomes all of its1,858 actual full-star provenance rows. Their inherited coefficients and demands are identical. No old obligation is dropped. All211,339 inherited positive-cost variables and all complete78-pair forest budgets remain.

The added positive-cost variables comprise59,962 mixed-three-port indicators and62,353 integer packing charges. Their coefficients are the complete exported features on the selected full-star classes, and **zero on other face classes**. This deliberate restriction is valid: dropping nonnegative payments only weakens coverage. It does not assume that unexported faces have no such incidence. All actual variants of every selected merged Q class are expanded consistently.

All5,693 inherited shared-transfer columns remain, including their conservative coefficients on unselected merged rows. The11,397 new labelled transfer columns use actual provenance coefficients on selected Q variants and retain their negative containing-simplex coefficients. They assign zero positive payment to other Q classes, again a valid weaker choice.

The4,969 actual unlabelled shared columns are omitted only after checking that each is exactly the sum of its labelled columns, including both positive provenance coefficients and negative simplex capacity. The earlier separate-dimension capacities are also dominated by shared capacity: `k*binomial(d+1,k) >= d+1`. No omitted redundant family supplies an additional restriction.

Every variable absent from all positive face entries is safely fixed to zero. Such a column either only spends vertex budget or has exclusively nonpositive signed entries. This reduction preserves the optimum. The target variable and all active charge variables are nonnegative.

## Exact certification

The [source](../scripts/fresh_sqrt5_labelled_packing_full_model_20261017.py) bounds the solver to one invocation. A saved invocation flag prevents accidentally restarting it. The saved status is success, objective12.908240166551062, four IPM iterations and elapsed47.568554seconds.

The returned primal charges were rounded to a common rational denominator. Exact face deficits were repaired by a common nonnegative increment to active inherited ordinary-pair and one-port variables only. Every face has a strictly positive coefficient sum in that set. Thus the repair increases every deficient face payment, and every full vertex budget was recounted exactly. Its exact maximum is the upper endpoint above, greater than12.

The returned face dual multipliers were rationalized and every signed-column deficit repaired by increasing its containing-simplex mass. Signed positive terms occur only on Q rows, so these repairs cannot create further signed deficits. All positive-cost column loads and all signed loads passed exact checks. This gives the strict lower endpoint above, rather than relying on a numerical infeasibility statement.

The separate [exact support diagnosis](fresh_sqrt5_labelled_packing_dual_support_20261017.md) finds5,638 positive face masses and758 positive forests. Although1,601 positive rows lie outside the newly exported feature classes, the restriction to the already exported classes alone still has exact lower bound12.6142151747, with no outside row restored by signed repair. Thus simply extending these same indicators onto the unexported classes cannot remove that restricted witness. Further structural incidence constraints remain a different possibility.

## Coverage and reproducibility

The231 untouched finite forests, fixed-discriminant stars and other missing faces remain outside the inherited model. Even an exact budget12 primal would be partial feasibility, with additional coverage required for any dimension theorem.

```sh
python3 scripts/fresh_sqrt5_labelled_packing_full_model_20261017.py --prepare
python3 scripts/fresh_sqrt5_labelled_packing_full_model_20261017.py --solve
python3 scripts/fresh_sqrt5_labelled_packing_full_model_20261017.py --verify
```

The [saved artifact](fresh_sqrt5_labelled_packing_full_model_20261017.json) stores the exact input hashes, every provenance-row map, new column-block maps, a digest of the complete reconstructed matrix and budgets, solver status, and any exact recovered certificates. The full inherited model is referenced rather than duplicated. Prepared model SHA256 is `85961c4af06135b227374dccc13ac58da3b071beeb59403175e55c7f96bde738`.

Generation30762 and fresh saved verification75913 both terminated successfully. The saved verification reconstructs the matrix digest and all coefficient identities, and checks both exact certificates without invoking the optimizer. Separate checks also confirmed that all saved dual masses are nonnegative and the primal target flag agrees with its exact maximum.

The root agent and construction_audit independently read the preservation, restricted-feature and exact-recovery logic and found no blocking issue. These are source/mathematical reviews, not independent solver runs. No frozen predecessor was changed.
