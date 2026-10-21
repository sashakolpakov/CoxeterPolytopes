# Exact obstruction to the current dimension-24 weight model

The normalized-A1 model cannot attain the required dimension-24 threshold 54 on its complete input catalogue. An exact rational dual gives a lower bound **111/2 = 55.5** for its maximum vertex weight. This is an obstruction to that weight model, not to compact polytopes or to more general weight rules.

The [source](../scripts/fresh_sqrt5_dimension24_weight_barrier_20261015.py) imports the frozen complete-face model read-only and rebuilds it at rank 24. The [certificate](fresh_sqrt5_dimension24_weight_barrier_20261015.json) records 21 nonzero dual terms: 18 distinct face inequalities and three vertex-forest inequalities. These reference 20 literal face rows because two pairs of rows have identical linear constraints. All coefficients and weak-duality slack are checked as rational numbers; numerical optimization only proposes the dual. No optimum is claimed.

The three spherical forests, each used with dual coefficient one third, are

- `E8 + E8 + E8`;
- `A1 + I2(5) + E7 + E7 + E7`;
- `A1 + I2(5) + D7 + D7 + D7`.

Their baseline totals are 42, 37, 37. The weighted face constraints force their average total weight to be at least 111/2. The dual involves 35 bonus variables, with 11 nonzero nonnegative variable slacks. The complete rebuilt model has 2,516 face rows, 11,573 distinct face inequalities and 5,196 bonus variables. It retains every dimension-feasible singular stopping case from the completed 5,172-state finite queue, alongside the older conservative input families.

The scope of the obstruction is precise: the q2/H4 baseline is fixed; connected-pair bonuses are nonnegative and depend on the entire connected marked-pair orbit and exact number of A1 components; A1-to-component-root bonuses are divided by the number of A1 components; all other disconnected marked pairs have weight zero. The certificate does not cover arbitrary forest-dependent weights, non-A1 cross weights, or lowering the baseline. Some input rows are necessary principal-star obligations rather than separately certified actual quadrilaterals, so this is also not a lower bound on every possible geometric charging scheme.

The sparse dual identifies a concrete enlargement to test: allow an entire I2(5) component to pay a cross angle with a root in another component. Dividing its rooted bonus by twice the number of I2(5) components makes the sum additive over roots in the other components. This changes the three-forest obstruction while preserving nonnegative baseline payments. Its global spectator contexts and vertex maximum require a separate calculation; no feasibility result is asserted here.

Generation and saved exact replay passed. The replay rebuilds the model and verifies each selected face/forest row, every dual sign, every column inequality, and the rational lower bound. No independent implementation review is claimed:

```sh
python3 scripts/fresh_sqrt5_dimension24_weight_barrier_20261015.py --verify
```
