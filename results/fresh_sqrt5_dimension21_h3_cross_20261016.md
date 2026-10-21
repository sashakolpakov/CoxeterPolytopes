# Dimension-21 H3 cross model: exact restricted obstruction

The [producer](../scripts/fresh_sqrt5_dimension21_h3_cross_20261016.py) and [rational certificate](fresh_sqrt5_dimension21_h3_cross_20261016_cap2.json) rebuild the frozen dimension-22 H3/I2(5)/A1 weight family in dimension 21. An exact dual certificate gives a lower bound greater than **40.25**, exceeding the required threshold **40**. This excludes this restricted weight family, not dimension-21 polytopes or richer weight rules.

The nonnegative q2/H4 baseline is unchanged. Bonuses depend on the exact A1 count, I2(5) count capped at two, and H3 presence. Connected bonuses and normalized A1, I2(5), and H3 crosses have the same definitions as in the frozen dimension-22 producer. All common external A1/I2(5)/H3 counts are enumerated; denominators use actual counts. All dimension-feasible singular rows of the completed finite queue are retained along with the inherited parameter-marked and dotted obligations.

The model has 1,979 face rows, 71,036 common-spectator cases, 33,353 distinct face inequalities, and 27,227 nonnegative bonus variables. The exact dual uses 1,128 constraints. Forest constraints are checked against actual rank-21 finite Coxeter forests; face constraints are rebuilt from the complete input model. Exact weak duality, rather than floating-point infeasibility, proves the restricted obstruction. No optimality claim is made.

The largest forest multipliers include `A1 + D6 + 2E7`, `D5 + 2E8`, `A1 + 2D7 + E6`, and `A1 + 2E8 + H4`. The heaviest face constraints couple E7/E8 and D-type corners, making exceptional-component contexts a concrete next enlargement to investigate.

Generation and saved replay passed. The source, JSON and report are frozen; no separate independent arithmetic review or dimension theorem is claimed.

```sh
python3 scripts/fresh_sqrt5_dimension21_h3_cross_20261016.py --hub-cap 2 --verify
```
