# Dimension-22 A2-presence diagnostic

The [producer](../scripts/fresh_sqrt5_dimension22_a2_context_20261016.py) and [exact dual certificate](fresh_sqrt5_dimension22_a2_context_20261016_cap2.json) give a restricted-model lower bound greater than **44.11**, above the dimension-22 target **44**. Its approximate value is 44.118737006486; the exact rational value is checked during replay.

All nonnegative bonuses depend on exact entire-A1 count, I2(5) count capped at two, and presence of an entire A2 component. The nonnegative q2/H4 baseline and existing A1/I2(5) normalized cross families are unchanged. No A2 cross family is added here. The full common spectator counts `(a,b,c)` are enumerated. Vertex forests fix the exact counts and exclude A1, I2(5) and A2 from the residual choices. Face component detection uses intrinsic `I2(3)`; the corresponding finite-model name is `A2`.

There are 2,151 face rows, 120,711 spectator cases, 42,543 distinct inequalities and 29,884 variables. The exact dual rechecks all selected face and finite-forest inequalities, dual signs and column constraints. It proves only an obstruction for this conservative model, with no claim of optimality, nonexistence, or independent implementation review.

Generation and saved replay passed. Source and JSON are frozen.

```sh
python3 scripts/fresh_sqrt5_dimension22_a2_context_20261016.py --hub-cap 2 --verify
```
