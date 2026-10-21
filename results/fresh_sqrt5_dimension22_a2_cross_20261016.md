# Dimension-22 normalized A2-hub cross diagnostic

The [producer](../scripts/fresh_sqrt5_dimension22_a2_cross_20261016.py) and [exact dual certificate](fresh_sqrt5_dimension22_a2_cross_20261016_cap2.json) give a restricted-model lower bound **greater than 44.05**, above the dimension-22 target **44**. The rational bound is approximately 44.052810277694.

This enlarges the A2-presence context model with normalized A2-to-component crosses. At a vertex with `c` entire A2 components, the new pair weight is `w/(2c)`; destination components exclude A1, I2(5) and A2. Thus these pairs are disjoint from the existing normalized A1 and I2(5) crosses. Summing over all two-root source components contributes `w` per destination root and preserves additivity for fixed exact counts. All nonnegative bonus variables depend on exact A1 count, I2(5) count capped at two, and A2 presence. The fixed nonnegative q2/H4 baseline is unchanged.

Every common spectator count `(a,b,c)` is included, and the finite-forest dynamic program fixes all three exact counts. Face A2 detection uses intrinsic `I2(3)`, while the finite-model name is `A2`. There are 2,151 face rows, 120,711 spectator cases, 42,703 distinct inequalities and 29,935 variables. Only 51 new variables occur relative to the context-only model, reflecting how few unpaid marked corners expose an entire A2 source.

Generation and saved rational replay passed. The dual verifies every selected face/forest row, all signs and all column inequalities. There is no optimality, independent-review, or geometric nonexistence claim. Source and JSON are frozen.

```sh
python3 scripts/fresh_sqrt5_dimension22_a2_cross_20261016.py --hub-cap 2 --verify
```
