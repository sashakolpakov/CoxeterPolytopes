# Dimension-21 exact E7/E8 counts: restricted obstruction

The [producer](../scripts/fresh_sqrt5_dimension21_e7e8_exact_20261016.py) and [rational certificate](fresh_sqrt5_dimension21_e7e8_exact_20261016_cap2.json) use exact numbers of entire E7 and E8 components in every bonus context, together with exact A1 count, I2(5) count capped at two, and H3 presence. An exact dual lower bound is greater than **40.23**, exceeding the threshold **40**. The result excludes this restricted weight family only.

Every common spectator count of the five distinguished component types is enumerated. Exact vertex counts are fixed and excluded from the residual finite-forest dynamic program before their component contributions are restored. The nonnegative q2/H4 baseline and normalized A1/I2(5)/H3 cross rules remain unchanged. The complete finite-queue singular obligations and inherited parameter-marked/dotted obligations are retained.

The model has 1,979 face rows, 86,723 common-spectator cases, 43,606 distinct face inequalities, and 190,589 variables. The exact dual uses 2,583 reconstructed face/forest constraints. Its rational weak-duality verification does not assert optimality. Presence and exact-count diagnostics therefore both leave the dimension-21 exclusion unresolved.

Generation and saved replay passed. The source, JSON and report are frozen. There is no separate independent arithmetic review and no improved dimension claim.

```sh
python3 scripts/fresh_sqrt5_dimension21_e7e8_exact_20261016.py --hub-cap 2 --verify
```
