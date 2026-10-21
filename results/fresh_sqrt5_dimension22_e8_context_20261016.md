# Dimension-22 E8-count-context diagnostic

The [producer](../scripts/fresh_sqrt5_dimension22_e8_context_20261016.py) and [exact rational dual](fresh_sqrt5_dimension22_e8_context_20261016_cap2.json) give a restricted-model lower bound greater than **44.13**, above the dimension-22 target **44**. This is not a nonexistence theorem or a barrier to broader weight rules.

Every nonnegative bonus depends on the exact entire-A1 count, the I2(5) count capped at two, and the exact entire-E8 count (zero, one or two). The nonnegative q2/H4 baseline and the normalized A1 and I2(5) cross families are inherited without change. All common spectator counts `(a,b,e)` are enumerated; the complete finite-forest optimization fixes all three counts and excludes A1, I2(5) and E8 from the residual component choices. Every leftover rank zero or at least two is represented conservatively.

The model has 2,151 face rows, 52,650 spectator cases, 29,456 distinct face inequalities and 44,826 variables. The exact weak-duality certificate has 1204 terms, including 215 finite-forest inequalities. Its rational lower bound is approximately 44.136536820338; the saved replay checks the rational signs, column constraints, and comparison with 44.

Generation and saved replay passed. No independent implementation review or optimality claim is made. The source and JSON are frozen.

```sh
python3 scripts/fresh_sqrt5_dimension22_e8_context_20261016.py --hub-cap 2 --verify
```
