# Dimension-21 A4 crosses with exact exceptional counts: obstruction

The [producer](../scripts/fresh_sqrt5_dimension21_a4_exceptional_20261016.py) and [rational certificate](fresh_sqrt5_dimension21_a4_exceptional_20261016_cap2.json) combine the normalized A4 cross family with exact E7/E8 counts. An exact dual lower bound exceeds **40.14**, above the dimension-21 threshold **40**. This is an obstruction to the stated weight family only.

The existing nonnegative baseline and normalized A1/I2(5)/H3 crosses remain unchanged. A4-to-other pairs receive `v/(4d)` outside the earlier hub types, and roots in distinct A4 components receive `(u_r+u_s)/(4(d−1))`. Bonuses depend on exact A1/E7/E8 counts, I2(5) and A4 counts capped at two, and H3 presence. Actual counts are used in every denominator and rank calculation.

All common external six-tuples are enumerated. Exact vertex counts of A1, I2(5), H3, E7, E8 and A4 are fixed, excluded from the residual finite-forest dynamic program, and then restored. The complete finite-queue singular obligations and inherited parameter-marked/dotted obligations are retained.

The model has 1,979 face rows, 120,996 common-spectator cases, 63,525 distinct face inequalities and 472,694 variables. Its exact dual uses 5,109 reconstructed face/forest constraints. Numerical multipliers were rounded to a common denominator, selected face terms were discarded, and the remaining face mass was shrunk before rational weak-duality verification. No optimality claim is made.

Generation and saved replay passed. The source, JSON and report are frozen. No separate independent arithmetic review or improved dimension bound is claimed.

```sh
python3 scripts/fresh_sqrt5_dimension21_a4_exceptional_20261016.py --hub-cap 2 --verify
```
