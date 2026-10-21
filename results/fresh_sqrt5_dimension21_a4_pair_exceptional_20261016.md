# Dimension-21 A4 pair bonuses with exact exceptional counts

The [producer](../scripts/fresh_sqrt5_dimension21_a4_pair_exceptional_20261016.py) and [rational certificate](fresh_sqrt5_dimension21_a4_pair_exceptional_20261016_cap2.json) combine independent A4 pair-orbit bonuses with exact E7/E8 component counts. An exact dual lower bound exceeds **40.07**, above the target **40**. This is a restricted-model obstruction only; it does not exclude dimension-21 polytopes.

The nonnegative baseline and normalized A1/I2(5)/H3/A4-to-other crosses remain unchanged. For `d≥2` A4 components, a marked pair in different components receives `z/(4(d−1))`, indexed by its unordered pair of rooted A4 orbits; unlisted pair types have weight zero. Its contribution per A4 component is one eighth of the sum over sixteen ordered root choices. Bonuses depend on exact A1/E7/E8 counts, I2(5) and A4 counts capped at two, and H3 presence. All denominators use actual counts.

Every common external six-tuple is included. Exact vertex counts of A1, I2(5), H3, E7, E8 and A4 are excluded from the residual finite-forest dynamic program and then restored. The completed finite-queue singular obligations and inherited parameter-marked/dotted obligations are retained.

The model has 1,979 face rows, 120,996 common-spectator cases, 63,525 distinct face inequalities and 472,694 variables. Its rational dual uses 5,218 reconstructed face/forest constraints. Generation and saved replay passed; there is no optimality claim or separate independent arithmetic review.

The source, JSON and report are frozen. No improved dimension bound is claimed.

```sh
python3 scripts/fresh_sqrt5_dimension21_a4_pair_exceptional_20261016.py --hub-cap 2 --verify
```
