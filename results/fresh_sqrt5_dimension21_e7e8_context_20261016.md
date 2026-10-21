# Dimension-21 joint E7/E8 presence: exact restricted obstruction

The [producer](../scripts/fresh_sqrt5_dimension21_e7e8_context_20261016.py) and [rational certificate](fresh_sqrt5_dimension21_e7e8_context_20261016_cap2.json) add both E7 and E8 presence to the dimension-21 H3/I2(5)/A1 bonus contexts. An exact dual lower bound exceeds **40.26**, above the target **40**. This excludes this restricted weight family, not dimension-21 polytopes or a rule with richer contexts.

The existing nonnegative baseline and all three disjoint normalized cross families are unchanged. Exact A1 count, I2(5) count capped at two, and H3 presence are retained. Every common external count of A1, I2(5), H3, E7 and E8 is included. Vertex separation fixes their exact counts, excludes all five component types from residual finite-forest choices, and then restores their contributions. The same complete finite-queue, parameter-marked and dotted obligations are included.

The model has 1,979 face rows, 86,723 common-spectator cases, 43,438 distinct face inequalities and 108,908 bonus variables. The exact dual uses 2,490 face/forest constraints, each rebuilt during replay. The dual was proposed numerically and repaired by discarding selected face multipliers and shrinking the remaining face mass; rational weak duality verifies the resulting lower bound. There is no optimality claim.

Generation and saved replay passed. The source, JSON and report are frozen. No separate independent arithmetic review or improved dimension theorem is claimed.

```sh
python3 scripts/fresh_sqrt5_dimension21_e7e8_context_20261016.py --hub-cap 2 --verify
```
