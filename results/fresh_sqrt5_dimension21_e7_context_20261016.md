# Dimension-21 E7 context: exact restricted obstruction

The [producer](../scripts/fresh_sqrt5_dimension21_e7_context_20261016.py) and [rational certificate](fresh_sqrt5_dimension21_e7_context_20261016_cap2.json) add the presence of an entire E7 component to every bonus context in the dimension-21 H3/I2(5)/A1 model. An exact dual lower bound is strictly greater than **40.36**, exceeding the target **40**. This rules out this restricted weight family; it is not a nonexistence theorem for dimension-21 polytopes or richer weights.

The existing nonnegative baseline and all three disjoint normalized cross families are unchanged. Contexts retain the exact A1 count, I2(5) count capped at two, and H3 presence. Every common external E7 count is included along with the three earlier spectator counts. Vertex separation fixes the actual counts of all four distinguished types, removes them from the residual finite-forest dynamic program, then restores their contributions. Capped bonus contexts never replace actual normalization denominators or rank counts.

There are 1,979 face rows, 80,416 common-spectator cases, 39,221 distinct face inequalities and 54,454 bonus variables. The input retains the dimension-feasible singular obligations from the completed finite queue and the inherited parameter-marked and dotted obligations. The exact dual uses 2,011 constraints. Saved replay reconstructs each used face or rank-21 forest constraint and verifies rational weak duality. No optimality claim is made.

Generation and saved replay passed. This producer package has no separate independent arithmetic review and makes no improved dimension claim. The source, JSON and report are frozen.

```sh
python3 scripts/fresh_sqrt5_dimension21_e7_context_20261016.py --hub-cap 2 --verify
```
