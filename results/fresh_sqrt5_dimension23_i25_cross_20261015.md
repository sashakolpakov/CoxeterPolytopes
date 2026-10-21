# Dimension-23 diagnostic for the same hub-context rule

The unchanged I2(5)-presence model does not reach the dimension-23 threshold `99/2 = 49.5` on the stated conservative catalogue. The [new source](../scripts/fresh_sqrt5_dimension23_i25_cross_20261015.py) and [exact dual certificate](fresh_sqrt5_dimension23_i25_cross_20261015_hub_context.json) prove a lower bound greater than **49.63**. This is a restricted-model obstruction, not a nonexistence theorem or a barrier to more general weights.

The model is the [dimension-24 construction](fresh_sqrt5_dimension24_i25_cross_20261015.md), rebuilt at rank 23 with its parity-correct threshold. It keeps the same nonnegative baseline, exact-A1 and I2(5)-presence contexts, and normalized A1/I2(5) cross families. All common spectator counts are enumerated again. There are 2,332 face rows, 54,993 face cases, 26,744 distinct inequalities and 10,605 variables.

The rational weak-duality certificate contains 287 terms, including 55 forest inequalities. Its forest contexts have I2(5) counts zero, one or two (17, 34 and four terms respectively). Each selected row is rebuilt from the exact face or finite-forest model, and all dual signs and column inequalities are checked over the rationals. Numerical optimization supplied the proposal; a rational normalization of forest coefficients and uniform downward adjustment of face multipliers produced the exact feasible dual. No optimality claim is made.

A concrete next distinction is to separate exactly one I2(5) component from two or more, or to add another normalized cross family. The current certificate does not rule out either. No such enlargement is tested here; the dimension-24 work and its independent review remain separate.

Generation and saved replay passed. Replay uses no numerical optimization:

```sh
python3 scripts/fresh_sqrt5_dimension23_i25_cross_20261015.py --split-presence --verify
```

No independent implementation review is claimed. The source, output and report belong to a new namespace; no dimension-24 or earlier proof file was changed.
