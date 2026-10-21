# Dimension-21 A4 pair-orbit bonuses: restricted obstruction

The [producer](../scripts/fresh_sqrt5_dimension21_a4_pairs_20261016.py) and [rational certificate](fresh_sqrt5_dimension21_a4_pairs_20261016_cap2.json) allow independent bonuses for the occurring unordered pairs of rooted A4 orbit types. The exact dual lower bound exceeds **40.10**, above the target **40**. This rules out this weight family, not dimension-21 polytopes.

For `d≥2` entire A4 components, two marked roots in different components receive `z/(4(d−1))`, where `z` depends on their unordered rooted-orbit pair. Unlisted pair types have zero bonus. There are 2 occurring pair orbits. Summing over all unordered pairs of components and all root choices gives `d/8` times the sum of `z` over the sixteen ordered root choices in one A4 diagram. Each A4 component therefore contributes one eighth of that ordered-pair sum, preserving additive finite-forest separation. The earlier separated-root rule has `z=u_r+u_s` on listed pairs; discarding unlisted pair payments only lowers its vertex cost while preserving these input face payments.

The nonnegative baseline and normalized A1/I2(5)/H3/A4-to-other crosses remain unchanged. Contexts are exact A1 count, I2(5) and A4 counts capped at two, and H3 presence. Actual counts are used for all denominators and rank checks. The complete common four-count enumeration has 102,917 cases from 1,979 face rows, producing 51,443 inequalities and 89,973 variables. Complete finite-queue singular obligations and inherited parameter-marked/dotted obligations are retained. The exact dual uses 2,620 reconstructed constraints; no optimality claim is made.

Generation and saved replay passed. The source, JSON and report are frozen. No separate independent arithmetic review or improved dimension bound is asserted.

```sh
python3 scripts/fresh_sqrt5_dimension21_a4_pairs_20261016.py --hub-cap 2 --verify
```
