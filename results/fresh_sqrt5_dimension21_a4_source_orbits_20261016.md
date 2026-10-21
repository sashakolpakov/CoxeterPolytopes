# Dimension-21 A4 source-root refinement: exact obstruction

The [producer](../scripts/fresh_sqrt5_dimension21_a4_source_orbits_20261016.py) and [rational certificate](fresh_sqrt5_dimension21_a4_source_orbits_20261016_cap2.json) separate the two source-root orbits in A4-to-component crosses. The exact dual lower bound is

\[
\frac{121786916150321}{3043000000000}
\approx 40.02199019>40.
\]

Thus this enlarged weight family still cannot meet the dimension-21 threshold. This is a restricted-model obstruction, not a nonexistence theorem or an improved dimension bound. The numerical value is not claimed to be optimal or directly comparable to the earlier nonoptimal certificates.

The model starts from the smaller A4 pair-orbit package. Its contexts remain exact A1 count, I2(5) and A4 counts capped at two, and H3 presence. The nonnegative baseline and every other cross family are unchanged. Every common external A1/I2(5)/H3/A4 count and every finite rank-21 forest is included. The input retains all dimension-feasible singular obligations from the completed finite queue and the inherited parameter-marked/dotted obligations.

For an A4 source component and another destination component outside A1, I2(5), H3 and A4, the new bonus is

\[
\frac{v_{\text{source root orbit},\,\text{destination root orbit}}}{4d},
\]

where `d` is the actual number of entire A4 components. A4 has two source-root orbits, its two endpoints and its two inner roots. Each orbit therefore supplies `2d` source roots, and summing its cross payments contributes exactly `v/2` per destination root. The finite-forest coefficients implement this factor `1/2` explicitly. Capping the bonus context never changes the actual denominator.

The previous 75 destination-root orbits split into 106 occurring source/destination combinations: 44 destinations occur with one source orbit and 31 with both. Missing combinations receive zero bonus. The literal `a4_cross_orbits` table includes both rooted diagrams, canonical source and destination codes, and source multiplicity two; `a4_source_root_orbits` separately records the two source classes. The existing two A4 self-pair orbits and their `z/(4(d−1))` normalization remain unchanged.

There are **1,979 face rows**, **102,917 common-spectator cases**, **51,443 distinct face inequalities**, and **93,045 variables**. The exact dual uses **2,804 constraints**. Generation and saved replay reconstruct the face and forest coefficients and verify rational weak duality. A separate literal weak-duality calculation also passed on 6,121 used columns. No independent implementation review or optimality claim is made.

Generation and saved replay passed. The source, JSON and report are frozen. No further family expansion was run.

```sh
python3 scripts/fresh_sqrt5_dimension21_a4_source_orbits_20261016.py --hub-cap 2 --verify
```
