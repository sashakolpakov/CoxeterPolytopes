# Dimension-21 normalized A4 crosses: restricted obstruction

The [producer](../scripts/fresh_sqrt5_dimension21_a4_cross_20261016.py) and [rational certificate](fresh_sqrt5_dimension21_a4_cross_20261016_cap2.json) add A4-to-component and A4-to-A4 bonuses to the dimension-21 H3/I2(5)/A1 model. An exact dual lower bound exceeds **40.08**, above the required threshold **40**. This is a restricted-model obstruction, not a dimension theorem.

The A4 count `d` is capped at two only in bonus contexts. Every face normalization and rank calculation uses the actual count. The new families are:

- Between an A4 root and a rooted destination component outside A1, I2(5), H3 and A4, assign `v/(4d)`. Summing all `4d` A4 source roots gives one `v` per destination root.
- For marked roots `r,s` in two different A4 components, assign `(u_r+u_s)/(4(d−1))`, for `d≥2`. Each root has `4(d−1)` partners in other A4 components. Summing over unordered root pairs therefore gives exactly one `u` per root.

These families are disjoint from the existing normalized A1, I2(5), and H3 crosses. The nonnegative q2/H4 baseline is unchanged. The new literal tables have 75 A4 destination-root orbits and 2 A4 self-root orbits. They were selected after unpaid A4 pairs, especially A4–A4, appeared prominently in the preceding dual support.

Every common external count of A1, I2(5), H3 and A4 is enumerated. Vertex separation fixes all four exact counts and excludes those component types from its residual finite-forest dynamic program. The input includes the completed finite-queue singular obligations and inherited parameter-marked/dotted obligations. There are 1,979 face rows, 102,917 common-spectator cases, 51,443 distinct face inequalities and 89,973 variables. Exact weak duality is verified on 2,187 reconstructed constraints. No optimality claim is made.

Generation and saved replay passed. The source, JSON and report are frozen. No independent arithmetic review or stronger dimension bound is asserted.

```sh
python3 scripts/fresh_sqrt5_dimension21_a4_cross_20261016.py --hub-cap 2 --verify
```
