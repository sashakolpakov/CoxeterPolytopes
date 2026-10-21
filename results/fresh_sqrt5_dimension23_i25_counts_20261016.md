# Dimension-23 exact I2(5)-count contexts

The [producer](../scripts/fresh_sqrt5_dimension23_i25_counts_20261016.py) and [rational certificate](fresh_sqrt5_dimension23_i25_counts_20261016_cap2.json) give an exact feasible weight rule with complete finite-forest maximum **485211/10000 = 48.5211**, strictly below the dimension-23 contradiction threshold **99/2 = 49.5**. This is a weight-model certificate; the global arithmetic conclusion also requires independent implementation review and complete face-coverage composition.

The preceding presence-only model had an exact lower bound above 49.63. This enlargement distinguishes an entire I2(5) count of zero, one, or at least two. Every nonnegative connected bonus, A1-to-component bonus and I2(5)-to-component bonus depends on the exact number `a` of entire A1 components and on `min(b,2)`, where `b` is the exact number of entire I2(5) components. The fixed nonnegative q2/H4 baseline is unchanged. Normalized cross-pair weights are `v/a` and `u/(2b)` respectively; the latter destinations exclude A1 and I2(5), so the two cross families do not overlap. All actual common spectator counts `(a,b)` are included at every quadrilateral.

The complete input includes the dimension-feasible singular stopping obligations from the completed 5,172-state finite queue, in addition to the earlier parameter-marked catalogue. There are **2,332 face rows**, **54,993 common-spectator cases**, **29,666 distinct inequalities**, and **16,474 variables**. Literal orbit tables contain 243 connected-pair orbits, 187 A1 cross-root orbits, and 102 I2(5) cross-root orbits. The saved variable order explicitly records `orbit`, `A1_count` and `I25_context`.

For every exact `(a,b)`, summing normalized cross bonuses over the source components makes the vertex weight additive over the other components. The exact finite-component dynamic program fixes both counts, excludes A1 and I2(5) from the residual component choices, and includes every remaining finite Coxeter component through rank 23. It checks **144 feasible exact-count contexts**. The largest weight occurs at `(a,b)=(9,2)` on `9A1 + A4 + A6 + 2I2(5)`.

Numerical optimization proposed the rule. Uniform rational repair and upward rounding to a denominator of 10,000 produced the saved feasible certificate. Every face inequality and every finite-forest maximum is then checked over the rational numbers. There is no optimality claim and no use of floating point in saved-certificate verification.

Generation, metadata refresh, and saved replay passed:

```sh
python3 scripts/fresh_sqrt5_dimension23_i25_counts_20261016.py --hub-cap 2 --verify
```

The source and JSON are frozen for independent review. This report does not itself assert a new global dimension bound or external mathematical validation. No prior frozen artifact was edited.
