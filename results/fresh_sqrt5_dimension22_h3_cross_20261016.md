# Dimension-22 feasible weights with normalized H3 crosses

The [producer](../scripts/fresh_sqrt5_dimension22_h3_cross_20261016.py) and [rational certificate](fresh_sqrt5_dimension22_h3_cross_20261016_cap2.json) give an exact feasible weight rule with complete finite-forest maximum **21977/500 = 43.954**, below the dimension-22 contradiction threshold **44**. This is a weight-model certificate. An arithmetic dimension conclusion still requires independent implementation review and complete face-coverage composition.

Every nonnegative bonus depends on the exact number `a` of entire A1 components, `min(b,2)` for the number `b` of entire I2(5) components, and the presence of an entire H3 component. The fixed nonnegative q2/H4 baseline is unchanged. In addition to the existing normalized A1 and I2(5) crosses, a new H3-to-component pair receives `w/(3c)`, where `c` is the exact number of entire H3 components. Its destination excludes A1, I2(5) and H3. The three cross families are therefore disjoint. Summing the new weights over all `3c` source roots contributes `w` once per destination root, preserving additive component optimization for fixed exact `(a,b,c)`.

Every common spectator count `(a,b,c)` is included. Residual rank zero or at least two is possible using components other than A1, I2(5) and H3; rank one is excluded. The exact finite-forest dynamic program fixes all three counts, excludes the three distinguished component types from residual choices, and then includes every remaining allowed finite Coxeter component through rank 22. Bonus contexts use capped counts, but the normalization denominators and all vertex/face count tests use actual counts.

The complete model has **2,151 face rows**, **89,255 common-spectator cases**, **40,294 distinct inequalities**, and **30,169 variables**. Literal tables contain 232 connected-pair orbits, 178 A1 cross-root orbits, 97 I2(5) cross-root orbits, and six H3 cross-root orbits. The last six have rooted destinations of types A3, A4, D4 and D5. The saved variable order explicitly records `orbit`, `A1_count`, `I25_context` and `H3_context`. The H3 table is named `h3_cross_orbits`.

The full finite-forest verification checks **410 feasible exact-count contexts**. Its maximum is attained, among others, by `A1 + D4 + 2D6 + H3 + I2(5)` and by `3A1 + A5 + 2A6 + I2(5)`. The input includes every dimension-feasible singular stopping obligation from the completed finite queue in addition to the earlier parameter-marked catalogue. A global coverage binder must additionally recheck all older closure branches directly with this rule.

Numerical optimization proposed the rule. Uniform rational repair and upward rounding of all bonuses to denominator 10,000 yielded the saved certificate, with 8,570 nonzero bonuses. Every face inequality and the complete vertex-forest maximum are then checked over the rationals. Saved replay performs no numerical optimization. There is no optimality claim.

Generation and saved replay passed. The source and JSON are frozen for independent review:

```sh
python3 scripts/fresh_sqrt5_dimension22_h3_cross_20261016.py --hub-cap 2 --verify
```

No external mathematical validation, novelty, or global dimension theorem is asserted by this producer package. No prior frozen artifact was edited. The prepared D7-context scout was not run because this H3 rule reached the target.
