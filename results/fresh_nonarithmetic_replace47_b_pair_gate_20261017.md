# Complete label gate for the paired B endpoints

This certificate tests a different topology while retaining the current47: a left cap handles both7–30 and7–47, and a right cap handles both30–52 and47–52. Both caps are rightS6, and their finite mutual contact is intended to create connecting endpoints on30 and47. **No pair in this family works, with all labels covered.**

The left7–30 endpoint has old diagramH4² and forces rightS6,7,30. Positivity at7–47 leaves precisely three contacts to47, with labels2,3,4. For the right cap, the B2 endpoint30–52 and orthogonal pair47–52 give precisely

```
(m30,m47,m52) = (2 or3, m, 2), m>=2;
             or (2, m, 3), m=2,3,4,5.
```

For right47 labels2–5, there are12 right profiles. Together with the three left profiles these yield30 norm branches. Fifteen branches have a strictly positive old-wall product; another falls in a strict consecutive-angle gap. The remaining three left and eleven right branches satisfy every64-wall native gate. Their33 pair products split into nine obtuse products, eight dotted products and sixteen strict consecutive-angle gaps. Thus none supplies a finite mutual Coxeter contact. All thresholds are enclosed exactly, and the gap between consecutive labels excludes every intermediate label without a finite-alphabet assumption.

For the two infinite right47 tails `m>=6`, finite rank3 classification forces the left cap to be right47 and the pair to be right to each other: an I2(m) with m at least6 has no connected spherical extension. There is only one surviving all-right left cap. The full frame `(n30,n47,n52,left)` then determines the right root from products `(-2cos(pi/m30),-c,0,0)`, where `c=2cos(pi/m)>=sqrt(3)` and `m30` is2 or3. Its norm-two equation is a quadratic in c. Both leading coefficients are positive, both linear coefficients are nonnegative, and the two polynomial values at sqrt(3) are strictly positive, respectively greater than0.86357 and2.09719. Both polynomials are increasing thereafter. No tail solution exists.

The preliminary [numerical bank](fresh_nonarithmetic_replace47_b_pair_20261017.json), with [source](../scripts/fresh_nonarithmetic_replace47_b_pair_20261017.py), sampled the isolated-label tails through128 plus192,256,512,1024. It found no finite pair, but that sample was not an all-label proof. The exact [gate source](../scripts/fresh_nonarithmetic_replace47_b_pair_gate_20261017.py) and [certificate](fresh_nonarithmetic_replace47_b_pair_gate_20261017.json) reconstruct all low profiles directly from the frozen geometry and prove the full tail inequalities. Generation and fresh saved replay pass; shared frozen geometry/arithmetic is used, with no separate independent arithmetic review claimed.

The right cap is assumed rightS6 in this family. Nonzero A2 profiles on that cap, different endpoint duties, three or more caps, and retained-wall replacements remain outside. The separate [two-edge gate](fresh_nonarithmetic_replace47_edges1252_3052_gate_20261017.md) does cover all A2 profiles for the distinct12–52/30–52 duty. No new wall is appended and no finite-volume conclusion follows.
