# Independent review of four H20 cusp-edge obstructions

**PASS.** A separate checker verifies all 418 finite contact profiles and all fifteen unrestricted `I₂(m)` tails, for every integer `m ≥ 7`. On these four specified supports, an acute cap giving a finite Coxeter endpoint cannot cut the saved strict escape. No new finite-volume example follows.

The [review source](../scripts/fresh_h20_e7_repaired_cusp_edgewalk_20261014_review.py) imports the older sparse squarefree-radical arithmetic, not the producer's tower arithmetic or search functions. Its [certificate](fresh_h20_e7_repaired_cusp_edgewalk_20261014_review.json) binds the frozen inputs and the [primary exact certificate](fresh_h20_e7_repaired_cusp_edgewalk_20261014_exact.json).

For each support, the review reconstructs the nineteen physical dual roots directly from the ambient form, checks positive definiteness, and rebuilds the orthogonal projection of the center. It verifies the strict signs at the opposite null ray for all 76 non-support walls and for the escape.

Profile completeness has an additional direct check independent of finite-diagram classification. The normalized inverse Gram matrices have nonnegative entries and diagonal entries at least `1/2`. Two contacts in one component already cost at least 2; four contacts do likewise. Every nonisolated root has inverse diagonal at least `2/3`, so a label at least 7 must be the sole contact at an isolated root. Three contacts must all have label 3. These exact inequalities reproduce the complete 418-profile census and fifteen symbolic tails.

For every finite profile, the review reconstructs the residual norm and both separating pairings at the saved rational parameter. Eight blocking inequalities involve the bridge radical; they pass in the independent arithmetic. For every symbolic tail it reconstructs both quadratic polynomials on `√3 ≤ 2cos(π/m) ≤ 2`, checks convexity and strict wall negativity, and checks strict positivity at the escape polynomial's minimum. Eleven minima are at the left endpoint and four at the right.

The cap orientation is forced: `ν = N₉ + N₈₈`, so acuteness gives `r·ν ≤ 0`. Equality is incompatible with the positive residual norm of a finite extension. Thus the parameter is positive. Every relevant pairing is strictly decreasing in that parameter, which turns each separator into an exclusion for all real cap positions.

Generation and saved-result replay passed:

```sh
python3 scripts/fresh_h20_e7_repaired_cusp_edgewalk_20261014_review.py --verify
```

This review covers only the four supports listed in the [primary report](fresh_h20_e7_repaired_cusp_edgewalk_20261014.md). It does not exclude non-cutting additions, other supports, ideal endpoint extensions, or changes to the wall set. The chamber still has 95 verified facets, at least 116 verified finite vertices, three specified bounded cusps, and an explicit strict spacelike escape proving infinite volume.
