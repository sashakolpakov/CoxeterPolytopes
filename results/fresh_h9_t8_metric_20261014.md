# Parametric metric tests for T-shaped case 8

Two continuous metric families of the saved T-shaped incidence are excluded. Neither result excludes the whole incidence. Both preserve genuine unbounded finite-label categories, rather than substituting the angle π/7.

The input is [case 8](fresh_h9_tree_cube_wedges_20261014_case_08.json) of the [tree catalogue](fresh_h9_tree_cube_wedges_20261014.json), with 462 actual vertices. The support `[1,2,4,5,7,8,12,13]` occurs in an actual vertex. Its first pair has label m≥7, the middle two pairs have label 4, and its last pair initially has label 4. Twelve high-label three-path tests and twelve relaxed affine `(4,3,4)` tests force every cross-block pair to be right. The resulting rank-eight positive subspace leaves at most two residual dimensions in H9.

The saved assignment has leaves 11 and 17 with own contacts `(3,3)` on the last block. They are mutually right. Literal eight-wall Schur budgets force both right to the first I2(m) block and permit at most one foreign π/3 contact to one of four other B2 roots. Wall 14 has a fixed projection: its pairings to roots 2,4,5 are −1/2 and its other support pairings vanish. Two actual seven-wall budgets force those vanishing contacts; two exact finite-vertex failures at label 4 restrict its contacts to the anchors to labels 2 or 3.

Thus there are 25 anchor-pattern pairs and 100 residual cases. Write t=cos²(π/m). The only parameter contribution is the third residual square `−√2/2−1/[4(1−t)]`. Each determinant is `D0−A/[4(1−t)]`, where the anchor determinant A is negative. It is increasing for 3/4<t<1. Every exact boundary value at t=3/4 is positive, with minimum `3/4+5√2/8`. This proves rank at least eleven for every finite integer m≥7, irrespective of the unassigned dotted distances.

The [first certificate](fresh_h9_t8_metric_20261014.json) has a [separate Fraction arithmetic review](fresh_h9_t8_metric_20261014_review.md). Both generation and saved replays pass.

## Two continuously varying internal labels

After inserting that gate, a [bounded continuation](fresh_h9_t8_metric_20261014_search_case_08.json) finds the next assignment in 857 nodes. It changes only the last internal label from 4 to 5. The [two-label extension](fresh_h9_t8_metric_20261014_two_labels.json) covers every last label n≥4 at once.

The common own-block projection cost is `k=1/[2(1−cos(π/n))]`, at least `1+√2/2`. All domain bounds persist: old Gram failures are monotone, the anchor budgets delete one own root, and wall 14's foreign-contact costs only increase. Each residual determinant is affine in k, with derivative `a1*q14+(ac−bc)²>0`: a1 is the negative slope of the anchor determinant, and q14 is negative. Therefore the first certificate's boundary remains the global minimum for all m≥7,n≥4. A [separate bivariate Fraction review](fresh_h9_t8_metric_20261014_two_labels_review.md) checks every coefficient identity and the guard transfer.

## Changed anchor profiles with two high labels

The [next continuation](fresh_h9_t8_metric_20261014_family_search_case_08.json) reaches another assignment in 932 nodes. Now both internal labels satisfy m,n≥7, and the final anchors each have own contacts `(2,3)`. This lies outside the preceding family. The [new certificate](fresh_h9_t8_metric_20261014_two_high.json) recomputes the changed literal guards, using the affine mutual-contact triangle at root 13 and budgets that retain 13 and delete 12.

Here `k=1/[4(1−cos²(π/n))]`. The exact polynomial `64u³−80u²+24u−1` vanishes at u=cos²(π/7), is increasing above 3/4, and has value −29/125 at 4/5. Consequently k>5/4 for every n≥7. All 100 determinant lower bounds remain positive at t=3/4,k=5/4, with minimum `(1+√2)/8`. The same monotonicity proves the entire two-high-label family impossible. A [separate Fraction review](fresh_h9_t8_metric_20261014_two_high_review.md) independently checks the changed literal domains, all 100 bivariate determinants, and the threshold through the alternative cos(2π/7) cubic. This is an additional exact family exclusion, not a finite-label sweep or a complete search of case 8.

All constructions and replays are under the separate `fresh_h9_t8_metric_20261014` prefix. The star prefix remains frozen. Future continuation should skip both proved T8 families and preserve any remaining category-7 parameters symbolically. No compact H9 example has been obtained.

```sh
python scripts/fresh_h9_t8_metric_20261014.py --verify
python scripts/fresh_h9_t8_metric_20261014_two_labels.py --verify
python scripts/fresh_h9_t8_metric_20261014_two_high.py --verify
python scripts/fresh_h9_t8_metric_20261014_search.py --case 8 --limit 100000 --verify
python scripts/fresh_h9_t8_metric_20261014_family_search.py --case 8 --limit 100000 --verify
```
