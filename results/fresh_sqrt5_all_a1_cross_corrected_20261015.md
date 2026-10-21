# Corrected rank-27 normalized-cross certificate

The corrected [producer](../scripts/fresh_sqrt5_all_a1_cross_corrected_20261015.py)
and [certificate](fresh_sqrt5_all_a1_cross_corrected_20261015.json) give the
same exact maximum **67009/1000 = 67.009**, below **143/2** in dimension 27.
They supersede the vertex-accounting calculation in the
[original frozen diagnostic](fresh_sqrt5_all_a1_cross_20261015.md).

The separate implementation review found that cross-orbit construction
called a rank-two component I2(3) or I2(4), while the original component
summation called it A2 or B2. The latter therefore omitted those rooted
bonuses. This source consistently uses the intrinsic component family in
both places. The original file remains unchanged as historical evidence;
its full table of context maxima must not be used as a theorem dependency.

In the existing rational rule, the affected maximum is at exact A1 count
20: it is **12359/250 = 49.436**, attained by 20A1+2A2+A3, rather than the
original 6687/200. Other rooted A2 bonuses are also included explicitly;
all rooted B2 bonuses in this rule are zero. The overall maximum remains
67009/1000. Re-running the corrected optimization produces the same feasible
rule and checks every context, rather than patching the final number.

The weighting definition, 1,343 parameter-marked graphs, 14,128 face cases,
6,696 variables and upward 1/1000 rounding are as specified in the preceding
report. The corrected calculation includes all 27 feasible exact A1 counts
from zero through 27, excluding 26 because one non-A1 dimension cannot be
filled. This is a feasible-weight certificate, not an optimality claim.

Generation and saved replay pass. The separate reviewer is recorded in its
own package. Complete finite attachment coverage remains a separate
requirement before any dimension theorem can follow.

```sh
python3 scripts/fresh_sqrt5_all_a1_cross_corrected_20261015.py --verify
```
