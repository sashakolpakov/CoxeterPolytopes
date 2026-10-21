# Exact global budget 40 on the current face domain

The [source](../scripts/fresh_sqrt5_global_discount_20261016.py) and
[exact rational rule](fresh_sqrt5_global_discount_20261016.json) achieve
**global vertex budget 40 in dimension 21**. This is a feasible weighting
proposal, not a dimension theorem: complete coverage of newly affected
baseline-one quadrilateral extensions is still missing.

This rule counts all unordered root pairs in a full rank-21 vertex forest.
There is no distinguished facet root. It retains the adjacent-edge discount,
the two-edge E6 exemption, actual H3 counts in all six ordinary cross families,
and E8 presence in H4 compensation. The distance-two discount is exactly zero.
The adjacent discount is an exact rational number satisfying
229/1000 < δ < 23/100, approximately 0.229636658054. Every bonus is nonnegative.

The former marked-forest cuts supplied 4,435 distinct complete forests, each
recounted literally and compared with additive component totals. Six fixed-40
solves added constraints from the complete forest oracle, ending with 4,762
literal cuts. Rational basis recovery solved 8,197 equations using 7,797
sparse pivots and a 400-variable FLINT remainder.

The final exact maximum is **40 across all 3,261 feasible count/presence
contexts** in the finite-forest DP. Every returned maximizing forest was
separately recounted as actual unordered root pairs. There are 380 chosen
context maximizers at 40; this is not a list of every tied forest. The rule
has 32,554 allocated coefficients and 7,160 nonzero coefficients. All 231
inherited coefficients requiring a facet mark are zero.

All 996 states of the complete triangle catalogue and all 37 known paid-quad
stopping cases were rebuilt without facet marks. Common external data include
the exact E8 count, shared by all corners. Together with 1,979 inherited
baseline-zero quad rows, this gives 62,817 distinct inequalities: 352,022
triangle contexts, 20,501 paid-quad contexts, and 227,408 inherited-quad
contexts. Every inequality is checked over the rationals. Generation and
fresh-process saved replay both pass.

Since δ < 1/2 and each discounted corner has baseline one, any quadrilateral
with baseline at least two still pays at least one. The unfinished coverage
therefore consists of baseline-one quadrilaterals whose paid corner can be
discounted. The frozen [frontier](fresh_sqrt5_adjacent_discount_frontier_20261016.md)
contains 745 potentially relevant minimal classes and 3,820 relevant first
extensions before isomorphism reduction. Its old arithmetic exclusions and
sealed stopping cases are useful input, but the complete finite and dotted
attachment closures must still be established. The E6 exemption alone is
not hereditary: an E6 component can grow into a discounted E7/E8 component.

A separate narrow source read checked restoration of full forests, unmarked
component totals, all unordered root pairs, E8 state handling, and the absence
of facet-mark subtraction. It was not an independent solver replay or a
complete quadrilateral-coverage review. The existing compact arithmetic
bound-21 proof candidate is unchanged; a compact bound of 20 is not asserted.

```sh
python scripts/fresh_sqrt5_global_discount_20261016.py --verify
```
