# Independent review of weights with capped I2(5) counts

The dimension-23 rule has independently verified maximum vertex weight
**485,211/10,000 = 48.5211**, below the required threshold **99/2 = 49.5**.
Generation and a fresh saved replay pass. The
[reviewer](../scripts/fresh_sqrt5_hub_count_weights_20261016_review.py) and
[certificate](fresh_sqrt5_hub_count_weights_20261016_review_n23_cap2.json)
use the frozen separate literal-graph engine, without importing a weight
producer or optimizer. A separate source review also checks the changed
normalization and spectator accounting.

## Complete weight-side checks

The review reconstructs all 2,332 dimension-feasible face rows and their
54,993 common spectator cases. These consist of 1,099 rows from the earlier
parameter-marked catalogue and 1,233 singular stops from the completed finite
queue. Every one of the queue's 2,070 singular stops is accounted for before
the dimension filter. Expanded singular states are included. All 2,397
original literal-input maps and their exact dotted parameters are rechecked;
1,839 fit dimension 23.

The rule has 16,474 nonnegative bonus variables: 243 connected-pair orbits,
187 A1-to-root orbits and 102 I2(5)-to-root orbits, with their count contexts.
All literal type identifications, intrinsic rank-two aliases, marked keys,
and variable ranges are checked separately. The minimum face slack is
1/160000. All 74 allowed finite component types and 144 possible exact
count pairs `(a,b)` at rank 23 are covered by the independent vertex DP.
Every saved maximizing forest is checked directly as well.

## Why the extra context is valid

Let `a` and `b` count entire A1 and I2(5) components. Every bonus depends on
`a` and `min(b,2)`, distinguishing zero, one, and at least two I2(5)
components. This cap selects the coefficient only. A face still divides
the A1 cross payment by the actual `a` and the I2(5) cross payment by the
actual `2b`. The same external counts are added at all four corners before
evaluating these quantities.

Each eligible destination root has exactly `a` A1 source roots and exactly
`2b` I2(5) source roots. The latter family excludes destinations in A1 and
I2(5), so the two families never overlap. Both denominator cancellations
are checked with rational arithmetic. Fixing exact `(a,b)` leaves additive
component costs; a reverse-type unbounded rank DP enumerates all residual
forests, excluding the two explicitly counted source-component types.
All bonuses are nonnegative, and the short-path/H4 baseline is unchanged.

This is a weight review. Complete geometric branch coverage is supplied by
the [separate binder](fresh_sqrt5_hub_count_coverage_20261016.md). Numerical
optimization merely proposed the coefficients. No optimality, external
validation or literature-priority claim is made.

```sh
python3 scripts/fresh_sqrt5_hub_count_weights_20261016_review.py --dimension 23 --hub-cap 2
python3 scripts/fresh_sqrt5_hub_count_weights_20261016_review.py --dimension 23 --hub-cap 2 --verify
```
