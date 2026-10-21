# Exact obstruction for the active/inactive compact-rank-21 support successor

The selected support successor has the exact lower bound

\[
z\geq \frac{40004479046115}{1000000000247}
 =40.004479036233896\ldots>40.
\]

Its excess is 4479036235/1000000000247. Thus the newly proved active/inactive
D5+2E8 charge rules, inserted into the old saved weight family, do not yet give
the target value 40. This is an obstruction to this particular charge model;
it is not an existence statement about polytopes or a new dimension bound.
It does not obstruct further face-dependent rules or a different weight family.

## Exact model and weak-duality proof

The certificate uses the 2,804 positive-dual-support rows in
`fresh_sqrt5_dimension21_a4_source_orbits_20261016_cap2.json`, with the old
D5+2E8 row replaced by its active version and one inactive row added. Every
other selected face and forest inequality is unchanged. There are 2,805 rows,
6,121 occurring columns and 17,387 nonzero coefficients. Variables are
nonnegative. The objective column is 93045 and all other objective coefficients
are zero. The geometric validity and precise scope of the new two-row rule
are established in the separately reviewed active/inactive geometry packet.

Write these literal rational inequalities as Ax <= b, with x >= 0 and objective
c.x = z. The certificate contains 2,716 negative row multipliers y. The replay
checks y <= 0 and every column of y^T A <= c^T, and sums y^T b exactly. Hence

\[
y^T b\leq y^T A x\leq c^T x=z.
\]

All columns are checked, including zero-primal columns; no floating-point
feasibility tolerance is used. Since these rows are a subset of the full
successor constraints, this exact lower certificate also excludes z <= 40 in
that full successor. A full primal-separation run is unnecessary for this
negative conclusion.

## Recovery without another optimizer

The preceding guarded packet retained the numerical proposal. This script
imports no numerical package and runs no optimizer. It rounds positive face
masses down and forest masses up on a grid of 10^12. Each forest mass supplies
nonnegative capacity in bonus columns; each face mass creates nonnegative
demand. For each overloaded column it decreases an incident face mass until
the excess is gone. Decreasing a face mass decreases demand in every column,
so this procedure cannot invalidate an already repaired column. It then
normalizes the total forest mass to one and independently replays the exact
column inequalities and objective sum.

Only one column was initially overloaded and only one face mass was clipped.
The objective loss before normalization was 15147/125000000000, and the
normalizer was 1000000000247/1000000000000. The compact JSON retains all exact
nonzero row multipliers and all nonzero forest masses. Source and input hashes
are checked by `--verify`.

## Which geometry matters next

The certificate is a mixture of forests; no single forest alone has been
proved to force the excess. Its largest positive forest masses include:

| Forest | Exact positive mass |
|---|---:|
| A1 + 2E8 + H4 | 37390344603/1000000000247 |
| A1 + D6 + 2E7 | 33113153373/1000000000247 |
| A1 + D7 + E8 + H3 + I2(5) | 28761456056/1000000000247 |
| A1 + E7 + E8 + H3 + I2(5) | 28761456056/1000000000247 |
| D5 + 2E8, active state | 26820895133/1000000000247 |
| A1 + 2D7 + E6 | 26820895133/1000000000247 |

The A1+D5+E7+E8 forest has **zero** mass in this witness. A new saving at only
that forest would therefore leave this exact obstruction unchanged. The first
row, A1+2E8+H4, is a better immediate target for a complementary incidence rule.
A saving in a positive-mass row lowers this particular dual objective by its
mass times the saving; this sensitivity statement does not establish that the
resulting successor would be feasible.

## Reproducibility and verification scope

```sh
python3 -B scripts/fresh_compact21_priority_support_dual_repaired_20261018.py --verify
```

Generation took 0.292 seconds with 38,404,096 bytes peak RSS; the fresh replay
also passed. Initial free disk was 197.85 GiB and available RAM was 9.32 GiB.
Both runs held the shared lock, used a 59-second alarm, streamed the large
source model, and left no temporary files. No Lean proof, external peer review,
or dimension improvement is claimed. The result is a fully rational,
internally replayed obstruction for the stated support successor.
