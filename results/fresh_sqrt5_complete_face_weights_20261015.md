# Exact weights including every finite-queue stopping state

The [producer](../scripts/fresh_sqrt5_complete_face_weights_20261015.py) now
incorporates every dimension-feasible singular state from the completed
5,172-state finite attachment queue, in addition to the earlier conservative
catalogue and its dotted obligations. The reoptimized rules satisfy:

| Dimension | Exact maximum vertex weight | Strict threshold | Finite-queue rows | All face cases |
|---|---:|---:|---:|---:|
| 25 | 59 | 60 | 1,488 | 21,897 |
| 26 | 31,017/500 = 62.034 | 65 | 1,623 | 24,609 |
| 27 | 67,009/1,000 = 67.009 | 143/2 = 71.5 | 1,767 | 27,523 |

The exact certificates are separate for
[dimension 25](fresh_sqrt5_complete_face_weights_20261015_n25.json),
[dimension 26](fresh_sqrt5_complete_face_weights_20261015_n26.json), and
[dimension 27](fresh_sqrt5_complete_face_weights_20261015_n27.json).
Generation and saved replay pass for all three. Independent implementation
review and the full mathematical coverage composition remain separate
requirements before asserting a dimension theorem.

## What the new inputs repair

The earlier partial rules missed real constraints in the finite queue:
173 dimension-25 cases and twelve cases each in dimensions 26 and 27 were
unpaid in live checks. Every singular terminal or intermediate stopping
choice in the completed queue is now added before optimization, including
all its dimension-feasible exact spectator-A1 counts. The input contains
2,709, 2,905, and 3,110 face rows respectively; those row counts include
duplicates of earlier catalogue graphs and are not counts of distinct
isomorphism classes. Duplicate linear inequalities are coalesced exactly.

The full finite queue contains 2,070 singular states. Rows whose corner
rank exceeds the target dimension are excluded from that target's weight
calculation. All remaining rows are included, regardless of whether they
are sealed, expanded, or already paid by an earlier rule. The completed
queue JSON is a direct hashed dependency. Its independent arithmetic and
completeness review is required separately by the coverage binder.

## Rule and verification

The baseline, exact-A1-count connected bonuses, and normalized A1-to-whole-
component cross bonuses are those of the
[corrected rank-27 construction](fresh_sqrt5_all_a1_cross_corrected_20261015.md).
The intrinsic family name is used consistently, including A2=I2(3) and
B2=I2(4). A cross angle gets u_a/a when there are exactly a entire A1
components; summing all such angles at a vertex contributes u_a per rooted
node of the other component. A1–A1 pairs remain zero.

Every face inequality and every vertex maximum is checked in rational
arithmetic. The latter enumerates all allowed finite forests by exact
unbounded rank dynamic programming for each feasible exact A1 count.
Upward rounding to a 1/1000 bonus grid is followed by another full exact
verification. Floating-point optimization supplies proposals only; no
optimality claim is made. Nonnegative bonuses preserve all reviewed
triangle and persistent short/H4 payments.

This calculation supplies numerical face and vertex inequalities for its
full stated input set. The separate coverage argument must still connect
every actual compact arithmetic quadrilateral to the enumerated branches,
check the older sealed/closure branches under these new rules, and combine
the strict h-vector inequality with the existing dimension-28/29 rule.

```sh
python3 scripts/fresh_sqrt5_complete_face_weights_20261015.py --dimension 25 --verify
python3 scripts/fresh_sqrt5_complete_face_weights_20261015.py --dimension 26 --verify
python3 scripts/fresh_sqrt5_complete_face_weights_20261015.py --dimension 27 --verify
```
