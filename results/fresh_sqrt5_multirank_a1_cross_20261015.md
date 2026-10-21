# Exact normalized-cross rules for dimensions 25 and 26

The [rank-specific producer](../scripts/fresh_sqrt5_multirank_a1_cross_20261015.py)
finds exact feasible rules at two further dimensions:

| Dimension | Exact vertex maximum | Strict face-count threshold | Face cases |
|---|---:|---:|---:|
| 25 | 59 | 60 | 11,625 |
| 26 | 60,799/1,000 | 65 | 12,846 |

The certificates for [dimension 25](fresh_sqrt5_multirank_a1_cross_20261015_n25.json)
and [dimension 26](fresh_sqrt5_multirank_a1_cross_20261015_n26.json) use the
corrected intrinsic component-family names described in the
[rank-27 erratum and corrected calculation](fresh_sqrt5_all_a1_cross_corrected_20261015.md).
Unfrozen earlier development outputs were replaced by complete corrected
runs. Generation and saved replay pass for both certificates.

The rule is the same class as the corrected rank-27 construction: retain the
distance-two/H4 baseline, let connected bonuses depend on the exact number
a of entire A1 components, and give each A1-to-component rooted pair bonus
u_a/a. The total cross contribution per non-A1 root is u_a. The coefficients
are optimized separately for each dimension; these are not claims that one
rule works unchanged in every dimension.

Every dimension-feasible spectator count in the conservative frozen face
catalogue is checked exactly. All vertex forests are separated by exact
rank dynamic programming at each exact A1 count. Floating-point proposals
are repaired with rational arithmetic, bonuses are rounded upward to a
1/1000 grid, and every face and vertex condition is checked again. There is
no optimality claim.

The same model did not produce certificates for dimensions 24 or 23. The
numerical relaxations reached 55.5 against 54, and 52 against 99/2,
respectively. These observations are numerical model obstructions only;
no exact infeasibility certificate or nonexistence conclusion follows.

All finite-queue zero-touch states still need direct evaluation under each
new rule, and separate implementation reviews are required. These partial
weight certificates alone establish no dimension bound.

```sh
python3 scripts/fresh_sqrt5_multirank_a1_cross_20261015.py --dimension 25 --verify
python3 scripts/fresh_sqrt5_multirank_a1_cross_20261015.py --dimension 26 --verify
```
