# Dimension-22 obstruction for exact A1/I2(5) count contexts

The [producer](../scripts/fresh_sqrt5_dimension22_i25_counts_20261016.py) proves that two enlargements of the existing nonnegative weight family still do not reach the dimension-22 threshold **44**. These are restricted-model obstructions on the conservative complete face catalogue, not dimension theorems or existence results.

| I2(5) context | Variables | Distinct inequalities | Exact certified lower bound |
|---|---:|---:|---:|
| Zero, one, at least two | 14,942 | 25,074 | greater than 44.11 |
| Every exact count | 62,849 | 29,258 | greater than 44.03 |

Both models include 2,151 face rows and 45,928 common spectator-count cases. The original nonnegative q2/H4 baseline is fixed. Connected bonuses, normalized A1 crosses and normalized I2(5) crosses depend on the exact entire-A1 count and on the stated I2(5) context. The latter hub cross has denominator `2b`; its destinations exclude A1 and I2(5). Thus the full exact-b model already distinguishes every possible pair of hub counts `(a,b)`.

The [capped certificate](fresh_sqrt5_dimension22_i25_counts_20261016_cap2.json) contains 566 terms, including 110 exact finite-forest inequalities. The [exact-count certificate](fresh_sqrt5_dimension22_i25_counts_20261016_exact.json) contains 819 terms, including 193 forest inequalities. Their rational weak-duality values are approximately 44.1164146868 and 44.0354441338 respectively. The full rational values are saved, and verification checks their comparison with 44 exactly.

The numerical LP only proposes dual multipliers. Recovery drops positive/tiny numerical noise, removes face terms with unsupported columns, normalizes the forest mass and shrinks face multipliers to meet every rational column constraint. Every selected face inequality is rebuilt from the stated catalogue and spectator counts; every forest inequality is rebuilt from its literal finite-component multiset. Saved replay runs no optimizer. There is no optimality claim or independent implementation-review claim.

Generation and saved replays passed:

```sh
python3 scripts/fresh_sqrt5_dimension22_i25_counts_20261016.py --hub-cap 2 --verify
python3 scripts/fresh_sqrt5_dimension22_i25_counts_20261016.py --hub-cap 0 --verify
```

The source and both JSON files are frozen. The next concrete enlargement is a third entire-component context, such as E8 or A2, while retaining exact common spectator enumeration and additive vertex optimization. The dual supports include substantial coefficients on forests `A1+D5+2E8`, `A1+3E7`, and `E6+2E8`, motivating an E8 split.
