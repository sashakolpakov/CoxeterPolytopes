# Silver Lambda22 strict hyperideal frontier through H<=150

For the exact arithmetic host

`<-(sqrt(2)-1)> + Lambda22`

in hyperbolic dimension 22, the fixed ray `z=(1,e1/3)` remains strictly
hyperideal after every possible reflective root of real Vinberg height at
most 150.  Its norm is `13/9-sqrt(2)>0`.  Thus this **bounded prefix** still
has infinite volume; this is not a global nonreflectivity statement.

## Complete necessary-superset scan

The exact scan has 58 control levels.  After spatial
minimum and pairing-gate filters, 72
binary Gram types at 18 levels remain.  Exact
PSD extension by `d=e1`, with all gate divisibilities imposed, leaves only
seven post-H14 types whose relaxed wall pairing is nonnegative.

They are excluded as follows:

1. `H=21+15sqrt(2)`: `u-v` is a nonzero S3 vector of norm 16, below 24.
2. `H=36+51sqrt(2)/2`: for all 17 types, either `u-v in S3` has norm at
   most 22 or `u-2v in S6` has norm at most 90.
3. `H=41+29sqrt(2)`: tangency needs an S2 norm-12 vector pairing 6 with
   `e1`; the complete shell has absolute pairings only `0,2,4`.
4. `H=42+30sqrt(2)`: the dangerous S6 difference has norm 12; in fact the
   whole q=12 packet has positive difference norm between 6 and 48.
5. `H=56+40sqrt(2)`, type `(124,52,80)`: `u-2v` would be an S2 norm-12
   vector pairing -6 with `e1`, excluded by the same shell census.
6. `H=56+40sqrt(2)`, type `(100,64,80)`: `u-v in S2` has norm 4, below 12.
7. `H=70+99sqrt(2)/2`: `u in S2`, but every possible `u^2` is 2 or 6
   modulo 8, whereas the four S2 radical classes have residues `0,4,4,4`.

An apparent H72 tangent disappears once the cross gate is imposed correctly:
`u in S6` and `v in S3` force `6 | u.v`, whereas every necessary type has
`u.v=525`.  As a redundant check, all 1,782 signed minimal S6 vectors have
absolute pairing `0,6,12` with `e1`, never the formerly required 18.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h150_frontier.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h150_frontier.py
```

Source: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/LAMBDA22.html
