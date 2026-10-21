# Silver Lambda22 replacement survives through real height 10000

For `<-(sqrt(2)-1)> + Lambda22` in `H22`, the fourth accepted cutter at
`2786+1970sqrt(2)` leaves the exact ray

`z=(1,(28/87)d)`,  `d^2=4`,

with norm `10705/7569-sqrt(2)>0`.  Starting strictly after that cutter, the
complete control scan through the rational identity-height bound `H=10000`
contains 129 controls and 1,555 gate-compatible necessary Gram types.

Only one relaxed type can pair nonnegatively with `z`.  It occurs at
`3196+2260sqrt(2)`, has coefficient Gram `(3194,1600,2260)`, and maximizing
pairings `(d.u,d.v)=(113,80)`.  An exhaustive integral PSD enumeration,
followed by the explicit `[-3,3]^3` short-combination gate, checks all 665
extensions not already contradicted by `min(Lambda22)=4` and
finds exactly one nonnegative extension, namely that maximizer.  Equality
in the `v` and cross Schur bounds
forces `v=20d`; then `e=u-28d` has norm two, contradicting
`min(Lambda22)=4`.  Thus the type is empty.

The complete bounded prefix through real height 10000 therefore retains a
strict hyperideal patch and is not finite-volume.  This is not Vinberg
termination and does not prove global nonreflectivity.

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h10000_replacement.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h10000_replacement.py
```
