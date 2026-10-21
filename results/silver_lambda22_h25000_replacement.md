# Silver Lambda22 replacement survives through real height 25000

For `<-(sqrt(2)-1)> + Lambda22` in `H22`, inherit the independently audited
strict hyperideal ray

`z=(1,(28/87)d)`, with `d^2=4`.

The complete interval strictly after rational identity height `10000` and
through inclusive height `25000` contains 301 reflective controls, 5,575
gate-compatible binary coefficient-Gram types, and 2,663,705 locally
admissible integral PSD extensions by `d` (counted with multiplicity).

Exactly two types have any nonnegative extension.

1. At `H=6392+4520sqrt(2)`, norm `q=4`, the Gram type
   `(12772,6400,9040)` has 437 extensions and exactly one nonnegative one:
   `(d.u,d.v)=(226,160)`, of determinant zero and margin
   `(64-44sqrt(2))/87>0`.  Cauchy equality forces `v=40d`; then
   `e=u-56d` has norm four.  The norm-four reflection gate gives `u in S2`,
   while `56d in S2`, so `e in S2`, contradicting `min(S2)=12`.

2. At `H=9588+6780sqrt(2)`, norm `q=6`, the Gram type
   `(28734,14400,20340)` has 327 extensions and exactly one nonnegative one:
   `(339,240)`, again of determinant zero, with margin
   `(32-22sqrt(2))/29>0`.  It forces `v=60d`; now `e=u-84d` has norm six
   and lies in `S3`, contradicting `min(S3)=24`.

There are no further nonnegative necessary types through the cutoff.  Thus
the exact ray remains strict against every possible wall through height
25000, so this bounded chamber prefix still contains an open hyperideal
patch.  This is not Vinberg termination and does not prove global
nonreflectivity.

A hostile replay importing no project search code independently reconstructs
the coefficient box, all 301 controls, all 5,575 Gram types, all 2,663,705
extensions, the retained-type stream digest, and both short-vector
reductions.  Its frozen verdict is `PASS`.

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h25000_replacement.py --freeze
python3 scripts/independent_hostile_audit_silver_lambda22_h25000_replacement.py --verify
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h25000_replacement.py
```
