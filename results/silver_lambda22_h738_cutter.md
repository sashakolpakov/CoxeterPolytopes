# Silver Lambda22: accepted H738 cutter, but a new ray survives

For `<-(sqrt(2)-1)> + Lambda22` in H22, the primitive reflective root
`r=(42+30sqrt(2),u+sqrt(2)v)`, with `u=14d+g`, `v=10d+g`,
`d=e1`, and `g=(0,-1,0^19,-1)`, has norm 2 and height
`738+522sqrt(2)`. It cuts the prior `13d/40` ray.

The exact 4-by-4 chronological Schur scan has 517 gate-compatible necessary
types and 23 nonnegative rows. Earlier positive rows are all empty: the
H21/H56/H274 rows by frozen minimum arguments, and H478 because its zero
Schur complement forces nonintegral vectors `(47d+4g)/3,(35d+4g)/3`.
The only remaining positive is the displayed root itself, so it is accepted.

The replacement `z''=(1,(32/99)d)` has norm
`13897/9801-sqrt(2)>0` and is strict against every possible wall through
H738. Therefore this bounded prefix still has infinite volume.

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h738_cutter.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h738_cutter.py
```
