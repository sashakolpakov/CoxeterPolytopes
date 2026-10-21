# Silver Lambda22: third accepted cutter at H1598

For `<-(sqrt(2)-1)> + Lambda22` in H22, let `d=e1` and
`e=-e3-e17`.  Then `d^2=4`, `e^2=10`, and `d.e=0`.  The root

`r=(62+44sqrt(2),(20d+e)+sqrt(2)(14d+e))`

is primitive, reflection-integral, has norm 2, and occurs at
`H=1598+1130sqrt(2)`.  It pairs positively with the previous `32d/99`
hyperideal ray, so it cuts that ray.

## Chronological acceptance

The exact 4-by-4 scan covers 284 controls and 1,134 gated post-H14
necessary Gram types.  It has 32 nonnegative rows, 34 extensions, and five
positive extensions.  The two H21 positives and the H274 positive are the
previously certified empty S3/S2 rows.  At H1598, the norm-two positive is
the displayed root.  The only other positive, a norm-four row, forces
`u-v=8d`; putting `h=v-20d` gives `h^2=20`, `e.h=14`, and therefore
`(h-e)^2=2`, impossible in Lambda22.  Thus the displayed wall is accepted.

Complete minimum-four extension scans show it is strictly acute to both H14
packets, and direct calculation gives pairings `-34-24sqrt(2)` and
`-37-26sqrt(2)` with the accepted H164 and H738 walls.

## Replacement ray

The ray `z'''=(1,(10/31)d)` has norm `1361/961-sqrt(2)>0` and pairs
strictly negatively with all three accepted cutters.  A complete gated
necessary scan has 1,138 retained types and only ten relaxed nonnegative
rows.  All ten are excluded by exact ambient/S2/S3/S6 minimum arguments,
the complete H14 shell bound, or the q4 reflection gate.  Hence this ray is
strict through H1598.

The bounded prefix therefore still has infinite volume.  This is not a
finite-volume example, Vinberg termination, or a global nonreflectivity
claim.

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h1598_cutter.py
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h1598_cutter.py --full
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h1598_cutter.py
```
