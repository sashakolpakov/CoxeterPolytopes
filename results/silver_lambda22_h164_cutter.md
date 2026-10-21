# Silver Lambda22: first accepted cutter at H164, but no closure

For the arithmetic host

`<-(sqrt(2)-1)> + Lambda22`

in hyperbolic dimension 22, the first explicit root found that cuts the
frozen ray `z=(1,d/3)` occurs at

`H=164+116sqrt(2)`.

Take `d=e1` and

`e=(-6,31,21,15,10,20,14,7,4,2,2,1,0^10)`.

Then `d^2=4`, `e^2=12`, `d.e=-2`, and

`u=8d-e`, `v=6d-e`, `r=(28+20sqrt(2),u+sqrt(2)v)`

has `r^2=4`.  It is primitive and reflection-integral.  Its pairing with
the old ray is `(-2+2sqrt(2))/3>0`, so it is a genuine cutter.

## Chronological acceptance

An exact necessary-superset scan uses every reflective control through H164,
every binary coefficient Gram type surviving the local gates and minima,
and every positive-semidefinite 4-by-4 extension by the fixed pair `(d,e)`.
There are nine pre-H164 nonnegative rows, all tangent, and no positive row.
At H164 there are two q=2 tangent rows and one positive q=4 extension.  The
positive extension has zero Schur complement and is exactly the displayed
root.  Hence no earlier wall rejects it.

## Replacement ray

The cutter still does not close the chamber.  The exact ray

`z'=(1,(13/40)d)`

has norm `569/400-sqrt(2)>0`.  A complete gated necessary-superset scan of
91 controls and 112 gate-compatible binary types through H164 leaves only four relaxed
nonnegative rows; all four are already excluded by the frozen H14/H150
shell and minimum arguments.  The explicit cutter pairs with `z'` by
`-19/20+9sqrt(2)/20<0`.

Thus the exact prefix through H164 still has infinite volume.  This is a
bounded-prefix result, not a global nonreflectivity theorem.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h164_cutter.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h164_cutter.py
```
