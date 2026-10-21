# Silver Lambda22: fourth accepted cutter at H2786

For `<-(sqrt(2)-1)> + Lambda22` in H22, put `d=e1` and
`e=e2+e22` in the catalogue basis.  Their Gram matrix is `[[4,2],[2,10]]`.
The primitive reflection-integral norm-two root

`r=(82+58sqrt(2),(26d+e)+sqrt(2)(18d+e))`

has height `2786+1970sqrt(2)`.  It pairs positively with the former
`(1,10d/31)` hyperideal ray, so it genuinely cuts that ray.

The complete chronological necessary scan contains 374 controls, 1,997
post-H14 gate-compatible Gram types, 23 nonnegative rows, 26 extensions, and
eight positive extensions.  Six are the previously certified S2/S3/S6
minimum contradictions.  The only new earlier positive is a q4 row at
`H=2404+1700sqrt(2)`.  Its Schur equations force `u-v=10d`; then
`w=9u-304d-13e` would have norm 18, be orthogonal to `d,e`, and satisfy
`w=2d-4e mod 9`.  An exact Smith/LLL reduction leaves 13,824 vectors inside
the dual coordinate bounds, and none has norm 18.  Thus that row is empty.
The last positive is exactly the displayed root, which is therefore accepted.

Complete H14 shell maximization and direct calculation show strict acuteness
to both H14 packets and the accepted H164, H738, and H1598 cutters.

The old `10d/31` ray is replaced by `z=(1,(28/87)d)`.  Its exact norm is
`10705/7569-sqrt(2)>0`, and it lies strictly behind all four accepted later
walls.  The complete 374-control necessary-wall replay retains 2,001 gated
types and only five relaxed threats.  The H14 threat is removed by the exact
S2 shell bound; three are the frozen S3/S6/S2 minimum contradictions; and
the H1202 q4 threat violates its reflection gate.  Thus `z` is strict against
every possible wall through H2786, so this complete bounded prefix still has
an open hyperideal patch.

This is not a finite-volume example, a Vinberg-termination claim, or a global
nonreflectivity theorem.

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h2786_cutter.py --freeze
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h2786_cutter.py --start 340 --stop 350
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h2786_cutter.py
```
