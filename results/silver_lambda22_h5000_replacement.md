# Silver Lambda22 replacement survives through real height 5000

For `<-(sqrt(2)-1)> + Lambda22` in `H22`, the exact ray

`z=(1,(10/31)d)`,  `d^2=4`,

has norm `1361/961-sqrt(2)>0`.  Starting immediately after the accepted
wall at `1598+1130sqrt(2)`, the complete exact control scan through the
rational field-element bound `H=5000` contains 69 new controls and 635
gate-compatible necessary Gram types.

Only seven relaxed types can pair nonnegatively with `z`.  Four force a
nonzero Lambda22 vector of norm two.  The other three force a vector in the
dyadic radical sublattice `S2` of norm four, below `min(S2)=12`.  Therefore
all seven are empty and `z` is strict against every possible wall through
the cutoff.

The bounded prefix still has infinite volume.  This is not Vinberg
termination and does not prove global nonreflectivity.

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h5000_replacement.py
```
