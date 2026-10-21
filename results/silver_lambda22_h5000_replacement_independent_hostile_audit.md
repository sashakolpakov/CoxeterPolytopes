# Silver Lambda22 H5000 replacement independent hostile audit

**Verdict: PASS.**  A clean-room exact replay confirms the bounded frontier
for the ray

`z=(1,(10/31)d)`, with `d^2=4`,

in `<-(sqrt(2)-1)> + Lambda22` in hyperbolic dimension 22.

The replay imports no project search or number-field module.  Its standalone
`Q(sqrt(2))` arithmetic reconstructs the reflective controls, binary
coefficient Gram types, local pairing gates, integral positive-semidefinite
extensions by `d`, and ray pairings.  An exact coefficient-box proof gives
`|x_0|<500` and `|x_1|<3`, hence `|m|<=251`, `|n|<=177` for
`x=m+n sqrt(2)`; the replay exhausts the larger symmetric box 252.

The interval semantics are exact: the lower endpoint
`1598+1130sqrt(2)` is excluded and the rational field element `5000` is
included, with both comparisons at the identity embedding.  The replay finds
exactly:

- 69 new controls;
- 635 retained gate-compatible necessary Gram types;
- seven types whose relaxed maximum can pair nonnegatively with `z`.

There is a subtle logical point here.  Ruling out only the maximizing
extension would not by itself rule out a different nonnegative extension of
the same Gram type.  The independent replay therefore enumerates **every**
locally admissible gated integral PSD extension.  Each of the seven types has
exactly one nonnegative extension.  Thus the production maximizers exhaust
the dangerous extensions, rather than merely bounding them.

For four extensions, `u-v-kd` has norm two in `Lambda22`, contradicting its
minimum four.  For the other three, `k` is even, so `u-v-kd` lies in the
dyadic pairing-radical lattice `S2`; it has norm four, contradicting the
independently certified `min(S2)=12`.  The relevant `k` values may be chosen
as `6,8,6,8,6,10,7` in chronological row order (the fifth row also admits
`k=7`).

Continuity is explicit: the replay checks that the preceding H1598 artifact
uses the same scale, certifies the ray hyperideal and strict through H1598,
and has negative exact pairings with all three accepted cutters.  It also
checks the prior independent lattice-minima artifact before using the
ambient and `S2` short-vector contradictions.

Therefore the production H5000 claim is correct: the same hyperideal ray is
strict against every wall in the necessary superset through identity height
5000.  This is still only a bounded-prefix result.  It proves neither global
nonreflectivity nor Vinberg termination.

Reproduce with:

```sh
python3 scripts/independent_hostile_audit_silver_lambda22_h5000_replacement.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h5000_replacement_independent_hostile_audit.py
```
