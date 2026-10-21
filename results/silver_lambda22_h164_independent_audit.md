# Independent silver Lambda22 H164 wall audit

The explicit vector

`r=(28+20sqrt(2), u+sqrt(2)v)`,

with `d=e1`, the frozen catalogue vector `e`, `u=8d-e`, and `v=6d-e`,
is a primitive reflective root of norm 4 at height `164+116sqrt(2)`.
Its pairing with the old hyperideal ray `z=(1,d/3)` is
`-2/3+2sqrt(2)/3>0`, so it genuinely cuts that ray.

## Chronological acceptance

The independent replay enumerates 91 exact controls and
108 gated necessary binary Gram
types after the complete H14 packet.  Exact 4-by-4 Schur-complement enumeration
finds no strictly positive pairing before H164.  There are
9 possible earlier tangencies, all harmless.
At H164 there are two possible q=2 tangencies and one positive q=4 extension.
The latter has Gram `(300,180,232)`, cross data `(34,26,-28,-24)`, and zero
Schur residual.  Hence it forces `u=8d-e,v=6d-e`: it is the root paired with
itself, with pairing 4.  Therefore the displayed root is accepted.

Against the two complete H14 shells, the exact pairing maxima are respectively
`-12-10sqrt(2)` and `-2-4sqrt(2)`, both strict.

## Replacement hyperideal ray

The ray `z'=(1,13d/40)` has norm `569/400-sqrt(2)>0`.  After all binary
S2/S3/S6 minimum and residue gates, 62
necessary types remain through H164.  Its largest possible pairing on that
superset is `(5-6sqrt(2))/40<0`.  In particular, its pairing with the new
H164 wall is `-19/20+9sqrt(2)/20<0`.  Thus the enlarged bounded prefix still
has infinite volume; this is not a global nonreflectivity result.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/independent_audit_silver_lambda22_h164.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h164_independent_audit.py
```

Source: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/LAMBDA22.html
