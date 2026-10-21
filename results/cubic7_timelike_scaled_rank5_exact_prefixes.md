# Exact timelike-scaled K7 rank-five prefixes

## Outcome

Neither timelike-scaled seed closes to a finite-volume Vinberg chamber in the
certified prefixes.  The frozen height-300 checkpoints are:

| lattice | candidates | accepted positive-height walls | total walls | finite rays | cusps | hyperideal rays |
|---|---:|---:|---:|---:|---:|---:|
| `T2 = diag(-2 theta,1,1,1,1)` | 45 | 22 | 26 | 36 | 0 | 33 |
| `T7 = diag(-theta(2-theta),1,1,1,1)` | 44 | 9 | 13 | 7 | 0 | 28 |

Thus these two determinant directions have not produced a new finite-volume
example.  This is a bounded exhaustive result, not a proof of
nonreflectivity: a later wall can cut off a current hyperideal cone ray.

## Why the mirror-norm lists are complete

Write `delta=2-theta` and let `s` be `2` or `delta`.  For a primitive vector
`r=(x,y1,...,y4)`, reflection integrality is equivalent to

`q(r) | 2 I(r)`, with `I(r)=(-s theta x,y1,...,y4)`.

The field has discriminant 49.  Its Minkowski bound is
`(3!/3^3)sqrt(49)=14/9<2`, so its class number is one.  The exact sign vectors
of the units `-1`, `theta`, and `theta^2-2` have rank three over `F_2`:

`(---), (+--), (--+)`.

Consequently the unit-signature map is onto, the narrow class number is one,
and every totally positive unit is a square.

The polynomial `x^3+x^2+1` has no root over `F_2`, so `(2)` is inert, while
`N(delta)=7`.  Away from primes dividing `2s`, primitivity makes `I(r)` the
local unit ideal.  At the squarefree prime dividing `s`, either a spatial
coordinate is a unit, or `x` is a unit and every spatial coordinate is in
that prime.  In the latter case the diagonal congruence gives valuation one
for `q(r)`.  This leaves exactly the following possible ideals modulo totally
positive unit squares:

- T2: `(1)`, `(2)`;
- T7: `(1)`, `(2)`, `(delta)`, `(2 delta)`.

The code still checks coordinatewise reflection integrality for every vector.
Because these ideals are squarefree, a vector of a canonical norm cannot have
a nonunit common divisor.

## Exact search gates

The height-zero search independently exhausts every normalized mirror and
finds precisely the 32 roots of `B4`: 8 norm-one coordinate roots and 24
norm-two signed-pair roots.  The positive-height enumerator uses certified
trace-dual coefficient boxes and exact signs at all three real embeddings.
It solves four squares only in the dominant chamber
`y1 >= y2 >= y3 >= y4 >= 0`.

All candidates at one exact height are first filtered against lower-height
walls; the survivors are proved mutually acute and then adjoined
simultaneously.  Extreme rays are computed by exact incremental double
description.  At height 50, an independent enumeration of every four-facet
support gives exactly the same projective rays for both seeds.

## Adaptive search record

The columns are candidate mirrors, accepted positive-height walls, total
walls, and `(finite,cusp,hyperideal)` cone rays.

| height | T2 | T7 |
|---:|---|---|
| 20 | `1 / 1 / 5 / (2,0,3)` | `3 / 1 / 5 / (2,0,3)` |
| 50 | `7 / 5 / 9 / (9,0,10)` | `11 / 5 / 9 / (5,0,13)` |
| 100 | `13 / 7 / 11 / (11,0,15)` | `16 / 6 / 10 / (6,0,15)` |
| 200 | `38 / 19 / 23 / (27,0,27)` | `39 / 8 / 12 / (6,0,24)` |
| 300 | `45 / 22 / 26 / (36,0,33)` | `44 / 9 / 13 / (7,0,28)` |
| 500 | `91 / 33 / 37 / (55,0,59)` | `92 / 15 / 19 / (11,0,44)` |
| 1000 | `223 / 56 / 60 / (97,0,101)` | `260 / 29 / 33 / (25,0,82)` |

The height-1000 passes were feasible and completed, but the compact frozen
artifact below deliberately contains the requested H100/H300 checkpoints.
Its per-run hashes cover every accepted root and a hyperideal witness.

## Independent live-AlVin comparison

After exact conversion from AlVin's cyclic RC7 basis to the power basis, the
first five positive-height walls `e5,...,e9` in each interrupted live job are
exactly the first five walls from this backend (unordered within equal-height
packets).  The frozen artifact records the source transcript hashes:

- T2 stdout: `1214da006c59186c1d7777bc6d34beb6dff1a6b7e043e828f18750c7e8d56540`;
- T7 stdout: `17ffce322e2495b5b88d755e0771f417dc5589b0b04274a43223329bc7d5419d`.

Those AlVin processes were externally interrupted, so this comparison is an
early-enumeration regression only.  It is not evidence of AlVin termination.

## Replay and frozen hashes

Run:

```text
PYTHONPATH=src python3 scripts/search_cubic7_timelike_scaled.py \
  --heights 100,300 \
  --output data/enumerations/k7_timelike_scaled_rank5_h100_h300.json
PYTHONPATH=src pytest -q tests/test_cubic7_timelike_search.py
```

The dedicated suite passes: `5 passed`.

- backend: `d5e291e3237b9be2e1030c710cad97eed911dd992eb1e148d5b0a62c2b5f53b7`;
- generator: `399d44fc9feecd0771140a28c4ca1c9de2b8aae161f73cc7e0bbc43b0d1c9bc9`;
- tests: `e77b0174bcb26ca6e5c779b3626c82a925387b9b792bd6a8e2b8e5ea067a7a87`;
- artifact file: `924949d1cab78c7816e1205662d7502e85c861fb353810b11838cd4fd63a9d34`;
- canonical artifact payload: `94f30e57c1edcd80de327a1be759b805c4942859c2bd28b5979e44bd5bba9fa5`.
