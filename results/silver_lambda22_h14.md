# Silver Lambda22 exact H14 prefix

The exact host

`<-(sqrt(2)-1)> + Lambda22`

in hyperbolic dimension 22 has a complete Vinberg prefix through
`H=14+10sqrt(2)`.  It is **not** finite-volume at this cutoff: for any
minimal `d in Lambda22`, the exact ray `z=(1,d/3)` has

`z^2 = 13/9-sqrt(2) > 0`

and satisfies every frozen wall strictly.  This is a prefix result, not a
global nonreflectivity theorem.

## Exact lattice and local gates

The vendored Nebe--Sloane Gram matrix has determinant 12, minimum 4,
kissing number 49,896, and Smith form `1^20,2,6`.  The complete primitive
mirror-norm alphabet is

`1, 2, 4+2sqrt(2), 4, 3, 6, 12+6sqrt(2), 12`.

The dyadic valuations are `0,2,3,4`; exponent 1 is excluded by evenness.
For `y=u+sqrt(2)v`, valuation 3 requires `G*u` even, valuation 4 requires
both `G*u,G*v` even, and the 3-primary classes require both divisible by 3.
Independent HNF+PARI checks give minima `min(S2)=12`, `min(S3)=24`, and
`min(S6)=96`.

## First simultaneous packet

At `H=14+10sqrt(2)` there are exactly two shell descriptions:

1. `q=2`, `x=6+4sqrt(2)`, and `y=(1+sqrt(2))w` for the complete ordinary
   `Lambda22` norm-10 shell.  Its exact cardinality is 467,596,800.
2. `q=4`, `x=8+6sqrt(2)`, and `y=(1+sqrt(2))u` for the complete norm-20
   shell of `S2={u:G*u even}`.  Its exact signed cardinality is 1,368,576.

The simultaneous packet therefore contains 468,965,376 roots.

All other necessary Gram types through the same height are killed by the
stated gate minima or mod-8 radical residues.  The two full components are
mutually Vinberg-acute.  In particular, cross dot 14 would give a forbidden
norm-2 vector, so the cross dot is at most 13.

For the witness, ordinary-shell pairings use `|d.w|<=5`; dyadic-shell
pairings use `|d.u|<=8`.  The respective worst exact pairings are
`-(1+sqrt(2))/3` and `-4/3+2sqrt(2)/3`, both negative.

## Frontier

The next sharp target is the first unresolved possible tangency at
`H=36+51sqrt(2)/2`, in the mixed norm `12+6sqrt(2)` class.  That continuation
is not part of this frozen H14 certificate.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h14.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22.py
```

Source: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/LAMBDA22.html
