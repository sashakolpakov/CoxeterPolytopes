# Silver prime-seven `D24+`: exact H=100 prefix

## Outcome

For

\[
L=\langle-(2\sqrt2-1)\rangle\oplus D_{24}^{+},
\]

the glue-complete Vinberg search through algebraic height 100 does **not**
produce a finite-volume chamber.  The full 55-wall prefix has an explicit
rank-24 hyperideal extreme ray.  Thus this coefficient is not a new compact
or finite-volume example at the searched bound.

## Complete primitive mirror norms

Write

```text
pi = sqrt(2),  epsilon = 1+sqrt(2),  delta = 3+sqrt(2),
a = 2sqrt(2)-1,  epsilon*a = delta,  Norm(delta)=7.
```

For a primitive vector `r=(x,y)`, reflection integrality gives

```text
(q(r)) | 2 I(r),   I(r)=(a*x,(y,D24+)).
```

Even unimodularity excludes support away from `(pi)` and `(delta)`.  At the
odd prime `(delta)`, a divisible norm forces `y in delta*D24+`; primitivity
then makes `x` a unit and

```text
q(r) = -a*x^2 mod delta^2,
```

so the valuation is exactly one.  At `(pi)`, evenness and primitivity leave
valuation zero or two, the latter being the ideal `(2)`.  Class number one
and the fact that every totally positive unit of `Z[sqrt(2)]` is a square
give the complete canonical list

```text
1, 2, 3+sqrt(2), 6+2sqrt(2).
```

For either prime-seven norm, the implementation first divides the spatial
target by `delta^2` and exhausts `y/delta` in `D24+`.  This is both the local
divisibility proof and a substantial exact-enumeration reduction.

## Glue and primitivity

The spatial search uses the genuine orthonormal model

```text
D24+ = D24 union ((1/2,...,1/2)+D24)
```

over `Z[sqrt(2)]`.  It exhausts the D24 root coset and all three nonzero
residues in `O_K/2O_K` for the spinor coset, with dominance
`y1>=...>=y23>=abs(y24)`.  The integral glue-basis Gram matrix is even,
unimodular, and has determinant one.

Primitivity is tested in that glue basis.  At each nonempty height there is
one raw norm-two vector divisible by `sqrt(2)`.  Its quotient is exactly the
norm-one mirror at the same height, so it is removed rather than counted as
a second root.

## Exact H=100 census

There are 25 exact control levels, grouped into 17 distinct heights.  Only
two packets are nonempty:

| height | raw reflective | imprimitive | primitive candidates | accepted | accepted spinor |
|---|---:|---:|---:|---:|---:|
| `5+4sqrt(2)` | 4 | 1 | 3 | 3 | 0 |
| `31+22sqrt(2)` | 30 | 1 | 29 | 28 | 5 |

Every control for norms `3+sqrt(2)` and `6+2sqrt(2)` is empty after the
required prime-seven divisibility and complete glue enumeration.  Every
height after the second packet through H=100 is also empty.

The first packet consists of

```text
q=1, x=1+sqrt(2), y=(1+sqrt(2),1+sqrt(2),0^22),
q=2, x=2+sqrt(2), y=((1+sqrt(2))^4,0^20),
q=2, x=2+sqrt(2), y=(2+2sqrt(2),0^23).
```

All three are mutually acute and survive simultaneously.  Exact double
description, with adjacency certified by rank 23 of common active
functionals, gives

| walls | rays | finite | ideal | hyperideal |
|---:|---:|---:|---:|---:|
| 27 | 113 | 4 | 0 | 109 |

The frozen projective ray-stream hash is

```text
b5066aaa308ccfc7a48fb8c5897f4fe4613e4683364c5d0fd2def8d81158d24d
```

## Full-prefix hyperideal certificate

Expanding all 28 walls of the second packet is unnecessary for the
finite-volume decision.  Two first-packet hyperideal extreme rays satisfy
every later half-space.  The archived certificate chooses

```text
z = ((20+12sqrt(2))/7;
     2+2sqrt(2), 2+2sqrt(2), 2, 2, 0^20),

z^2 = (48-16sqrt(2))/7 > 0.
```

Its pairings with all 28 later accepted walls are nonpositive; five are
zero.  Adding those zero facets leaves its active-functional rank equal to
24, so `z` remains an extreme hyperideal ray of the full 55-wall H=100
cone.  This proves that the full enumerated prefix has infinite volume.

## Reproducibility and scope

- Backend: `src/coxeter_polytopes/search/silver_p7_d24_plus.py`
- Audit: `scripts/audit_silver_p7_d24_plus_h100.py`
- Tests: `tests/test_silver_p7_d24_plus.py`
- Archive: `data/enumerations/silver_p7_d24_plus_h100.json`

The archive contains every control level, glue count, imprimitive vector,
accepted wall, all 113 first-packet rays, and the 28 exact pairings in the
full-prefix witness.  This is an exhaustive H=100 obstruction, not a global
nonreflectivity theorem.
