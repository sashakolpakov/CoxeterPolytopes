# Cubic-seven p7 `D24+`: exact H=300 prefix

## Outcome

For

\[
L=\langle-\theta(2-\theta)\rangle\oplus D_{24}^{+},
\]

the glue-complete exact Vinberg prefix through height 300 still has just one
positive-height wall.  Together with the 24 height-zero `D24` simple roots,
it gives 25 walls and 25 extreme rays, classified exactly as

```text
finite / cusp / hyperideal = 2 / 0 / 23.
```

There is no finite-volume closure through H=300.  This is strictly an
exhaustive bounded-prefix result; it is not a proof of nonreflectivity.

## Complete mirror-norm gate

Put `delta=2-theta`, the totally positive generator of the ramified prime of
norm seven, and `a=theta*delta`.  For a primitive vector `r=(x,y)`, set

```text
I(r) = (a*x, (y,D24+)).
```

Because `D24+` is even unimodular, reflection integrality gives
`(q(r)) | 2I(r)`.  Primitivity and self-duality exclude all primes away from
`(2)` and `(delta)`.  The valuation at the squarefree odd prime `(delta)` is
zero or one.  At the inert prime `(2)`, evenness and the same divisibility
give valuation zero or one.  Since the field has class number one and every
totally positive unit is a square, the complete canonical candidate norms
are

```text
1, 2, delta, 2*delta.
```

For `q=delta` or `q=2*delta`, the spatial reflection condition forces
`y in delta*D24+`: in the first case 2 is a unit modulo `delta`, and in the
second case the factor 2 cancels.  The enumerator therefore writes
`y=delta*z` and exhausts `z in D24+` at square target `y^2/delta^2`.  If that
target is not integral, the control is certified empty before spatial
enumeration.  This is an exact reduction, not a heuristic cutoff.

## Exact census and cone

There are 51 compatible control levels through H=300, broken down by norm as
follows:

| norm | controls |
|---|---:|
| `1` | 4 |
| `2` | 10 |
| `delta` | 9 |
| `2*delta` | 28 |

Only one of the 51 levels contains a reflective dominant candidate, and it
is accepted.  It is the norm-one root-coset affine wall

```text
x = 2 + 3theta + theta^2
h = 7 + 22theta + 10theta^2
y = (1 + 2theta + theta^2)(e1+e2).
```

Every other control through H=300 has zero reflective spatial candidates;
in particular, all 20 controls added beyond the prior H=100 census are
empty.  No spinor-coset wall occurs.  Equal-height candidates are processed
simultaneously, with zero positive pairings among survivors.

Exact double description of the resulting affine cone gives:

| walls | rays | finite | cusp | hyperideal |
|---:|---:|---:|---:|---:|
| 25 | 25 | 2 | 0 | 23 |

Every extreme ray has active-functional rank 24.  Adjacency uses exact rank
23 of common active functionals.  The canonical ray-stream SHA256 is
`d6feb8d1ddfac35443db8332bcfc06f55e93046e63905ee989b78944fbb8cc15`.

The last enumerated control below the bound is

```text
h = 89/2 + 125theta + 111/2 theta^2
q = 4 - 2theta
x = 4 + 9theta + 4theta^2
y^2 = 71 + 187theta + 83theta^2.
```

## Reproducibility

- Backend: `src/coxeter_polytopes/search/cubic7_d24_plus_coefficients.py`
- Audit: `scripts/audit_cubic7_p7_d24_plus_h300.py`
- Test: `tests/test_cubic7_d24_plus_coefficients.py`
- Archive: `data/enumerations/cubic7_p7_d24_plus_h300.json`

The archive records all 51 exact controls, every height packet, the accepted
wall, complete norm proof, exact frontier, cone counts, and ray hash.
