# Cubic-seven `D24+` small-coefficient screen through H=100

## Outcome

No compact or finite-volume H24 chamber was found.  All three genuinely new
small ideal squareclasses remain open after their exact compared packets.
The formal winner by the requested hyperideal-ray count is

```text
a = theta*(2-theta): 25 walls, F/C/H = 2/0/23.
```

It is a sparse winner, not evidence of approaching closure: its affine wall
is the only accepted positive wall through H=100.  The denser
`a=2*theta*(2-theta)` lane grows from 23 to 109 and then 398 hyperideal rays.

## No duplicate unit screen

The two apparent unit coefficients are explicitly time-coordinate isometric:

```text
(theta^2-2)^2 * theta = theta-1,
Norm(theta^2-2)=1.
```

Thus `<-theta>+D24+` and `<-(theta-1)>+D24+` differ by an integral unit change
of the time coordinate.  The former was not redundantly recomputed.

## Complete primitive mirror norms

Put `delta=2-theta`, the totally positive generator of the ramified prime of
norm seven.  For a primitive vector `r=(x,y)`, let

```text
I(r) = (a*x, (y,D24+)).
```

Because `D24+` is even unimodular, reflection integrality gives
`(q(r)) | 2 I(r)`.  Away from primes dividing `2a`, primitivity and
unimodularity make `I(r)` a unit.  At the odd squarefree prime `(delta)`,
either `y` is primitive and the norm has valuation zero, or `y` is divisible
by `delta`, `x` is a unit, and

```text
q(r) = -a*x^2 mod delta^2,
```

so the valuation is exactly one.  At inert `(2)`, evenness gives `y^2 in
(2)`.  When `a` is odd, the norm has 2-adic valuation zero or one.  When
`v_2(a)=1`, the entire form is even and the same congruence shows that every
reflective norm has valuation exactly one.  Class number one and `U^+=U^2`
remove totally positive unit-square ambiguity.  Therefore the complete
canonical lists are:

| coefficient | ideal support | complete norms |
|---|---|---|
| `2*theta` | `(2)` | `2` |
| `theta*delta` | `(delta)` | `1, 2, delta, 2*delta` |
| `2*theta*delta` | `(2)(delta)` | `2, 2*delta` |

In particular, the older `{1,2}` list for the even `2*theta` host was safe
but overinclusive; norm one cannot occur.

## Exact H=100 census

The spatial search uses the genuine orthonormal model

```text
D24+ = D24 union ((1/2,...,1/2)+D24)
```

and exhausts the root coset plus all seven nonzero residues in `O_K/2O_K`
for the spinor coset.  Dominance is
`y1>=...>=y23>=abs(y24)`.  No spinor-coset wall is accepted in any seed.
Equal-height survivors are inserted simultaneously and are pairwise acute.

| coefficient | controls | nonempty packets | accepted walls | accepted spinor |
|---|---:|---:|---:|---:|
| `2*theta` | 5 | 4 | 10 | 0 |
| `theta*delta` | 31 | 1 | 1 | 0 |
| `2*theta*delta` | 8 | 3 | 8 | 0 |

The first packet of every seed is one affine highest-root wall.  Its exact
cone has 25 rays and `F/C/H=2/0/23`.

| coefficient | packet | exact height | candidates / accepted | walls | rays | F/C/H |
|---|---:|---|---:|---:|---:|---:|
| `2*theta` | 1 | `2theta+theta^2` | 1 / 1 | 25 | 25 | 2/0/23 |
| `2*theta` | 2 | `3+9theta+4theta^2` | 3 / 2 | 27 | 113 | 4/0/109 |
| `theta*delta` | 1 | `7+22theta+10theta^2` | 1 / 1 | 25 | 25 | 2/0/23 |
| `2*theta*delta` | 1 | `2+7theta+3theta^2` | 1 / 1 | 25 | 25 | 2/0/23 |
| `2*theta*delta` | 2 | `7+22theta+10theta^2` | 2 / 2 | 27 | 113 | 4/0/109 |
| `2*theta*delta` | 3 | `12+34theta+15theta^2` | 5 / 5 | 32 | 404 | 6/0/398 |

For `2*theta`, the complete H=100 wall census also has later packets of
2, 5, and 0 accepted walls, but cone geometry was intentionally stopped
after packet two because this coefficient was already the poor control.

Every reported cone ray has an active-functional rank of 24.  Double
description uses the exact rank-23 common-active adjacency test; no numerical
geometry or full-dimensional floating-point hull is used.

## Reproducibility and scope

- Backend: `src/coxeter_polytopes/search/cubic7_d24_plus_coefficients.py`
- Glue DP: `src/coxeter_polytopes/search/cubic7_dn.py`
- Audit: `scripts/audit_cubic7_d24_plus_coefficient_screen.py`
- Tests: `tests/test_cubic7_d24_plus_coefficients.py`
- Archive: `data/enumerations/cubic7_d24_plus_coefficient_screen_h100.json`

The archive stores every exact control and accepted root, packet DD profiles,
causal counts, and canonical ray-stream hashes.  Its conclusion is only that
these exhaustive H=100 prefixes have not closed.  It is not a proof of global
nonreflectivity.
