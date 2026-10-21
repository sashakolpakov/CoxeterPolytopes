# Cubic-seven (p_7) `D24+`: exact H=1000 prefix

## Outcome

For

\[
K=\mathbb Q(\theta),\qquad \theta=2\cos(2\pi/7),\qquad
L=\langle-\theta(2-\theta)\rangle\oplus D_{24}^{+},
\]

the glue-complete Vinberg search through exact height (H=1000) has 179
walls: 24 height-zero `D24` walls and 155 positive-height walls.  The old
H=300 cone's 23 hyperideal rays are all cut by the new walls.  Nevertheless,
the refined cone has the exact extreme ray

\[
v=\left(\frac{23+5\theta-3\theta^2}{14},1,1,1,0^{21}\right),
\]

whose Lorentz norm is

\[
(v,v)=\frac{131-54\theta-39\theta^2}{28}>0.
\]

Its active-functional rank is 24 and its exact pairing census against all 179
walls is

```text
negative / zero / positive = 154 / 25 / 0.
```

Thus the exact (H\le1000) prefix still does not have finite volume and is
not compact.  This is strictly a bounded-prefix result: it neither proves
global nonreflectivity nor rules out closure by later Vinberg walls.

## Complete arithmetic gate

Put `delta=2-theta`, the totally positive generator of the ramified prime
above seven, and `a=theta*delta`.  The field has class number one, `(2)` is
inert, `(7)=(delta)^3`, and every totally positive unit is a square.  For a
primitive vector `r=(x,y)`, reflection integrality gives

```text
(q(r)) divides 2 I(r),   I(r)=(a*x,(y,D24+)).
```

Self-duality excludes primes away from `(2)` and `(delta)`.  At each supported
prime the valuation is zero or one.  The complete canonical mirror-norm list
is therefore

```text
1, 2, delta, 2*delta.
```

For `q=delta` or `q=2*delta`, reflection integrality is equivalent to
`y in delta*D24+`.  Writing `y=delta*z` reduces the spatial enumeration
exactly to `z in D24+` at square `y^2/delta^2`; a nonintegral quotient is a
proof that the control is empty.

The `D24+` enumerator uses an even integral Gram matrix of determinant one.
It exhausts the `D24` root coset and all seven nonzero residues of
`O_K/2O_K` in the spinor coset, with the common nonzero parity condition in
all 24 coordinates and both dominant last-coordinate signs.  An independent
hostile audit also reran height zero over all four norm classes and found only
the classical 24 `D24` chamber walls.

## Exact H=1000 census

There are 94 compatible controls, all at distinct exact heights:

| mirror norm | controls |
|---|---:|
| `1` | 6 |
| `2` | 20 |
| `delta` | 18 |
| `2*delta` | 50 |

Only three height packets are nonempty:

| height | candidates | accepted | accepted cosets |
|---|---:|---:|---|
| `7+22theta+10theta^2` | 1 | 1 | `D: 1` |
| `78+220theta+98theta^2` | 77 | 77 | `D: 72`, spinor residues `1: 4`, `2: 1` |
| `122+342theta+152theta^2` | 77 | 77 | `D: 72`, spinor residues `1: 4`, `7: 1` |

The two new packets both have norm 2.  Their controls and spatial-square
targets are respectively

```text
x2 = 6+14theta+6theta^2,   y2^2 = 158+440theta+196theta^2,
x3 = 6+18theta+8theta^2,   y3^2 = 246+684theta+304theta^2.
```

Exact acuteness checks give:

- within packet two: all `C(77,2)=2926` pairings are negative, with maximum
  `1-theta-theta^2 < 0`;
- within packet three: all 2926 pairings are negative, with maximum
  `-2-2theta < 0`;
- packet three against packet two: 5927 negative and 2 zero pairings;
- both new packets have no positive pairing against any earlier wall.

Thus simultaneous insertion of every equal-height survivor is exact.

## Hyperideal-ray audit

Every one of the 23 hyperideal rays of the H=300 affine cone has a positive
pairing with packet two, and likewise with packet three, so none survives.
This was checked by the complete exact `23 x 154` pairing matrix.

The replacement ray (v) above is active on facets

```text
0, 1, 3, 4, ..., 23, 62, 142
```

when facets `0..23` are the height-zero walls, facet 24 is the affine wall,
facets `25..101` are packet two, and facets `102..178` are packet three.  Its
blockwise exact pairing censuses are

| wall block | negative | zero | positive |
|---|---:|---:|---:|
| height zero | 1 | 23 | 0 |
| affine wall | 1 | 0 | 0 |
| packet two | 76 | 1 | 0 |
| packet three | 76 | 1 | 0 |

The 25 active functionals have exact rank 24, so (v) spans an extreme ray.
A floating-point LP was used only to propose its active set.  The archived
certificate stores the field vector itself; replay verifies every pairing,
the rank, and the norm exactly in (K).  No full double-description census
of the 179-wall cone is claimed.

## Reproducibility

- Backend: `src/coxeter_polytopes/search/cubic7_d24_plus_coefficients.py`
- Generator: `scripts/audit_cubic7_p7_d24_plus_h1000.py`
- Tests: `tests/test_cubic7_p7_d24_plus_h1000.py`
- Archive: `data/enumerations/cubic7_p7_d24_plus_h1000.json`

The archive records every exact control, all 155 accepted roots, glue-coset
labels, same-height pairing checks, the full old-ray cut audit, the exact
replacement-ray certificate, and the precise H=1000 frontier.

Two independent reviews returned `PASS` with no defect: a clean exact
pairing/witness audit reproduced the complete `23 x 154` old-ray matrix and
the replacement-ray certificate, while a hostile local-arithmetic audit
independently checked the mirror-norm theorem, delta reduction, glue coverage,
94-control census, and all within/cross-packet acuteness counts.
