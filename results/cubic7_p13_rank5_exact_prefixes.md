# Exact rank-five conductor-seven searches at the primes over 13

## Outcome

The exact weighted-diagonal Vinberg backend now covers the three forms

\[
\operatorname{diag}(-\theta,1,1,1,t),\qquad
t\in\{-\theta+2\theta^2, 1-\theta+\theta^2, 1+\theta^2\},
\]

where \(\theta=2\cos(2\pi/7)\).  These are totally positive generators of
the three primes over 13.  All three searches are exact through height 300;
the third is exact through height 1000.  None closes, and none has stabilized,
so this lane has produced no new finite-volume polytope yet.

The literal height-20--300 artifact is
[`data/enumerations/k7_nonunimodular_rank5/p13_adaptive_h20_h300.json`](../data/enumerations/k7_nonunimodular_rank5/p13_adaptive_h20_h300.json).
The height-500/1000 summary retained after the large in-memory pools were
released is
[`data/enumerations/k7_nonunimodular_rank5/p13_r8_h500_h1000.json`](../data/enumerations/k7_nonunimodular_rank5/p13_r8_h500_h1000.json).

## Complete mirror-norm alphabet

For a primitive vector \(r\), put

\[
I(r)=(r_0,r_1,r_2,r_3,t r_4).
\]

Reflection integrality is equivalent to \((r,r)\mid2I(r)\).  Because
\((t)\) is a prime of norm 13, \(I(r)\) is either \(O_K\) or \((t)\).  In
the latter case the first four coordinates are divisible by \((t)\), while
\(r_4\) is a local unit, and

\[
(r,r)\equiv t r_4^2\pmod {t^2}.
\]

Thus the \((t)\)-valuation of the norm is exactly one.  The rational prime 2
is inert and occurs to exponent at most one.  Consequently every primitive
mirror, modulo unit scaling, has one of the four norms

\[
\boxed{1,\ 2,\ t,\ 2t}.
\]

Narrow class number one and squarehood of every totally positive unit make
these element representatives complete.  The implementation also checks
reflection integrality coordinatewise for every candidate.

Each height-zero system was independently exhausted: exactly 20 roots and
type \(B_3\times A_1\).  Roots at equal exact height are tested against only
lower-height walls, verified pairwise nonpositive, and accepted
simultaneously.

## Prefix census

Entries in the last column are `(finite, cusp, hyperideal)` exact cone rays.

| prime generator | height | candidates | accepted + height-zero | rays |
|---|---:|---:|---:|---:|
| \(-\theta+2\theta^2\) | 20 | 12 | 3 + 4 | (7, 0, 5) |
|  | 50 | 23 | 4 + 4 | (7, 0, 9) |
|  | 100 | 48 | 9 + 4 | (11, 0, 21) |
|  | 200 | 106 | 13 + 4 | (12, 0, 33) |
|  | 300 | 171 | 20 + 4 | (13, 0, 54) |
| \(1-\theta+\theta^2\) | 20 | 10 | 3 + 4 | (7, 0, 5) |
|  | 50 | 21 | 4 + 4 | (7, 0, 9) |
|  | 100 | 46 | 10 + 4 | (8, 0, 29) |
|  | 200 | 97 | 18 + 4 | (9, 0, 55) |
|  | 300 | 161 | 21 + 4 | (9, 0, 72) |
| \(1+\theta^2\) | 20 | 11 | 3 + 4 | (7, 0, 5) |
|  | 50 | 23 | 6 + 4 | (8, 0, 14) |
|  | 100 | 44 | 8 + 4 | (8, 0, 21) |
|  | 200 | 101 | 16 + 4 | (11, 0, 40) |
|  | 300 | 169 | 18 + 4 | (13, 0, 43) |
|  | 500 | 316 | 27 + 4 | (17, 0, 70) |
|  | 1000 | 729 | 52 + 4 | (28, 0, 158) |

At height 1000 the accepted norm histogram is

\[
47\times 2,qquad 3\times1,qquad 2\times(1+\theta^2),
\]

with no accepted \(2(1+\theta^2)\)-wall.  The newest accepted height is

\[
\frac{289}{2}+\frac{811}{2}\theta+\frac{361}{2}\theta^2.
\]

This continued wall growth and the 158 exact hyperideal rays prove only that
the height-1000 prefix is not finite-volume.  They do not prove eventual
nonreflectivity.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/search_cubic7_nonunimodular_seeds.py \
  --heights 20,50,100,200,300 --seed p13 \
  > data/enumerations/k7_nonunimodular_rank5/p13_adaptive_h20_h300.json

PYTHONPATH=src pytest -q tests/test_cubic7_nonunimodular_seed_search.py
```

The focused test result is `7 passed`.

## Cheapest clean mixed-support continuation

The next seed whose complete norm proof follows immediately is an odd
squarefree product \(t=\delta p_{13}\), with norm 91.  For each of its three
conjugate choices the eight canonical norms are

\[
1,2,\delta,p_{13},\delta p_{13},2\delta,2p_{13},2\delta p_{13}.
\]

The coefficient of norm 56 supported on \((2)(\delta)\) is smaller, but the
factor 2 in the reflection divisibility condition interacts with the
2-adic valuation of the diagonal coefficient.  Its complete norm alphabet
needs a separate local audit and is not covered verbatim by the odd
squarefree-prime proof.

