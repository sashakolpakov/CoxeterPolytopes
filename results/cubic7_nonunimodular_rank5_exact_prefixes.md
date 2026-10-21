# Exact rank-five searches over the conductor-seven cubic field

## Outcome

Two first non-unimodular arithmetic lattices over

\[
K=\mathbf Q(\theta),\qquad \theta=2\cos(2\pi/7),\qquad
\theta^3+\theta^2-2\theta-1=0
\]

have now been searched exactly through Vinberg height 300:

\[
L_2=\operatorname{diag}(-\theta,1,1,1,2),\qquad
L_7=\operatorname{diag}(-\theta,1,1,1,\delta),\quad
\delta=2-\theta.
\]

There is **no new finite-volume chamber yet**.  Neither prefix satisfies
Vinberg stopping and neither wall sequence has stabilized.  This is an exact
bounded result, not a claim that either lattice is nonreflective.

The frozen record is
[`data/enumerations/k7_nonunimodular_rank5/adaptive_h20_h300.json`](../data/enumerations/k7_nonunimodular_rank5/adaptive_h20_h300.json).
Its file SHA-256 is

```text
c24f8ca90a0705c4923b5e3dfae56460e8eeb90f7d3375ef65ea01ebfead318d
```

and its canonical payload SHA-256 is

```text
7ee9fd3c7b26b942582fe79db5708f69faca87e2dd98d15759978c29b7f40612
```

## Complete primitive mirror norms

For a primitive vector (r=(r_0,\ldots,r_4)), let

\[
I(r)=(r_0,r_1,r_2,r_3,t r_4),
\]

where (t=2) or (delta).  Its reflection preserves the standard
\(O_K\)-lattice precisely when

\[
(r,r)\mid 2I(r).
\]

The prime 2 is inert in (K), and ((\delta)) is the norm-seven ramified
prime.  Away from the support of (t), primitivity makes (I(r)) the unit
ideal.  At a prime occurring once in (t), either (I(r)=O_K), or the first
four coordinates are divisible by that prime and (r_4) is a local unit. In
the latter case the diagonal norm equation modulo the square of the prime
shows that its valuation in ((r,r)) is exactly one.  Thus the possible norm
ideals are

\[
L_2:\quad (1),(2),
\qquad
L_7:\quad (1),(2),(\delta),(2\delta).
\]

The narrow class number of (K) is one and every totally positive unit is a
square.  Therefore every primitive mirror has, up to scaling its normal by a
unit, exactly one of the element norms

\[
L_2:\quad 1,2,
\qquad
L_7:\quad 1,2,\delta,2\delta.
\]

These lists are complete; they are not guessed small-norm cutoffs.  Every
enumerated vector is additionally checked coordinate by coordinate for
reflection integrality.

## Exact enumeration and same-height audit

The three real embeddings give simultaneous positive-definite coordinate
bounds at the two conjugate places and a prescribed height bound at the
Lorentzian place.  A certified trace-dual coefficient box enumerates all
algebraic integers inside those embedding boxes.  All sign decisions,
weighted sum-of-squares equations, reflection tests, and wall pairings are
over (K), without floating-point acceptance tests.

Height zero was independently exhausted before choosing a sector.  Both
seeds have exactly 20 normalized height-zero roots and finite system

\[
B_3\times A_1.
\]

Positive-height roots are grouped by their **exact** height.  Every member of
a height packet is first tested only against lower-height walls; all survivors
are then checked pairwise nonpositive and accepted simultaneously.  Through
height 300 there are no same-height conflicts.  In particular, the result is
independent of the secondary tuple order used for serialization.  The largest
accepted packets have sizes four for (L_2) and three for (L_7).

## Adaptive checkpoints

The ray triples below are `(finite, cusp, hyperideal)` for the exact current
wall cone.  A hyperideal ray proves only that the displayed prefix has not
closed; a later mirror can remove it.

| seed | height | candidates | walls | cone rays |
|---|---:|---:|---:|---:|
| (L_2) | 20 | 11 | 8 | (9, 0, 6) |
| (L_2) | 50 | 26 | 12 | (13, 0, 12) |
| (L_2) | 100 | 58 | 17 | (19, 0, 24) |
| (L_2) | 150 | 92 | 22 | (30, 0, 32) |
| (L_2) | 200 | 122 | 26 | (38, 0, 42) |
| (L_2) | 300 | 203 | 32 | (49, 0, 47) |
| (L_7) | 20 | 15 | 8 | (9, 0, 7) |
| (L_7) | 50 | 33 | 10 | (10, 0, 14) |
| (L_7) | 100 | 72 | 13 | (18, 0, 19) |
| (L_7) | 150 | 114 | 15 | (21, 0, 22) |
| (L_7) | 200 | 156 | 17 | (21, 0, 27) |
| (L_7) | 300 | 250 | 21 | (23, 0, 48) |

The exact ray census uses incremental double description.  It was checked
facet-for-facet and ray-kind-for-ray-kind against the independent brute
four-facet cofactor enumerator for both seeds at heights 3, 7, and 20.

## Negative stabilizer probe so far

The natural Pell isometries on the ((r_0,r_4))-plane do not yet provide an
infinite chamber stabilizer.  For (L_7),

\[
u=v=2\theta+\theta^2,qquad u^2-\theta\delta v^2=1,
\]

folds back to the identity through the prefix reflection group.  For (L_2),

\[
u=11+28\theta+12\theta^2,qquad
v=6+18\theta+8\theta^2,qquad u^2-2\theta v^2=1,
\]

does likewise, by the exact wall-reflection word

```text
[7, 0, 1, 2, 1, 0, 7]
```

in the height-300 wall ordering.  Hence these smallest obvious Pell elements
lie in the reflection group; they do not prove nonreflectivity.

## Replay

```bash
PYTHONPATH=src python3 scripts/search_cubic7_nonunimodular_seeds.py \
  --heights 20,50,100,150,200,300 --seed both \
  > data/enumerations/k7_nonunimodular_rank5/adaptive_h20_h300.json

PYTHONPATH=src pytest -q \
  tests/test_cubic7_exact_and_discovery_slices.py \
  tests/test_cubic7_nonunimodular_seed_search.py
```

The focused test result at freeze time was `11 passed`.

