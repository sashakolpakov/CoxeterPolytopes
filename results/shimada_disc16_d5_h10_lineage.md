# Lineage of the Shimada D5 involution chamber in H10

## Result

The noncentral D5 involution produces a genuine finite-volume, noncompact
hyperbolic Coxeter polytope in `H10` with 30 facets.  Its lineage is now
exact:

- the saturated fixed lattice is `2*I_(1,10)`;
- all 30 facet mirrors are primitive integral reflecting roots;
- those 30 roots generate the entire fixed lattice;
- the ambient maximal root-reflection group is Vinberg's classical `P10`,
  whose chamber has 12 facets;
- the 30-wall group has exact index

  ```text
  [W(P10) : W30] = 12,142,080
                    = 2^9 * 3^2 * 5 * 17 * 31;
  ```

- the literal 30-wall weighted diagram is absent from the pinned 1,214-file
  ordinary CoxIter corpus and the separate 2,959-file reflection-doubling
  archive.

This is a positive algebraic construction and a possible priority candidate
for the *literal subgroup chamber*.  It is **not** a new lattice, maximal
reflection group, arithmetic commensurability class, or first example in
dimension 10.  The bounded corpus exclusion is **not** a worldwide priority
claim.

The replayable certificate is
[`shimada_disc16_d5_h10_30_lineage.json`](../data/enumerations/shimada_disc16_d5_h10_30_lineage.json).

## Hyperbolic-polytope gate

Two exact finite-volume tests are available upstream.

1. The GMP rational double-description has 1,540 extreme rays:
   1,333 finite vertices, 207 ideal vertices, and zero hyperideal rays.
   The cone is pointed and full-dimensional.
2. The literal Vinberg audit enumerates all 10,030 elliptic rank-9 edge
   diagrams.  Every edge has exactly two finite/ideal endpoints, with

   ```text
   FF = 4,423,  FI = 4,484,  II = 1,123.
   ```

   It also enumerates all 436 connected affine components; every full
   zero-neighborhood is a valid rank-9 cusp neighborhood.  There are no bad
   edges and no bad parabolic neighborhoods.

Thus the chamber is finite-volume and noncompact, not merely a locally valid
Coxeter diagram.

The complete face vector, vertices first, is

```text
(1540, 10030, 25206, 34857, 30452,
 17784,  7072,  1890,   321,    30, 1).
```

## Saturated divisor-fixed lattice

Let `M` be the archived noncentral involution acting on Shimada's integral
rank-20 Neron--Severi lattice.  Row vectors are fixed by `M`, so the fixed
lattice is

```text
Fix(M) = ker_Z(M^T - I).
```

The Smith kernel is primitive and has rank 11.  In the reduced saturated
basis its Gram matrix `G` has

```text
signature       (1,10)
determinant     2048 = 2^11
Smith diagonal  (2,2,2,2,2,2,2,2,2,2,2).
```

Every entry of `G` is even.  The half-scaled form `H=G/2` is integral, odd,
unimodular, and Lorentzian.  Hence it is the unique odd unimodular lattice
`I_(1,10)`.  This is not only a classification argument: the following
determinant-minus-one matrix is an explicit isometry.  Its columns are the
standard `I_(1,10)` basis in the reduced fixed basis, and it satisfies

```text
Q^T (G/2) Q = diag(1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1).
```

```text
Q =
[ 24, 1, 0,  0, 1, 1, 0, -3,  0, -6, -6]
[ 28, 0, 1, -1, 0, 0, 0,  0,  5, -7, -7]
[-21, 0, 0,  1, 0, 0, 0,  0, -4,  5,  5]
[ 24, 0, 0,  0, 1, 1, 0, -3,  0, -6, -6]
[  6, 0, 0,  0, 0, 0, 0, -1,  0, -2, -1]
[ 24, 0, 0,  0, 0, 1, 0, -3,  0, -6, -6]
[  4, 0, 0,  0, 0, 0, 1,  0,  0, -1, -1]
[ -5, 0, 0,  0, 0, 0, 0,  0, -1,  1,  1]
[ -8, 0, 0,  0, 0, 0, 0,  0, -1,  2,  2]
[ -4, 0, 0,  0, 0, 0, 0,  0, -1,  1,  1]
[ -8, 0, 0,  0, 0, 0, 0,  1,  0,  2,  2].
```

Multiplying a quadratic form by 2 changes neither its integral orthogonal
group nor the reflection formula.  Therefore the fixed lattice and
`I_(1,10)` have the same integral reflection group.

## Primitive mirrors and reflection integrality

Each retained facet is initially represented by a sum of Shimada dual-wall
vectors.  The certificate converts that covector into fixed-lattice
coordinates and replaces it by the primitive integral vector defining the
same hyperplane.  The result is

| lattice normalization | square/divisibility | number |
|---|---:|---:|
| original fixed lattice `G` | `(-2,2)` | 18 |
| original fixed lattice `G` | `(-4,2)` | 12 |
| half-scaled `I_(1,10)` | `(-1,1)` | 18 |
| half-scaled `I_(1,10)` | `(-2,1)` | 12 |

For every primitive root `r`, the matrix

```text
s_r = I - 2 r (r^T G)/(r^T G r)
```

is integral, involutive, has determinant `-1`, and preserves `G`.  In
addition, the Smith diagonal of the `11 x 30` target root matrix is

```text
(1,1,1,1,1,1,1,1,1,1,1).
```

So the roots span the whole fixed lattice integrally; there is no hidden
root-span index or reflection-stable overlattice correction.

## Exact identification of the maximal chamber

Use the standard basis `e0,e1,...,e10` of `I_(1,10)`, with `e0^2=1` and
`ei^2=-1` for `i>0`.  A 12-root system for `P10` is

```text
r_i  = e_i - e_(i+1),                   1 <= i <= 9,
r_10 = e_10,
r_11 = e_0 - e_1 - e_2 - e_3,
r_12 = 3e_0 - e_1 - ... - e_10.
```

In `2*I_(1,10)` these consist of ten roots of type `(-4,2)` and two of type
`(-2,2)`.  All twelve reflections are integral.  Their normalized Coxeter
relations match the pinned CoxIter record

```text
graphs/group units qf/10-1,1,1,1,1,1,1,1,1,1,1.coxiter
```

entry by entry in the displayed root order; the facet permutation is the
identity.

There are two further independent checks.

- Alexeev's exact matrix for the 2-elementary lattice `(11,11,1)` is the same
  weighted root Gram matrix.  His source identifies the lattice as
  `I_(1,10)(2)`, gives 12 roots and two isotropic orbits.  The one-based map
  from the roots above to Alexeev's root order is

  ```text
  (3,4,5,6,7,8,9,11,1,12,2,10).
  ```

- Scharlau--Kirschmer `List10/Gram10` entry 1 is `I_(10,1)()` in their
  opposite-sign convention.  Its finite face counts, facets through finite
  vertices, are

  ```text
  (12,65,210,450,672,714,540,285,99,17),
  ```

  and it has two cusps.  The signed lattice Gram matrix matches entrywise
  after moving the time coordinate from first to last.

These data identify the ambient full root-reflection chamber as the classical
Vinberg `P10`.  They do not turn the 30-wall chamber into a new maximal group.

Sources:

- [Vinberg, groups of units of quadratic forms](https://doi.org/10.1070/SM1972v016n01ABEH001346)
- [Alexeev, reflective hyperbolic 2-elementary lattices](https://arxiv.org/abs/2209.09110)
- [Scharlau--Kirschmer reflective Lorentzian census](https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html)
- [CoxIter repository](https://github.com/rgugl/CoxIter)

## Solomon index and covolume

The connected elliptic component census has 861 components, with rank
histogram

```text
rank:   1  2  3   4   5   6   7   8  9
count: 30 45 90 126 114 120 108 132 96.
```

Combining mutually orthogonal components gives every finite special
subgroup.  The finite-special subset counts by rank are

| rank | `W30` | `W(P10)` |
|---:|---:|---:|
| 0 | 1 | 1 |
| 1 | 30 | 12 |
| 2 | 321 | 65 |
| 3 | 1,890 | 210 |
| 4 | 7,072 | 450 |
| 5 | 17,784 | 672 |
| 6 | 30,452 | 714 |
| 7 | 34,857 | 540 |
| 8 | 25,206 | 285 |
| 9 | 10,030 | 99 |
| 10 | 1,333 | 17 |

Solomon's formula

```text
chi(W) = sum_{T finite special} (-1)^|T| / |W_T|
```

then gives

```text
chi(W30)  = -527/483840,
chi(P10)  = -1/11147673600,
index     = chi(W30)/chi(P10) = 12,142,080.
```

The quotient is an index because the explicit primitive-root coordinates
first prove `W30 < W(P10)`; it is not an index inferred from an unexplained
integer Euler ratio.

Gauss--Bonnet in dimension 10 gives

```text
vol(P30) = -(32*pi^5/945) chi(W30)
         = pi^5 * 527/14288400.
```

## Bounded diagram search and novelty boundary

The pinned ordinary CoxIter corpus has 1,214 records.  Its `H10` facet-count
histogram is

```text
11:14, 12:8, 13:4, 14:1.
```

The separate reflection-doubling archive has 2,959 records.  Its `H10`
histogram is

```text
12:3, 13:4, 14:3, 15:1, 16:3,
19:1, 20:3, 28:2, 44:2.
```

Neither source contains an `H10` 30-facet header, so there is no candidate on
which an exact weighted-diagram isomorphism could fire.  This proves absence
only in those two pinned corpora.

The correct status is therefore:

| question | answer |
|---|---|
| new hyperbolic Coxeter polytope found by this search? | yes, a certified positive 30-wall algebraic construction |
| new lattice? | no, `2*I_(1,10)` |
| new maximal reflection group? | no, Vinberg `P10` |
| new arithmetic commensurability class? | no |
| literal 30-wall chamber printed in the bounded CoxIter corpora? | no |
| worldwide priority for the literal chamber settled? | no |
| remaining plausible novelty | divisor-involution realization and possibly the literal index-12,142,080 subgroup chamber |

## Replay

From the repository root:

```bash
PYTHONPATH=src:. python3 scripts/certify_shimada_disc16_d5_h10_lineage.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_shimada_disc16_d5_h10_lineage.py
```

The first command recomputes the Smith kernel, explicit isometry checks, all
30 primitive mirrors, all 861 connected elliptic components, both Solomon
sums, the three-source `P10` identity, and both corpus scans.

