# Silver `A1^24` / Golay H24 exact screen

This screen studies

\[
  \langle-a\rangle\oplus N(A_1^{24})
  \quad\text{over }\mathbf Z[\sqrt2]
\]

for the three admissible coefficients
`a = 1+sqrt(2), sqrt(2), 2sqrt(2)`.  The Niemeier lattice uses its full
extended-Golay Construction-A glue; no root-lattice-only shortcut is used.

## Complete mirror norm classes

- `a=1+sqrt(2)`: canonical norms `{1,2}`.
- `a=sqrt(2)`: canonical norms `{2,2+sqrt(2)}`.
- `a=2sqrt(2)`: canonical norms `{2,4+2sqrt(2)}`.

Write `pi=sqrt(2)`, so `(2)=(pi)^2`.  Local unimodularity away from `pi`
eliminates all other primes.  At `pi`, primitivity, evenness, and
`q(r)|2I(r)` leave valuations `{0,2}`, `{1,2}`, and `{2,3}`, respectively.
For the last coefficient, the valuation-three class additionally requires
the spatial vector to lie in `pi*N`; this condition removes square-equation
solutions which are not actual mirrors.  The unit `1+sqrt(2)` has sign
pattern `(+,-)`, and every totally positive unit is a square, giving the
listed canonical generators.

## First complete packets

Each coefficient has exactly one height-zero root orbit, giving the 24
`A1` chamber walls.  Its first nonempty positive packet consists of exactly
24 coordinate caps.  The dominant affine section is a 24-cube, so its
`2^24=16,777,216` vertices are the complete projective cone-ray census.

| `a` | first height | cap coordinate | finite | cusp | hyperideal |
|---|---:|---:|---:|---:|---:|
| `1+sqrt(2)` | `2+2sqrt(2)` | `1` | 25 | 0 | 16,777,191 |
| `sqrt(2)` | `1+sqrt(2)` | `sqrt(2)/2` | 25 | 0 | 16,777,191 |
| `2sqrt(2)` | `1+sqrt(2)` | `1` | 25 | 0 | 16,777,191 |

Thus none of these first packets has finite volume.

## Unit coefficient: second complete packet

Through exact height 20 the only nonempty packets for `a=1+sqrt(2)` are:

1. the 24 coordinate caps above;
2. at `x=2+sqrt(2)`, height `7+5sqrt(2)`, the multiset
   `[(3+sqrt(2))/2]^1 [ (1+sqrt(2))/2]^7 [0]^16`.

The latter is labelled by an octad and a distinguished point of that octad,
so it contains exactly `759*8=6,072` walls.  Exact Golay intersection cases
show that all same-height pairings are nonpositive (one same-octad case is
tangent).  Their coefficient rows have rank 24.  The resulting symmetric
extreme ray has coordinate

\[
  c=\frac27+\frac1{14}\sqrt2
\]

and norm

\[
  \frac{167}{49}+\frac{47}{49}\sqrt2>0.
\]

It is therefore still hyperideal.

## Exact relaxed H=100 obstruction

For every complete control shell through height 100, an exact half-integral
square-partition DP was filtered by the joint Golay residue table and the
first coordinate caps.  Later-wall acuteness was deliberately omitted,
which enlarges the candidate set.  Among 210 resulting unlabelled patterns,
the maximum symmetric slope is

\[
  \max \frac{\sum_i d_i}{x}=7-\frac32\sqrt2
\]

at height `41+29sqrt(2)`.  Its square is strictly below the finite/cusp
threshold `12(1+sqrt(2))`.  Consequently

\[
  \left(1;\left(\frac{10}{89}+\frac{17}{178}\sqrt2\right)^{24}\right)
\]

lies in the actual H=100 prefix cone and has norm

\[
  \frac{3815}{7921}+\frac{239}{7921}\sqrt2>0.
\]

This proves that the exact unit-coefficient prefix is not finite-volume
through H=100.  It is not a global nonreflectivity claim.

## Reproduction

Run:

```bash
PYTHONPATH=src python3 scripts/search_silver_a1_golay_h24.py
PYTHONPATH=src pytest -q tests/test_silver_a1_golay.py tests/test_silver_d24_plus.py
```

Artifact:
`data/enumerations/silver_a1_golay_h24_h100.json`

SHA-256:
`b1affb892cda04db90f39e9ce99ff3eb15b2e00dcf715437769c5a7f1c8d668e`

