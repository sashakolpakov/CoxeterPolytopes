# RC7 fraction-state and spatial-R1 audit

This archive pins AlVin commit
`d92006a842cb4668e7c89e8c8e2d30337e4b0f26` in image
`coxeter-alvin:rc7-spatial-r1-d92006a`, immutable image ID
`sha256:b995880fcd16c6940af225e0c6e5d8734be495f434d8e66640774a4714252896`.

The image carries two narrow RC7 enumerator repairs. Neither patch pads a
search range or weakens an exact acceptance filter.

## 1. Vinberg-fraction state restoration

`RCyclotomic7Integer_VFs::computeNextAlVinFractions()` applied the order-three
Galois operation `conjugate(2)` to `rciE` twice while computing its two
nonidentity embeddings, but restored only `rciMinusAlpha0`. It then formed
the fraction from `sigma^2(rciE)` and left the shared possible-norm entry
rotated across later batches. The first patch adds the mathematically missing
third `rciE->conjugate(2)`, returning the norm to its canonical embedding
before division and storage.

This repairs the hostile form `[-1,0,0],1,[-2,-5,-7]`: its normalized missing
wall has control `[0,1,-1]` and norm `[-6,-4,-4]`. The complete basis and
order-14 audit is in the preceding fraction-fix archive's `AUDIT.md`.

## 2. Recursive nonintegral-R1 spatial bound

For the p7 form, AlVin sorts the diagonal as

```text
[-1,0,0],[-3,-2,-2],1,1,1.
```

The required sixth wall is

```text
([0,1,-1], [0,-1,-1], 0, 0, 0),
```

with norm `[-6,-4,-4]` and height approximately `1.52445867`. In the power
basis with `theta=2*cos(2*pi/7)`, its relevant exact data are

```text
[0,-1,-1] = 1 + theta,
[-3,-2,-2] = 2 - theta,
[-2,-2,-3] = 1 + theta + theta^2,
[-1,-2,-2] = 2 + theta.
```

At recursion index one, the prepatch trace retained the correct current
coefficient and residual sum and computed the correct interval:

```text
q = [-3,-2,-2]
S = [-2,-2,-3]
R1 = sqrt(S/q) = 1+theta in
     [2.246979603717466, 2.246979603717469]
b range = -1..0
target prefix a=0,b=-1: c range = 0..0
```

Thus neither `q` nor `S` was mutated, conjugation order was already correct,
and GAOL enclosed the quantity it was asked to compute. The defect was the
quantity supplied to GAOL in the nonintegral branch: the source recomputed
`sqrt(S/[-1,-2,-2]) = sqrt(S/(2+theta))` instead of using the already
computed `R1=sqrt(S/q)` for the current diagonal coefficient.

For `a=0,b=-1`, the upper expression is exactly `theta-1`. The wrong lower
bound is approximately `theta-1-1.082`, whose ceiling is zero and therefore
excludes the required `c=-1`. The correct identity is

```text
(theta - 1) - R1 = (theta - 1) - (1 + theta) = -2.
```

The repaired lower endpoint is consequently `gaol_max-gaol_R1`; its ceiling
is `-2`, so `c=-1` is enumerated. This is precisely the interval analogue of
the adjacent integer-R1 branch, which subtracts `iR1`. All crystallographic,
principal-square, conjugate-square, component-order, and bilinear wall
filters remain unchanged and exact. A source-wide search found no other live
hardcoded `RCyclotomic7Integer(-1,-2,-2)` occurrence.

## Regression gates

All three mandatory gates use the pinned image above.

- `hostile-m14/result.json`: `acceptance_passes=true`; AlVin terminates with
  four walls, the exact order-14 Coxeter matrix, compactness, and a passing
  independent exact replay.
- `control/result.json`: `acceptance_passes=true`; the known eight-wall
  compact chamber is unchanged and passes independent exact replay.
- `p7-spatial-r1-prefix/result.json`: `acceptance_passes=true`; a deliberately
  capped `-maxv 6` regression matches all six exact roots in order, sees norm
  class `[-6,-4,-4]`, places the repaired wall at `e6`, and contains no
  `Unknown weight`. This bounded run makes no termination or volume claim.

The separate unbounded `ramified-7-rank5` capture crossed the former crash,
printed the repaired wall as `e6`, and continued through `e22` before an
intentional runner interrupt. Its `complete_termination=false` is expected;
it is evidence only for crash removal, not a discovery result.

The exact backend regression
`test_l7_spatial_r1_wall_is_enumerated_and_selected_exactly` proves in field
arithmetic that the basis-equivalent target has norm `2*(2-theta)`, height
`(-1+2*theta+theta^2)/2`, occurs in the exhaustive pool through height two,
and survives Vinberg selection. It passes independently. The nine existing
AlVin verifier regressions also pass.

No candidate, reflectivity, termination, finite-volume, or compactness claim
is made for the p7 seed beyond the exact bounded prefix.
