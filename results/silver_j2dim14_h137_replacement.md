# Silver J2dim14 H28 continuation through H=137+97sqrt(2)

## Outcome

The old strict hyperideal ray `z=(1,(3/2)h)`, with `h=G^(-1)e_13`,
survives the complete unit-coefficient chronology through `H=103+73sqrt(2)`.
Its first genuine accepted cutter occurs at

\[
H=137+97\sqrt2.
\]

The cutter does not close the chamber: `z'=(1,(5/4)h)` is strictly
hyperideal and survives every possible wall through that height. Thus this
bounded prefix remains infinite-volume and supplies no compact example.

## Exact short-shell and CVP gate

PARI enumerates the complete norm-eight and norm-ten shells and gives
`|d_13|<=2` on each. For norms twelve and fourteen, Cauchy permits at most
four. An exact rational-LDL affine-slice enumeration visits respectively
331
and 9647
nodes. It finds exactly two norm-twelve vectors with coordinate `+4`, but no
norm-fourteen vector with coordinate `+4`. Hence the certified bounds are

| norm | 8 | 10 | 12 | 14 |
|:---|---:|---:|---:|---:|
| maximum `|d_13|` | 2 | 2 | 4 | 3 |

These bounds are imposed on every short combination of the binary Gram
vectors, together with exact positive semidefiniteness of the augmented Gram.
The resulting census is a necessary superset of all actual walls.

## Complete chronology

There are 40 arithmetic controls at
27 heights through the endpoint. Exactly
4 controls retain 9
necessary binary types.

| height | binary Gram `(u^2,v^2,u.v)` | max `h.(u+sqrt(2)v)` | margin at `3h/2` | margin at `5h/4` |
|:---|:---|:---|:---|:---|
| `14 + 10*sqrt(2)` | `[10, 10, 10]` | `2 + 2*sqrt(2)` | `3 + 1*sqrt(2)` | `7/2 + 3/2*sqrt(2)` |
| `82 + 58*sqrt(2)` | `[86, 40, 58]` | `7 + 5*sqrt(2)` | `7/2 + 5/2*sqrt(2)` | `21/4 + 15/4*sqrt(2)` |
| `82 + 58*sqrt(2)` | `[82, 42, 58]` | `8 + 6*sqrt(2)` | `2 + 1*sqrt(2)` | `4 + 5/2*sqrt(2)` |
| `103 + 73*sqrt(2)` | `[104, 52, 73]` | `8 + 6*sqrt(2)` | `4 + 2*sqrt(2)` | `6 + 7/2*sqrt(2)` |
| `103 + 73*sqrt(2)` | `[100, 54, 73]` | `8 + 6*sqrt(2)` | `4 + 2*sqrt(2)` | `6 + 7/2*sqrt(2)` |
| `137 + 97*sqrt(2)` | `[140, 68, 97]` | `11 + 8*sqrt(2)` | `3/2 + 1*sqrt(2)` | `17/4 + 3*sqrt(2)` |
| `137 + 97*sqrt(2)` | `[136, 70, 97]` | `14 + 10*sqrt(2)` | `-3 - 2*sqrt(2)` | `1/2 + 1/2*sqrt(2)` |
| `137 + 97*sqrt(2)` | `[132, 72, 97]` | `10 + 8*sqrt(2)` | `3 + 1*sqrt(2)` | `11/2 + 3*sqrt(2)` |
| `137 + 97*sqrt(2)` | `[128, 74, 97]` | `9 + 7*sqrt(2)` | `9/2 + 5/2*sqrt(2)` | `27/4 + 17/4*sqrt(2)` |

At `H=82+58sqrt(2)` and `H=103+73sqrt(2)`, the difference `u-v` has
norm eight or ten. The shell bound eliminates the raw Cauchy threats and
leaves strictly positive old-ray margins. The first surviving threat is the
Gram type `(136,70,97)` at `H=137+97sqrt(2)`.

## Explicit accepted cutter

Take the integral vectors `u,v` frozen in the certificate. They satisfy

\[
u^2=136,\quad v^2=70,\quad u\mathbin{\cdot}v=97,
\quad (u-v)^2=12,\quad (u_{13},v_{13})=(14,10).
\]

Then

\[
r=(8+5\sqrt2,u+\sqrt2v),\quad r^2=2.
\]

The root is primitive and reflection-integral. It cuts the old ray by
`3 + 2*sqrt(2)>0`.

It is nevertheless chronologically accepted. Exhaustion of all
374,976 first-packet walls gives maximum spatial
pairing `27 + 20*sqrt(2)`, below the acute
threshold `36 + 26*sqrt(2)` by the strict
margin `9 + 6*sqrt(2)`.

## Replacement ray

For `z'=(1,(5/4)h)`,

\[
(z')^2=5/2-(1+\sqrt2)=3/2-\sqrt2>0.
\]

The normalized first-packet margin is
`-1/2 + 2*sqrt(2)`.
The worst margin over every necessary wall type through `H=137+97sqrt(2)` is
`1/2 + 1/2*sqrt(2)`. Both are strictly positive,
so an open hyperideal patch remains.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/audit_silver_j2dim14_h137_replacement.py --verify
PYTHONPATH=src pytest -q tests/test_silver_j2dim14_h137_replacement.py
```
