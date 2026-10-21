# Record 0, compact `a=2`, source 29: the triple-123 fibre budget

## Outcome

Source 29 remains **OPEN**.  The strict link at `123` nevertheless gives a
second exact analytic cut after the
[`012` cut](h9_14_record0_a2_source29_triple012_cut.md): its six normalized
branches refine into 144 finite outside-pair configurations, each with an
explicit budget on the remaining killer fibres.  No timeout or unsuccessful
search is used in this result.

## Exact fibre model

Let `T=123`.  The nine possible outside vertices are

```text
Q = {4,5,6,7,8,9,10,11,12}.
```

The 84 allowed six-sets through `T` are `T union U`, with
`U in C(Q,3)`.  There are 36 candidate A5s containing `T`, in bijection with
outside pairs `A in C(Q,2)`: the block is `T union A`.  Every fibre `T union
U` also has exactly three non-`T` killer blocks, `(T-{t}) union U` for
`t in T`; these give 252 distinct candidates.

The pre-addition link vector is `(g1,g2,g3)=(2,3,4)`.  If `d` selected A5s
contain `T` and `m` six-fibres are missing, then

```text
g2 = 3-d,
g3 = 4+5d-m.
```

Strictness gives the complete table:

| `d` | selected outside pairs | fibres already covered | allowed `m` | outside-fibre budget |
|---:|---|---:|---:|---:|
| 0 | none | 0 | 0..4 | 4 |
| 1 | one pair | 7 | 7..9 | 2 |
| 2 | two intersecting pairs | 13 | 13..14 | 1 |
| 2 | two disjoint pairs | 14 | 13..14 | 0 |

The last column is decisive: a non-`T` killer supported on a fibre already
covered by a selected outside pair costs nothing, while its support on a new
fibre consumes exactly one unit.

The pair census is also exact.  Among the 630 unordered choices of two
outside pairs, 252 intersect and 378 are disjoint.

## Symmetry refinement of all six `012` branches

In the zero-support `012` branch, the stabilizer relevant to `T` acts as
`S8` on `5,...,12` with vertex `4` fixed.  The cases `d=0,1,2` have
respectively `1,2,5` pair-configuration orbits: 8 refined branches.

In every active `012` branch, the normalized free fibre is `567`, and the
relevant action is `S3 x S5` on `567` and `8,...,12`, with `4` fixed.  The
cases `d=0,1,2` have `1,5,24` orbits: 30 refined branches before accounting
for a fixed `123` killer.

The five active fibre patterns are

```text
{227}, {229}, {227,229}, {229,230}, {227,229,230}.
```

Mask `230` is the `123` killer `(12) union 567`.  In the final two patterns
it forces fibre `567` to be missing.  Seven disjoint-`d=2` orbits do not
cover `567`; their outside-fibre budget is zero, so they are impossible.
Those two active patterns therefore have 23 rather than 30 refinements.
Altogether:

```text
8 + 3*30 + 2*23 = 144.
```

This is an exhaustive configuration refinement, not 144 claimed survivors.
Each row still needs its bounded killer-fibre support and global shadows
solved.

## Why the frozen near seed misses by one

The one-defect seed has the two `T`-blocks

```text
2094 = 123 union {5,11},
302  = 123 union {5,8}.
```

Their outside pairs intersect, so they cover 13 fibres and allow at most one
new fibre.  The seed selects two outside killers:

```text
1578 on U={5,9,10},
342  on U={4,6,8}.
```

Thus `m=13+2=15`, while strictness requires `m<=14`; this is exactly its
single `g3=-1` defect.  If the two `T`-blocks are held fixed, at least one of
these masks must be removed or moved onto a base-covered fibre.

## Next finite completion step

For each of the 144 pair configurations, quotient the support of at most
`4,2,1,0` outside killer fibres by its residual stabilizer and impose the
exact `(M6,M7,M8)=(147,420,629)` shadows plus all strict links.  Only a strict
seed should proceed to blocker purity, 176 facets, and ridge degree `0/2`.

Artifacts:

- [archive](../data/enumerations/h9_14_record0_a2_source29_triple123_cut.json)
- [certifier](../scripts/certify_h9_14_record0_a2_source29_triple123_cut.py)
- [tests](../tests/test_h9_14_record0_a2_source29_triple123_cut.py)

No source-29 exclusion, compact incidence witness, Coxeter labeling,
Lorentzian Gram matrix, or hyperbolic polytope is asserted.
