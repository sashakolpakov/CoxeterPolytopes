# Hostile audit: source 29, triple `123`

## Verdict

**PASS, with the original scope unchanged.**  An independent enumerator with
no imports from the project mathematics reproduces the triple-`123` fibre
model, strict budgets, symmetry quotients, mask-`230` cut, total `144`, and
the near-seed defect.  It finds no arithmetic or stabilizer error.  The 144
objects are pair-configuration branches, **not** strict-link survivors;
source 29 remains **OPEN**.

## Independent reconstruction

The audit starts only from the frozen source data

```text
missing edges:  0-13, 1-13
missing A4s:    0123, 0124
T:              123
Q:              4,...,12
```

It directly counts the link extensions of `T` as `(1,9,36,84)`.  The
six-dimensional link transform gives

```text
h=(1,3,6,10),  g=(2,3,4).
```

The minimal-ground A5 census is 1,600.  Exactly 36 contain `T`, bijectively
labelled by pairs in `C(Q,2)`.  The 84 six-fibres `T union U` each have
exactly three non-`T` killers, and their union has size 252.

For `d=0,1,2`, the Macaulay inequalities independently give the same
missing-six intervals and residual support budgets:

| pair configuration | base coverage | allowed `m` | new-fibre budget |
|---|---:|---:|---:|
| none | 0 | `0..4` | 4 |
| one pair | 7 | `7..9` | 2 |
| two intersecting pairs | 13 | `13..14` | 1 |
| two disjoint pairs | 14 | `13..14` | 0 |

The raw two-pair census is 252 intersecting and 378 disjoint configurations.

## Symmetry and the fixed masks

Burnside averaging, separate from the project's generator-orbit routine,
gives:

```text
zero support, S8 on 5..12:          d=0,1,2 -> 1,2,5
active support, S3(567)xS5(8..12):  d=0,1,2 -> 1,5,24
```

The residual groups are correct.  Once `T=123` is fixed, neither the
`S_{0,1}` nor the `S_{3,4}` factor preserves `T`.  The active mask patterns
are each invariant under `S3(567)xS5(8..12)`, so no pattern-specific action
on the outside pairs was silently lost.

Among the 24 active two-pair orbits, the independent profile is

| pair intersection | covers `567` | orbits |
|---:|:---:|---:|
| 0 | no | 7 |
| 0 | yes | 4 |
| 1 | no | 10 |
| 1 | yes | 3 |

Of masks `227,229,230`, only `230={1,2,5,6,7}` is a `T=123` killer on the
fibre `567`.  It occurs in exactly the last two normalized active patterns.
Their disjoint-pair budget is zero, so precisely the seven disjoint orbits
not covering `567` are removed from each.  Thus the exact configuration
count is

```text
8 + 30 + 30 + 30 + (30-7) + (30-7) = 144.
```

This accounts for fixed masks before quotienting the still-undecided
outside killer-fibre supports by each configuration's residual stabilizer.
That later quotient is not part of the 144 count and remains to be solved.

## Near seed

Direct mask decoding finds the two `T`-blocks `2094` and `302`, labelled by
the intersecting pairs `{5,11}` and `{5,8}`.  They cover 13 fibres.  Among
the other 19 masks, exactly two lie on new `T`-fibres:

```text
342  -> {4,6,8}
1578 -> {5,9,10}
```

Consequently `m=13+2=15`, while `d=2` gives `g2=1` and
`g3=4+5*2-15=-1`.  The lower-bound defect is exactly one.

Artifacts:

- [`independent_h9_14_record0_a2_source29_triple123_hostile_audit.py`](../scripts/independent_h9_14_record0_a2_source29_triple123_hostile_audit.py)
- [`test_h9_14_record0_a2_source29_triple123_hostile_audit.py`](../tests/test_h9_14_record0_a2_source29_triple123_hostile_audit.py)

The guarded replay peaks at about `0.026 GiB`; no timeout or unknown solver
status is used.
