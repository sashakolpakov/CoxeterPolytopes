# Source 29 after `012/123`: exact residual-support quotient

## Outcome

Source 29 remains **OPEN**.  The 144 pair-configuration branches now have an
exact second finite layer: 25,407 pattern-weighted orbits of the remaining
new `T=123` killer-fibre supports.  This is a finite search reduction, not a
count of strict survivors and not an incidence, Coxeter, Gram, or hyperbolic
witness.

## What is quotiented

For each normalized pair configuration `A` at `T=123`, let `B(A)` be the
six-fibres already missing because they contain a selected outside pair.
The strict local budget says that the new support

```text
S subset C({4,...,12},3) - B(A)
```

has size at most `4,2,1,0` for `d=0`, `d=1`, intersecting `d=2`, and
disjoint `d=2`, respectively.  If mask `230` is fixed and `567` is not
already in `B(A)`, then `567` is required to lie in `S`.

The certifier computes the literal stabilizer of `A` inside `S8` (zero
support) or `S3(567) x S5(8,...,12)` (active support), then quotients all
bounded supports by that stabilizer.  Every row is checked twice:

1. literal orbit deletion on the bounded support sets;
2. an independent Burnside cycle-index coefficient count.

The two counts agree for every configuration and support size.

## Exact counts

The 144 pattern branches collapse to 61 quotient-row types: 8 zero-support
rows, 30 active rows without fixed `230`, and 23 active rows with fixed
`230`.  Their support-orbit totals are:

| `012` family | fixed `230` | `d=0` | `d=1` | `d=2` | subtotal |
|---|:---:|---:|---:|---:|---:|
| zero | no | 335 | 55 | 22 | 412 |
| active | no | 6,867 | 729 | 215 | 7,811 |
| active | yes | 588 | 139 | 54 | 781 |

The unweighted row-type total is `9,004`.  Restoring the three active mask
patterns without `230` and the two with `230` gives the actual finite layer

```text
412 + 3*7,811 + 2*781 = 25,407 support branches.
```

As useful cross-checks, the largest `d=0` strata split by total support size
as follows:

```text
zero:                 1, 2, 8, 43, 281       (sizes 0..4)
active, no 230:       1, 7, 59, 581, 6219    (sizes 0..4)
active, fixed 230:    1, 6, 53, 528          (sizes 1..4)
```

## Necessary gate and honest lack of a cut

For every support representative and every literal `012` mask pattern, the
certifier tries all choices of one non-`T` killer per newly supported fibre.
It requires all monotone vertex/edge/triple incidence caps and checks that
the resulting minimal order-6/7/8 shadows do not already exceed
`(147,420,629)`.  All 25,407 branches pass this necessary gate.  Thus this
gate is sound but supplies **zero additional exclusions**; the exact global
shadows and full triple `g3` conditions remain the next problem.

## Solver-polarity incident

An exploratory CP model initially passed **missing-six** indicators to a
strict-link formula that expects **surviving-six** indicators.  That polarity
error produced spurious instant `INFEASIBLE` statuses.  It was caught before
freezing or ledger promotion, corrected, and every such status is explicitly
discarded.  With the correct polarity, two 60-second exploratory passes
returned `UNKNOWN`; these are also zero mathematical evidence.  No solver
outcome is used in this result.

## Scope and next step

The 25,407 objects specify pair configurations and new killer-fibre supports,
not the choice of every killer block or every remaining A5.  The completion
step is to impose the exact `(M6,M7,M8)=(147,420,629)` shadows and all strict
links within these finite branches, retaining only terminal exact proofs or
literal replayed seeds.  Source 29 remains **OPEN** until that is done.

Artifacts:

- [certifier](../scripts/explore_h9_14_record0_a2_source29_triple123_residual_supports.py)
- [test](../tests/test_h9_14_record0_a2_source29_triple123_residual_supports.py)
- [independent first-layer hostile audit](h9_14_record0_a2_source29_triple123_hostile_audit.md)

The frozen replay used a `0.5 GiB` watchdog, completed in `37.49 s`, and
peaked at `0.346 GiB` (`0.348 GiB` on the preceding transcript run).
