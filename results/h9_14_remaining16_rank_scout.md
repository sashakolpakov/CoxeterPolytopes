# H9/14: scouting rank for the sixteen compact-facing shadows

## Subsequent status

Records 44, 42, 35, 43, and 32, the first five ranked targets below, have since
been frozen as whole-shadow compact-Coxeter no-gos; see
[`h9_14_record44_whole_shadow_compact_no_go.md`](h9_14_record44_whole_shadow_compact_no_go.md),
[`h9_14_record42_whole_shadow_compact_no_go.md`](h9_14_record42_whole_shadow_compact_no_go.md),
[`h9_14_record35_whole_shadow_compact_no_go.md`](h9_14_record35_whole_shadow_compact_no_go.md),
[`h9_14_record43_whole_shadow_compact_no_go.md`](h9_14_record43_whole_shadow_compact_no_go.md),
and
[`h9_14_record32_whole_shadow_compact_no_go.md`](h9_14_record32_whole_shadow_compact_no_go.md).
Record 18 has now also been frozen as a whole-shadow compact no-go; see
[`h9_14_record18_whole_shadow_compact_no_go.md`](h9_14_record18_whole_shadow_compact_no_go.md).
The current compact-facing queue therefore has ten open shadows, exactly
`0, 1, 3, 4, 5, 7, 9, 14, 30, 31`, and starts with records 14 and 30.  The
compact-facing ledger is `132 + 10`; ordinary incidence remains `123 + 19`.
The original sixteen-shadow counts are
retained below as the dated scouting input, not as the current theorem ledger.

## Scope and nonclaim

After the whole-shadow compact closures of records 17, 12, and 47, the
sixteen positive graph shadows still open are

```text
0, 1, 3, 4, 5, 7, 9, 14, 18, 30, 31, 32, 35, 42, 43, 44.
```

This note ranks theorem-search routes.  It is **not itself a certificate**.
The original record-42, record-35, and record-43 scouting counts below have
since been superseded by the linked production certificates and independent
no-import replays.  The dated route descriptions are retained to record how
the theorem searches were selected.

Every conclusion in this note is about the **compact** problem.  A
finite-volume noncompact chamber has ideal faces and must instead satisfy
Vinberg's literal elliptic/parabolic endpoint and full zero-neighborhood
tests.  None of the completion censuses below decides that problem.

## Search-priority ranking

The following is a work ranking, not a mathematical invariant.  `a-range`
lists the Dehn--Sommerville branches left after only the strict total
vertex-link bound `sum g4 >= 14`.  `a0=(b,f8)` gives the number of independent
five-nonfaces and facets in the `a=0` schedule.

| priority | record | `a-range` | `a0=(b,f8)` | frozen selected compact-label obstruction? | reason for placement |
|---:|---:|---:|---:|---|---|
| 1 | 44 | `0..1` | `(4,90)` | yes | the `a=1` local shadow gate is empty and the `a=0` target is unique |
| 2 | 42 | `0..2` | `(5,102)` | yes, at `a=1` | scouting closes `a=2`, makes `a=1` the selected product, and leaves one `a=0` orbit which dies at the local label gate |
| 3 | 35 | `0..2` | `(5,102)` | yes, at `a=0` | only one new `a=1` incidence orbit survives; it has a bounded disjoint-block label/inertia gate |
| 4 | 43 | `0..1` | `(4,90)` | yes, for one of two `a=0` orbits | shortest dual/example route: the second orbit has 9,216 locally valid label states and a genuine dotted coupling |
| 5 | 32 | `0..3` | `(6,114)` | yes, for displayed `a=0` and `a=2` nerves | four schedules and two already-obstructed calibration completions; likely a six-row census |
| 6 | 18 | `0..4` | `(8,110)` | yes, at `a=0` | five schedules, but larger symmetry and only 145 single-`A4` shadow survivors |
| 7 | 14 | `0..4` | `(8,110)` | yes, at `a=2` | same schedule size as 18 but 369 single-`A4` survivors and less symmetry |
| 8 | 30 | `0..5` | `(8,132)` | yes, at `a=2` | six branches, but active automorphism order 48 and a product calibration nerve |
| 9 | 31 | `0..5` | `(9,122)` | yes, at `a=0` | six branches and a selected no-go, but weaker graph symmetry |
| 10 | 9 | `0..7` | `(12,130)` | no | eight branches; the `a=1` gate has only three symmetry orbits |
| 11 | 7 | `0..7` | `(12,130)` | no | eight branches and seven surviving single-`A4` orbits |
| 12 | 5 | `0..8` | `(13,142)` | no | nine branches and no selected unrestricted label gate |
| 13 | 4 | `0..10` | `(16,150)` | no | eleven branches |
| 14 | 3 | `0..13` | `(20,170)` | no | fourteen branches, partly offset by active automorphism order 8 |
| 15 | 1 | `0..13` | `(22,160)` | no | fourteen branches and a larger ridge core than record 3 |
| 16 | 0 | `0..16` | `(26,180)` | no | seventeen surviving schedules and the largest ridge core, 880 |

For comparison, the frozen ridge-core sizes down this list are respectively
`190,226,226,190,262,310,310,326,346,430,430,466,550,670,760,880`.
The bottom half should not be attacked by unstructured facet fishing: it
needs a parametrized stellar/multiwedge classification, a new link identity,
or a diagram-first obstruction shared across many schedules.

## Route 1: record 44, unique completion

The four numerical schedules are

```text
(a,b,f8,sum g4) = (0,4,90,24), (1,10,88,18),
                    (2,16,86,12), (3,22,84,6).
```

Thus the total link sum removes `a=2,3`, not `a=1`.  The exact local route is

```text
a=1: 690 independent A4s -> 378 interval-feasible -> 0
     strict forced-shadow/Macaulay-feasible.
```

For `a=0`, the current exhaustive four-row calculation has 12 active-row
configurations, 2,782 paths, 2,214 exact mask clutters, and 502
`S4-row x S9-isolate` quotient states.  Exactly one state has 90 facets, and
it is the frozen selected record-44 nerve.  That nerve already has an
unrestricted compact-label no-go.

**Stopping criterion.**  Freeze the all-schedule census, independently replay
the `a=1` failure and the lossless `a=0` quotient, prove literal equality of
the unique target with the selected nerve, and chain to the frozen label
certificate.  No new Gram search is then needed.

## Route 2: record 42, one new `a=0` orbit

The five schedules have

```text
(a,b,f8,sum g4) =
(0,5,102,30), (1,11,100,24), (2,17,98,18),
(3,23,96,12), (4,29,94,6).
```

Strict total positivity removes `a=3,4`.  The original exact scout gave:

- `a=2`: `148,785` `A4` pairs, `131,355` family-interval feasible, and
  zero forced-shadow/Macaulay-feasible;
- `a=1`: 546 `A4`s, 126 strict-shadow feasible, one
  graph/isolate orbit.  For representative `(2,3,4,5)`, the forced shadow
  has ten five-nonfaces and the vertex equations force the unique additional
  minimal `A5`, `(6,7,8,9,10)`.  This is exactly the frozen record-42 product
  nerve, already compact-label impossible;
- `a=0`: 82 active configurations give 1,098,455 row/column paths.  After
  28,805 repeated-row rejects there are 1,069,650 valid paths.  Thirty-five
  paths have 102 five-transversals, giving 15 exact five-face families.  Five
  pass the full selection and link-homology audit, and all five lie in one
  `Aut(C5) x S9-isolate x S5-row` orbit.

A representative of the last orbit has five minimal nonfaces

```text
(0,3,4,5,6), (1,3,4,5,6), (2,3,4,5,6),
(2,7,8,9,10), (7,8,9,10,12).
```

Its blocker is pure with 102 members and it has 102 actual rank-nine
vertices.  Joining complete Lanner placement tables and checking every
actual vertex gives passing-state counts

```text
420 -> 96 -> 0
```

on the first three common-four supports.

**Stopping criterion.**  A production census and a genuinely separate
replay must reproduce the full path/target/orbit counts, the unique `a=1`
orbit, the empty `a=2` gate, and the `420 -> 96 -> 0` label transcript.  At
that point every compact branch is closed.

## Route 3: record 35, bounded cross-block gate closed

Record 35 has the same schedule triples surviving total positivity as record
42:

```text
(a,b,f8,sum g4) = (0,5,102,30), (1,11,100,24), (2,17,98,18).
```

For `a=0`, the now-frozen five-row census gives 375,582 paths, 7,239
repeated-row rejects, and 368,343 valid paths.  Twenty paths attain 102
facets, giving 18 exact families.  Five pass the full selection/homology
audit and form three full graph-row-isolate orbits.  Complete local label
joins kill all three orbits; their decisive passing counts are

```text
orbit 0: 420 -> 96 -> 0
orbit 1: 420 -> 96 -> 1,152 -> 2,304 -> 0
orbit 2: 420 -> 96 -> 0.
```

The nonzero-`a` shadow census is also small after quotienting:

- `a=2`: 119,805 raw pairs, 111,699 family-feasible, only 28
  forced-shadow/Macaulay-feasible, one orbit.  Its forced shadow has 15
  five-nonfaces and requires two additions; all 540 incidence-cap-feasible
  added pairs have a transversal of order at most four;
- `a=1` has two `A4` orbits.  The `(2,3,4,11)` orbit requires three added
  `A5`s: among 16,328,620 raw triples, 1,794,605 meet vertex caps and 840
  have no small transversal, but none has the scheduled 100 facets.  The
  `(2,3,4,5)` orbit requires two added `A5`s: 31,626 raw pairs reduce to
  5,626 cap-feasible, then one no-small-transversal pair, which is a valid
  100-facet incidence complex.

The last candidate has minimal higher nonfaces

```text
A4: (2,3,4,5)
A5: (6,7,8,9,10), (3,4,5,11,12).
```

The local Lanner/actual-vertex join leaves `35,280 = 84*420` states.  The
six-wall overlapping block and the disjoint five-wall Lanner block have no
dotted cross pair.  The frozen exhaustive gate tests all 30 cross pairs:
`1,058,400` one-`m=3` trials leave 5,040 states, with exactly one cross edge
in each.  Labels `m=4,5` leave no state.  The surviving `m=3` states fall
into two exact block types, with multiplicities 4,320 and 720, and their
ten-wall principal Gram blocks both have inertia `(8,2,0)`.  By interlacing,
neither can occur in the Gram matrix of a compact hyperbolic 9-polytope.

**Stopping criterion (achieved).**  The production census and a separate
no-import replay reproduce every schedule, incidence orbit, label join, and
exact inertia obstruction.  Thus record 35 is excluded across the whole
compact shadow; see
[`h9_14_record35_whole_shadow_compact_no_go.md`](h9_14_record35_whole_shadow_compact_no_go.md).

## Route 4: record 43 as the dual/Gram-candidate track

The schedules left by total positivity are `(a,b,f8,sum g4)=(0,4,90,24)`
and `(1,10,88,18)`.

For `a=1`, 28 four-faces form one orbit after the local shadow gate.  A
representative has eight forced five-nonfaces and needs two additions.  Of
586,986 raw pairs, only 540 satisfy the vertex incidence caps, and every one
has a transversal of order at most four.  Thus this branch appears empty.

For `a=0`, the existing exploratory quotient has 1,073 paths, 887 exact
clutters, 394 row/isolate quotient states, and 234 states after also using
the active graph automorphisms.  Exactly two inequivalent states have 90
facets, and both pass the full incidence and link-homology audit.  The first
is the frozen selected multiwedge nerve and already has an unrestricted
compact-label no-go.  The second has missing five-faces

```text
(1,3,4,5,6), (3,4,5,6,11),
(2,7,8,9,10), (7,8,9,10,12).
```

Its two overlapping-support blocks each leave 96 label states, so the local
gate leaves `96*96=9,216` states.  Across the two six-wall blocks there are
36 pairs, one of which, `(1,12)`, is dotted.  Of the
`9,216*35=322,560` finite one-`m=3` trials, 4,032 pass all actual-vertex
tests.  Unlike record 35, the unavoidable dotted coupling can change the
inertia of the combined block.  This is therefore the most credible current
route to a surviving Gram candidate, rather than the shortest expected
no-go.

**Stopping criterion.**  Enumerate the 4,032 one-edge survivors through
labels 4 and 5 and compatible multi-edge patterns, retaining the dotted
variable symbolically.  Either:

1. exact Sturm/minor arguments give negative index at least two for every
   dotted value greater than one, yielding a compact no-go; or
2. a state reaches the full fourteen-wall Gram equations with rank 10,
   inertia `(9,1,4)`, all six dotted entries strictly below `-1`, and all
   actual vertices elliptic.  That would be a genuine Gram candidate, not
   yet a polytope example; it would still require the full compact Vinberg
   endpoint/completeness audit and an exact independent replay.
