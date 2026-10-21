# Record 30: exact positive-schedule compact scout

Status: all positive missing-`A4` schedules have now passed through an exact
symmetry quotient.  Schedule `a=5` is impossible by a forced edge-link
obstruction; schedules `a=6,7` are impossible by the summed strict
vertex-link inequality.  Subsequent exact certificates close `a=1,2,3,4`;
only schedule `a=0` remains open on the compact Record-30 track.

This is compact-only.  It says nothing about finite-volume noncompact
chambers.

## Graph and schedules

The missing graph is `3K2 + 8K1`, with edges

```text
0--11, 1--12, 2--13.
```

There are 806 independent four-sets and 1,372 independent five-sets.  If
`a` is the number of missing independent four-sets, Dehn--Sommerville gives

```text
b = 8 + 6a,      facets = 132 - 2a.
```

The global `g4` is `8-a`, while the sum of the fourteen vertex-link `g4`
values is `6(8-a)`.  Strict positivity therefore excludes `a=6,7` because
the sums are only 12 and 6.

## Exact positive quotient

The local vertex/edge/triple interval table reduces the 806 independent
four-sets to 742 individual candidates.  The exact quotient uses
`Aut(3K2) x S8` together with permutations of the missing-`A4` rows.

| `a` | raw unordered families | interval-feasible | strict/aggregate | full orbits |
|---:|---:|---:|---:|---:|
| 1 | 742 | 742 | 742 | 3 |
| 2 | 274,911 | 270,375 | 97,699 | 21 |
| 3 | 67,811,380 | 64,513,204 | 524,552 | 32 |
| 4 | 12,528,152,455 | 11,245,592,071 | 28,896 | 7 |
| 5 | 1,849,155,302,358 | 1,510,284,730,038 | 672 | 1 |

The `a=5` traversal is an independent row/column quotient.  It visits
125,783,953 isolate-column paths, removes 37,660,554 equal-active-row
duplicates and 2,103,609 repeated-row paths, and leaves 85,285,134 exact
`S8` orbits on the interval ground set.  Only 24 `S8` orbits, representing
672 labelled families, pass the vertex aggregate gate; all form one full
graph orbit.

## Exact `a=5` obstruction

The unique orbit has forced missing-`A5` shadow 36, so its schedule would
require two additional minimal `A5`s.  Before either can be chosen, a
meeting-pair link has

```text
h = (1,3,4,4),       g = (2,1,0).
```

The top strict `g` entry is already zero.  An additional missing `A5`
containing this pair can only decrease it; one not containing the pair
leaves it zero.  Hence no member of the orbit can be the nerve of a compact
hyperbolic Coxeter polytope—or even pass the necessary compact link gate.

## Exact scope

This particular certificate closes schedules `a=5,6,7`; subsequent exact
certificates close `a=1,2,3,4` as well.  The archive deliberately keeps
`whole_shadow_compact_coxeter_excluded=false`, because `a=0` remains open.
The complete `a=2` theorem is recorded in
`results/h9_14_record30_a2_compact_schedule_no_go.md`, and the complete
`a=1` theorem in `results/h9_14_record30_a1_compact_schedule_no_go.md`.  No finite-volume
noncompact conclusion is made.

## Artifacts

- `scripts/explore_h9_14_record18_positive_orbits.cpp` (`-DH9_RECORD30`)
- `scripts/independent_h9_14_record30_a5.cpp`
- `scripts/certify_h9_14_record30_positive_scout.py`
- `data/enumerations/h9_14_record30_positive_scout.json`
- `tests/test_h9_14_record30_positive_scout.py`

Frozen certificate SHA-256:

```text
02a2238457b33d0cba59bb9fb85579a5e27c4520ea86f17a4ce0614fcabb111e
```
