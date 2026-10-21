# Record 14: complete compact `a=2` schedule no-go

## Theorem and exact scope

No compact hyperbolic Coxeter 9-polytope with fourteen facets can have the
Record-14 missing-edge graph and exactly two minimal missing four-faces.

Together with the separate upper-schedule certificate, schedules
`a=2,3,4,5,6` are now closed.  This is not yet a whole-shadow theorem:
`a=0,1` remain open.  It is also compact-only and says nothing about a
finite-volume noncompact chamber, whose ideal links and parabolic diagrams
would require a different Vinberg audit.

## Complete source census

The frozen positive-`A4` quotient contains exactly 36 full graph-action
orbits in schedule `a=2`.  Every literal key, orbit size, and representative
is consumed in archived order.  For each representative, the program derives
the full forced A5 shadow, the remaining minimal-A5 candidates, and all local
addition intervals.

The exact incidence model includes:

- four or five or three additional A5 rows, as dictated by the forced shadow;
- destruction of every base transversal of order at most four;
- an equivalence variable for each surviving order-five facet complement;
- exactly 106 surviving facet complements; and
- ridge degree zero or two in every six-set.

The result is

```text
32 source orbits -> 0 literal targets
 4 source orbits -> 6 literal targets each
total             24 literal targets
```

The four target-bearing source indices are exactly `1,5,13,16`.

## Independent literal replay

The production enumeration uses Z3.  A separate source file rebuilds the
graph, local Macaulay arithmetic, candidate families, facet equivalences, and
ridge equations without importing the production completion engine, then
enumerates with OR-Tools CP-SAT.  Every terminal model is independently
blocker-dualized.

The two solvers agree on the complete literal target set for each of all 36
source orbits—not merely on the counts.  Both give exactly the same 24 pure
106-facet families.

## Exact target-orbit quotient

No target is discarded by an assumed symmetry.  For each target-bearing
source representative, all `2*9!` graph actions are scanned and only those
preserving its two A4 rows are retained.  Their orders and target partitions
are

| source orbit | fixed-A4 action order | literal target groups |
|---:|---:|:---|
| 1  | 432 | `(0,1,2,3,4,5)` |
| 5  | 432 | `(0,1,2,3,4,5)` |
| 13 | 288 | `(0,1,2,3,4,5)` |
| 16 | 288 | `(0,1,2)`, `(3,4)`, `(5)` |

Thus the 24 literal families form exactly six fixed-source target orbits.

## Unrestricted compact-label gate

For one representative of each of the six target orbits, the verifier joins
all 108 labelled order-four Lanner placements and all 420 labelled order-five
Lanner placements on every minimal nonface.  After every join it tests every
actual rank-nine vertex restriction for ellipticity.  The successive sparse
state counts are

```text
source 1 / target 0:  108,84,48,216,0
source 5 / target 0:  108,84,0
source 13 / target 0: 108,84,12,24,0
source 16 / target 0: 108,84,12,0
source 16 / target 3: 108,0
source 16 / target 5: 108,0
```

All six terminate at zero.  This is unrestricted over the compact finite
alphabet: unassigned meeting labels begin at two, and increasing any to
`3,4,5` increases the nonnegative cosine adjacency matrix entrywise, so a
failed elliptic principal restriction cannot be repaired by later labels.

## Frozen artifacts

- `scripts/certify_h9_14_record14_a2_schedule_no_go.py`
- `scripts/explore_h9_14_record14_fixed_orbit_z3.py`
- `scripts/independent_h9_14_record14_a2_audit.py`
- `data/enumerations/h9_14_record14_a2_schedule_no_go.json`
- `tests/test_h9_14_record14_a2_schedule_no_go.py`

Certificate digest:

```text
d47d7c062ae1b40a3dca8bb1f3ffe6c3f87b5444e5c0d78a8515e2caa4aab7a1
```

Reproduction:

```console
PYTHONPATH=src:. python3 scripts/certify_h9_14_record14_a2_schedule_no_go.py
PYTHONPATH=src:. caffeinate -i python3 \
  scripts/certify_h9_14_record14_a2_schedule_no_go.py --fresh \
  --write data/enumerations/h9_14_record14_a2_schedule_no_go.json
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_h9_14_record14_a2_schedule_no_go.py
```
