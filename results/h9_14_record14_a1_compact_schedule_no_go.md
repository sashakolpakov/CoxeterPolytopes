# Record 14: complete compact `a=1` schedule no-go

## Theorem and scope

No compact hyperbolic Coxeter 9-polytope with fourteen facets can have the
Record-14 missing-edge graph and exactly one minimal missing four-face.

Together with the previously frozen `a=2` and upper-schedule certificates,
this closes schedules `a=1,2,3,4,5,6`.  Schedule `a=0` is not part of this
certificate.  The theorem is compact-only: it makes no inference about
finite-volume noncompact chambers or ideal links.

## Complete source census

The exact positive-`A4` quotient has five full missing-graph orbits in this
schedule.  Two independently implemented Boolean models consume all five:

- production uses Z3;
- the replay rebuilds the graph, Macaulay bounds, candidate rows, blocker
  equivalences, facet count, ridge equations, and symmetry slice in CP-SAT.

Both models impose exact addition cardinality and local intervals, destroy
every base transversal of order at most four, retain exactly 108 order-five
facet complements, and require degree zero or two in every six-set.  Each
also imposes only generator lex inequalities.  This is lossless because the
global lexicographic minimum of every finite stabilizer orbit satisfies all
generator inequalities.

The two solvers terminate and agree on the literal generator-lex target sets:

| source orbit | generator-lex pure targets | full target orbits |
|---:|---:|---:|
| 0 | 11 | 6 |
| 1 | 0 | 0 |
| 2 | 0 | 0 |
| 3 | 8 | 4 |
| 4 | 0 | 0 |

Thus the coarser generator slice has 19 targets.  Exact quotienting by the
full stabilizer of the fixed `A4` row gives ten target orbits.

## Exact full-symmetry quotient

The Record-14 graph fixes vertices `2,12,13`, permits the twin swap
`0 <-> 1`, and has a symmetric isolate block `3,...,11`.  For each source,
the verifier enumerates every permutation inside and outside its fixed
four-row and includes the twin swap exactly when it preserves that row.  The
two target-bearing stabilizers have orders 4,320 and 10,080.  Their exact
partitions of the generator-lex targets have sizes

```text
source 0: (2,2,2,2,2,1)
source 3: (2,2,2,2)
```

No target is discarded by a guessed symmetry.

## Unrestricted compact-label gate

For one representative of each of the ten full target orbits, the verifier
joins all 108 labelled order-four Lanner placements and all 420 labelled
order-five placements on every minimal nonface.  After every join it checks
all 108 actual rank-nine vertices for ellipticity.  Every target terminates
with zero sparse states.

This is unrestricted over finite labels `2,3,4,5`.  Unassigned meeting pairs
start right-angled.  Increasing an entry to `3`, `4`, or `5` increases the
nonnegative cosine adjacency matrix entrywise, so Perron--Frobenius
monotonicity cannot repair an already nonelliptic actual-vertex restriction.
No Gram-signature stage is needed after the state set becomes empty.

## Frozen artifacts

- `scripts/certify_h9_14_record14_a1_schedule_no_go.py`
- `scripts/explore_h9_14_record14_a1_symmetry_z3.py`
- `scripts/independent_h9_14_record14_a1_audit.py`
- `data/enumerations/h9_14_record14_a1_schedule_no_go.json`
- `tests/test_h9_14_record14_a1_schedule_no_go.py`

The certificate digest is

```text
7c579a3ac2260275a504d7a5a2fa2cf851cb2eb379c040069e130a07cc1da4f3
```

Reproduce the frozen proof with

```console
PYTHONPATH=src:. caffeinate -i python3 \
  scripts/certify_h9_14_record14_a1_schedule_no_go.py --fresh \
  --write data/enumerations/h9_14_record14_a1_schedule_no_go.json
PYTHONPATH=src:. python3 \
  scripts/certify_h9_14_record14_a1_schedule_no_go.py
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_h9_14_record14_a1_schedule_no_go.py
```
