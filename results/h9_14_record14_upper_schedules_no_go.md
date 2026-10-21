# Record 14: exact compact closure of schedules `a=3,4,5,6`

## Theorem and scope

For the `H9/14` Record-14 missing-edge graph

```text
M???????????o?w??
```

no compact Coxeter nerve can occur in Dehn--Sommerville schedules
`a=3,4,5,6`, where `a` is the number of minimal missing four-faces.

This is a genuine schedule-level theorem, but not a whole-shadow theorem.
Schedules `a=0,1,2` remain outside this certificate.  The result is
compact-only: it neither constructs nor excludes a finite-volume noncompact
chamber, because ideal links, parabolic cusp diagrams, and Vinberg endpoint
extensions are not modeled.

## Every schedule

Direct Dehn--Sommerville calculation from the five missing edges gives the
complete list

| `a` | all missing independent A5s | facets | sum of vertex-link `g4` |
|---:|---:|---:|---:|
| 0 | 8  | 110 | 42 |
| 1 | 14 | 108 | 36 |
| 2 | 20 | 106 | 30 |
| 3 | 26 | 104 | 24 |
| 4 | 32 | 102 | 18 |
| 5 | 38 | 100 | 12 |
| 6 | 44 | 98  | 6  |

Every compact vertex link is spherical and has strict `g4 >= 1`; the sum
over fourteen links is therefore at least fourteen.  This excludes `a=5,6`.
Global `g4` is nonpositive for `a>=7`, so the displayed table is complete.

## Literal source-orbit coverage

The independently frozen positive-`A4` quotient supplies exactly 59 source
orbits at `a=3` and 13 at `a=4`.  The certificate consumes all 72 literal
keys, sizes, and representatives in their archived order.  No representative
selection or additional symmetry assumption is introduced downstream.

For each source orbit, the production traversal derives the complete forced
A5 shadow and every remaining minimal-A5 candidate, enforces the exact local
addition intervals, and destroys every base transversal of order at most
four.  It then imposes the scheduled order-five facet count and blocker
purity.  Across all 72 branches it visits 6,566,617 search states.  There are
30,399 families after the small-transversal gate, but none has the scheduled
pure blocker.  Hence zero incidence target survives.

## Structurally independent replay

The second implementation does not replay the recursive cover traversal.  It
builds one Boolean model per source orbit with:

- exact addition cardinality and every vertex interval;
- one destruction clause for every order-at-most-four transversal;
- an iff literal for every surviving order-five transversal;
- the exact scheduled facet count; and
- degree zero or two in every six-set.

It then blocks and enumerates all solutions, with literal blocker
dualization of every model.  All 72 systems are UNSAT before blocking any
model: independent model total `0`, independent pure-target total `0`.
Thus the two implementations agree on the terminal empty target set while
using different state spaces and branching mechanisms.

## What remains

The theorem-facing boundary has moved from seven open schedules to three:

```text
closed:    a=3,4,5,6
open here: a=0,1,2
```

Residuals already found during the continuing `a=2` census are streamed to
the unrestricted compact-label gate separately.  They do not weaken this
upper-schedule theorem and are not silently included in it.

## Frozen artifacts

- `scripts/certify_h9_14_record14_upper_schedules.py`
- `scripts/explore_h9_14_record14_addition_cover.cpp`
- `scripts/explore_h9_14_record14_fixed_orbit_z3.py`
- `data/enumerations/h9_14_record14_upper_schedules_no_go.json`
- `tests/test_h9_14_record14_upper_schedules.py`

The frozen certificate digest is

```text
9392e15c3d1c446150febb215573a7ecbe1a0843a31e6ed541b0d9896bd7fa8c
```

Reproduction:

```console
PYTHONPATH=src:. python3 scripts/certify_h9_14_record14_upper_schedules.py
PYTHONPATH=src:. caffeinate -i python3 \
  scripts/certify_h9_14_record14_upper_schedules.py --fresh \
  --write data/enumerations/h9_14_record14_upper_schedules_no_go.json
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_h9_14_record14_upper_schedules.py
```
