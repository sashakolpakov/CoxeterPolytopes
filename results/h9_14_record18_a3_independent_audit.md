# Independent all-orbit audit of the record-18 `a=3` schedule

## Result and scope

All forty positive-\(A_4\) orbits in the record-18 `a=3` schedule are excluded
at a necessary incidence gate by a fresh no-project replay.  The independently
generated orbit list splits as

```text
29 orbits requiring two A5 additions,
10 orbits requiring three A5 additions,
 1 orbit  requiring four A5 additions.
```

The complete exact aggregate is

| quantity | count |
|---|---:|
| raw unordered addition families | 14,933,639,143 |
| vertex-interval-feasible families | 2,392,528,542 |
| families destroying every transversal of order at most four | 6,011 |
| pure blockers with the scheduled 104 facets | 0 |

Thus the `a=3` schedule is incidence-impossible.  This is compact-only and
partial: schedules `a=0,1,2` remain open here, so this is not a whole-shadow
Coxeter theorem and has no finite-volume noncompact implication.

## Independent implementation

The input representatives are reconstructed from a fresh run of the independent
incidence-matrix positive-orbit quotient; neither production representatives nor
production addition tables are read.  For each representative the replay derives
the complete forced five-shadow and all vertex Macaulay intervals.

Candidates are ordered by decreasing small-transversal coverage.  At every
search node the algorithm pivots on an uncovered transversal and distinguishes
the earliest selected candidate destroying it.  This enumerates unordered
families without the production driver's chosen-set memo.  Cap-feasible totals
are independently obtained from a degree-signature generating recurrence.  For
each surviving family the complete blocker is generated through the three
minimal vertex-cover cores of the double-star missing-edge graph and checked for
five-uniformity and the scheduled cardinality.

The unique four-addition orbit visits 23,273,233 nodes and produces 4,626 of the
6,011 no-small-transversal families.  Every one fails the final pure-104 blocker
test.  The per-orbit counts, interval arrays, blocker histograms, and independent
positive keys are archived, so the aggregate is not standing in for the forty
individual decisions.

## Replay

```console
python3 scripts/certify_independent_h9_14_record18_a3_audit.py
python3 scripts/certify_independent_h9_14_record18_a3_audit.py --fresh
PYTHONPATH=.:src pytest -q tests/test_h9_14_record18_a3_independent_audit.py
```

The transcript is
`data/enumerations/h9_14_record18_a3_independent_audit.json`.
