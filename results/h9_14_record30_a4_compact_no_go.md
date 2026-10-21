# H9/14 Record 30: the compact `a=4` schedule is impossible

## Certified conclusion

All seven positive-A4 symmetry orbits in Record 30's compact `a=4`
Dehn--Sommerville schedule are incidence-impossible.  Consequently the whole
`a=4` schedule is excluded.

The scheduled data are

```text
f = (14, 88, 328, 802, 1340, 1540, 1184, 558, 124),
minimal A4 nonfaces = 4,
minimal A5 nonfaces = 32.
```

Six representatives require two additional A5s after their forced shadows;
each retains a transversal of order at most four.  The remaining
representative requires three additions.  Exactly 840 addition triples remove
all small transversals, but their numbers of order-five minimal transversals
are

```text
92 (90 families), 97 (120), 102 (480), 107 (150).
```

None has the scheduled 124 facets.  No ridge, topology, or Coxeter-label
candidate survives.

The production route literally traverses all 16,328,620 increasing triples.
An independent route uses a pivoted exact-cover recursion and canonical
memoization.  They agree orbit by orbit.  The frozen projection and
certificate digests are

```text
e11696f63b8cf3105292295e14370a598594bf96dd702eb1a415c31cfc2630ac
12918a797758ec9b14ed02657c8a98060c77b54abaa93d0e627dc31d0efb3bb5
```

Artifacts:

- `data/enumerations/h9_14_record30_a4_compact_no_go.json`
- `scripts/explore_h9_14_record30_a4_additions.cpp`
- `scripts/independent_h9_14_record30_a4_pivot_audit.cpp`
- `scripts/certify_h9_14_record30_a4_no_go.py`
- `tests/test_h9_14_record30_a4_no_go.py`

This is compact-only.  It gives no finite-volume noncompact conclusion and no
new hyperbolic Coxeter example.  Subsequent complete theorems also close
`a=2` and `a=1`, so only compact schedule `a=0` remains open.
