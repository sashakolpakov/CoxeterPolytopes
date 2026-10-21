# Independent audit of the record-18 positive branches

## Result and scope

A standalone incidence-matrix quotient independently reproduces the complete
positive-
\(A_4\) scout for record 18, including every one of the seventy canonical
orbit keys and literal orbit sizes.  Independent downstream addition searches
then exclude all thirteen `a=4` orbits at a necessary incidence gate.

This is a **compact-only partial result**.  Schedules `a=0,1,2,3` remain open in
this certificate.  It is not a whole-shadow Coxeter theorem and says nothing
about finite-volume noncompact chambers.

## Independent reduction and quotient

The active missing-edge graph is the double star

```text
0--12--13--2
1--/     \--3
```

with automorphism group of order eight; vertices `4,...,11` are isolates.  The
independently derived vertex-, edge-, and triple-link table has 410 bounds.  In
particular, the two centre-vertex bounds have upper endpoint zero, so every
minimal missing four-face avoids `12,13`.  Thus

```text
735 independent A4s -> C(12,4) = 495 interval-admissible A4s.
```

Instead of scanning literal families, the replay represents an unordered
family by its binary row-incidence matrix.  Its four leaf columns are quotiented
by the order-eight active symmetry; the eight isolate columns are stored by
their complete row-pattern multiplicities, exactly quotienting by `S8`; every
row permutation is tested explicitly.  Orbit--stabilizer weights recover
`C(495,a)` before any prune in each branch.

| `a` | raw | interval-feasible | strict/aggregate | full orbits |
|---:|---:|---:|---:|---:|
| 1 | 495 | 495 | 145 | 3 |
| 2 | 122,265 | 120,969 | 9,144 | 14 |
| 3 | 20,092,215 | 19,449,199 | 188,944 | 40 |
| 4 | 2,471,342,445 | 2,310,730,275 | 14,812 | 13 |

Comparison with the production archive gives zero producer-only and zero
replay-only `(canonical key, literal size)` pairs in all four branches.  The
standalone quotient takes about 0.31 seconds after compilation on the audit
machine.

## Independent `a=4` addition gate

Twelve orbits require two additional five-faces.  A direct reverse-order pair
enumeration gives the exact classes

| orbit count | candidates | raw pairs | cap-feasible | no-small-transversal |
|---:|---:|---:|---:|---:|
| 3 | 252 | 31,626 | 3,501 | 0 |
| 2 | 246 | 30,135 | 3,475 | 0 |
| 5 | 126 | 7,875 | 1,375 | 0 |
| 2 | 56 | 1,540 | 540 | 0 |

The exceptional size-28 orbit needs four additions.  Its exact vertex
intervals consist of ten copies of `[0,2]`, two copies of `[1,2]`, and two
copies of `[0,0]`.  A separate coverage-descending, earliest-pivot transversal
search visits 845,919 nodes and independently obtains 3,960 addition families.
Complete blocker dualization gives

```text
82:160, 83:1920, 85:240, 87:120,
88:900, 92:260, 94:360, 102:0.
```

Hence no exceptional family has the scheduled pure 102-facet blocker, and the
entire `a=4` schedule is incidence-impossible.  This replay has a different
candidate order and search tree from production.

## Necessity and caveats

The local intervals are consequences of strict compact spherical-link
positivity and the Macaulay inequalities.  The upper shadow consists of every
missing five-set already forced by a missing four-set.  A base transversal of
order at most four must be made disjoint from at least one newly added
five-face; otherwise the completed blocker has a facet complement smaller than
five and cannot be the pure boundary of a compact nine-polytope.  Finally, the
complete blocker must consist of exactly the scheduled number of five-sets.

These are necessary incidence conditions, not Vinberg's finite-volume test.
Ideal links and parabolic endpoint data are absent.

## Replay

```console
python3 scripts/certify_independent_h9_14_record18_positive_audit.py
python3 scripts/certify_independent_h9_14_record18_positive_audit.py --fresh
PYTHONPATH=.:src pytest -q \
  tests/test_h9_14_record18_positive_independent_audit.py
```

The frozen audit is
`data/enumerations/h9_14_record18_positive_independent_audit.json`.
