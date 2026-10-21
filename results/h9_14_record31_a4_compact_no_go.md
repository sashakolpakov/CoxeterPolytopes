# Record 31: the compact `a=4` schedule is impossible

## Theorem and scope

Every completion of H9/14 graph-shadow Record 31 in the compact
Dehn--Sommerville schedule with exactly four missing four-faces is
incidence-impossible.

This closes schedule `a=4` in addition to the previously excluded schedules
`a=5,6,7`.  The subsequent exact `a=3` theorem leaves only schedules
`a=0,1` open for other completions of Record 31; a later exact certificate
closes `a=2`.  This is not a whole-shadow
theorem and makes no statement about
finite-volume noncompact chambers, ideal faces, or parabolic diagrams.

## Upstream orbit theorem

The independently replayed positive-A4 certificate gives exactly twenty full
`Aut(P5 disjoint union 9I) = C2 x S9` orbits in schedule `a=4`.  The schedule
requires 33 missing five-faces and 114 facets.  Their forced A5 shadows split
as follows:

```text
shadow 32:  3 orbits, requiring 1 addition
shadow 31: 14 orbits, requiring 2 additions
shadow 30:  3 orbits, requiring 3 additions.
```

The upstream certificate is
`data/enumerations/h9_14_record31_positive_a5.json`, digest

```text
05a19d64d4784848ce5c59a1a7581e6bf59068828cd7fd92db6032f1a3176725
```

## One- and two-addition branches

For each of the seventeen low-addition representatives, the production
search constructs every independent A5 outside the forced shadow, applies
the exact residual vertex-link Macaulay upper bounds, and enumerates every
unordered addition family literally.

The candidate counts range from 6 to 449.  After the vertex caps, the exact
numbers of addition choices range from 6 to 27,693.  None destroys every
transversal of order at most four.  Consequently all seventeen
representatives fail before facet, ridge, topology, or Coxeter-label tests.

## Three-addition branches

The three remaining representatives are production-orbit indices `2,18,19`.
Their literal raw triple spaces and exact outcomes are:

| orbit | candidates | raw triples | cap-feasible triples | no-small-transversal triples | exact blocker-count histogram |
|---:|---:|---:|---:|---:|:---|
| 2  | 764 | 74,032,364 | 3,951,530 | 0 | empty |
| 18 | 462 | 16,328,620 | 1,794,605 | 270 | `105:150, 109:120` |
| 19 | 462 | 16,328,620 | 2,819,841 | 1,121 | `94:160, 95:120, 96:340, 99:80, 100:260, 102:150, 104:10, 115:1` |

The production search obtains the terminal families through an exact
minimum-remaining-values transversal-cover recursion.  Every completion must
destroy each current transversal of order at most four; branching on all
candidates that destroy a least-supported uncovered transversal is lossless.
Residual vertex-incidence caps are enforced at every node.

The 1,391 terminal families are then dualized literally on all `2^14`
subsets.  None has the scheduled blocker count 114.  Notice in particular
that the single 115-facet family in orbit 19 is not interchangeable with the
required 114-facet schedule.

## Independent replay

The independent verifier does not use the pivot recursion.  It directly
loops over every cap-feasible literal combination for all twenty
representatives: all singles, all pairs, and all

```text
3,951,530 + 1,794,605 + 2,819,841 = 8,565,976
```

cap-feasible triples.  It reconstructs every forced shadow, residual cap,
small-transversal set, and blocker count without importing project proof
modules.  It reproduces all candidate counts, all seventeen low-branch
zeros, the terminal counts `0,270,1121`, and both nonempty blocker
histograms exactly.

## Why downstream gates are empty

A scheduled `a=4` nerve must have a pure order-five blocker of cardinality
114 before it can have the scheduled face vector or enter ridge and topology
audits.  The exact target count here is zero.  Therefore:

```text
scheduled 114-facet targets:       0
ridge/topology targets:            0
unrestricted compact-label targets: 0
```

No Coxeter-label search is being silently replaced by a sparse or capped
search: there is simply no combinatorial nerve to label in this schedule.

## Current Record-31 compact frontier

```text
a=7,6: excluded by strict compact vertex-link positivity
a=5:   excluded by the independently replayed positive-A4/Macaulay gate
a=4:   excluded here by exact additions and blocker count
a=3:   excluded subsequently by exact completion/topology/label gates
a=2:   excluded by the later strict-link/incidence/compact-label theorem
a=1,0: open for other completions
```

The next branches are therefore `a=1,0`.

## Frozen artifacts

- `scripts/explore_h9_14_record31_a4_additions.py`
- `scripts/explore_h9_14_record31_a4_add3.cpp`
- `scripts/independent_h9_14_record31_a4_additions_audit.cpp`
- `scripts/certify_h9_14_record31_a4_no_go.py`
- `data/enumerations/h9_14_record31_a4_compact_no_go.json`
- `tests/test_h9_14_record31_a4_no_go.py`

The certificate digest is

```text
160d4b6dec78bea2c42246ac8b5e14b388eac8397e26448f9c7b05946115c3e0
```

Replay with

```bash
python3 scripts/certify_h9_14_record31_a4_no_go.py
python3 scripts/certify_h9_14_record31_a4_no_go.py --fresh
PYTHONPATH=.:src pytest -q tests/test_h9_14_record31_a4_no_go.py
```
