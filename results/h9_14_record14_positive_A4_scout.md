# Record 14: compact positive-`A4` scout

Status: the Dehn--Sommerville schedules, every vertex/edge/triple link
interval, and every positive-`A4` strict/aggregate orbit have been derived
directly from the record-14 missing-edge graph.  The summed strict-link test
immediately excludes schedules `a=5,6`.  None of the positive surviving
schedules `a=1,2,3,4` closes at this scout gate.

This is a compact-only necessary-incidence scout.  It proves neither that a
survivor is a sphere nor that it has a Coxeter labeling or hyperbolic Gram
matrix.  It does not test finite-volume noncompact chambers.

## Graph-derived arithmetic

For

```text
M???????????o?w??
```

the missing edges are

```text
0--12, 1--12, 0--13, 1--13, 2--13.
```

There are 86 independent pairs, 310 independent triples, 735 independent
four-sets, and 1,212 independent five-sets.  Dehn--Sommerville symmetry and
the strict generalized-lower-bound/Lanner condition give exactly seven
schedules:

| `a` | missing independent five-sets | facets | global `g4` | summed vertex-link `g4` | strict sum |
|---:|---:|---:|---:|---:|:---:|
| 0 | 8  | 110 | 7 | 42 | yes |
| 1 | 14 | 108 | 6 | 36 | yes |
| 2 | 20 | 106 | 5 | 30 | yes |
| 3 | 26 | 104 | 4 | 24 | yes |
| 4 | 32 | 102 | 3 | 18 | yes |
| 5 | 38 | 100 | 2 | 12 | no |
| 6 | 44 | 98  | 1 | 6  | no |

Strict `g4>=1` in each of fourteen vertex links requires the summed value to
be at least 14, so `a=5,6` are impossible before any search.

## Exact local table

The source independently constructs 14 vertex-link, 86 meeting-edge-link,
and 310 independent-triple-link constraints.  The 410 upper endpoints have
histogram

```text
0:66, 1:146, 2:57, 3:95, 4:37, 5:9.
```

In particular, the vertex-13 interval is `[0,0]`, the vertex-12 interval is
`[0,1]`, the vertex-0 and vertex-1 intervals are `[0,3]`, and the vertex-2
interval is `[0,4]`.  Applying the complete table to a single missing
four-face reduces the 735 independent four-sets to exactly 579 candidates.
The frozen archive records every one of the 410 rows, not only this histogram.

## Positive-`A4` census

The order-two automorphism group on the five active vertices and the full
permutation group on the nine isolates are both quotiented.  A canonical key
stores the active row masks and every complete binary isolate-column pattern;
all row permutations and both active automorphisms are minimized explicitly.

| `a` | raw families | interval-feasible | strict/aggregate | full orbits |
|---:|---:|---:|---:|---:|
| 1 | 579 | 579 | 369 | 5 |
| 2 | 167,331 | 160,245 | 60,810 | 36 |
| 3 | 32,183,329 | 28,293,741 | 580,272 | 59 |
| 4 | 4,634,399,376 | 3,571,339,779 | 28,854 | 13 |

For `a=4`, the collision-free incidence-matrix traversal has 1,365 active-row
multisets and 469,543 column-multiplicity paths.  Removing 159,572 duplicates
from equal active rows and 19,829 repeated-row paths leaves 290,142 distinct
matrix orbits.  Exactly 102,295 are orbits on the 579-element `A4` ground set,
and their orbit--stabilizer weights recover all 4,634,399,376 raw families.
After the local table, 66,860 isolate-column orbits remain; 18 pass the
strict/aggregate shadow gate and merge to 13 full graph orbits of total raw
weight 28,854.

Thus the first immediate closure is only the schedule cut `a=5,6`.  The
positive schedules need additional-five-face, blocker-purity, topology, and
eventually unrestricted Coxeter-label audits.

## Artifacts

- `scripts/explore_h9_14_record14_positive_orbits.cpp`
- `scripts/certify_h9_14_record14_positive_A4_scout.py`
- `data/enumerations/h9_14_record14_positive_A4_scout.json`
- `tests/test_h9_14_record14_positive_A4_scout.py`
