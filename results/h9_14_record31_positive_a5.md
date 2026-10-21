# Record 31: compact positive-`A4` census and `a=5` branch no-go

## Theorem and exact scope

For H9/14 graph-shadow record 31, the compact Dehn--Sommerville schedule
with exactly five missing four-faces is impossible.

This is a **compact-only partial theorem**.  It closes schedule `a=5`; later
exact certificates also close `a=4`, `a=3`, and now the complete `a=2`
schedule, so the current open schedules are `a=0,1`.  It neither proves a
whole-shadow compact no-go nor addresses
finite-volume noncompact chambers.  No parabolic or Vinberg endpoint audit is
performed, and no new hyperbolic Coxeter polytope is claimed.

## Graph, symmetry, and all schedules

The missing-edge graph is

```text
11--0--13--1--12    disjoint union    9 isolates,
```

with graph6 string `M?????????O?O?o??`.  Its full automorphism group is

```text
Aut(P5 disjoint union 9I) = C2 x S9,
```

of order `2 * 9! = 725,760`.

There are 87 independent pairs, 319 independent triples, 771 independent
four-sets, and 1,296 independent five-sets.  If `a` denotes the number of
minimal missing four-faces, Dehn--Sommerville gives exactly the following
compact schedules before strict vertex-link positivity:

| `a` | missing five-faces | facets | global `g4` | sum of vertex-link `g4` | strict compact sum |
|---:|---:|---:|---:|---:|:---:|
| 0 | 9  | 122 | 8 | 48 | yes |
| 1 | 15 | 120 | 7 | 42 | yes |
| 2 | 21 | 118 | 6 | 36 | yes |
| 3 | 27 | 116 | 5 | 30 | yes |
| 4 | 33 | 114 | 4 | 24 | yes |
| 5 | 39 | 112 | 3 | 18 | yes |
| 6 | 45 | 110 | 2 | 12 | no |
| 7 | 51 | 108 | 1 | 6  | no |

Every one of the fourteen compact vertex links has strict `g4 >= 1`, so its
sum is at least fourteen.  This excludes `a=6,7` immediately and leaves
exactly `a=0,...,5`.

## Exact local necessary conditions

The production and certifying code derive, rather than import, all 14 vertex,
87 meeting-edge, and 319 independent-triple link intervals.  Their 420 upper
endpoints have histogram

```text
0:30, 1:156, 2:84, 3:103, 4:2, 5:36, 6:9.
```

All lower endpoints are zero.  Applying this complete table to a single A4
reduces the 771 independent four-sets to 663 admissible rows.  The frozen
archive stores all 420 rows, not merely the histogram.

The forced A5 shadow of an A4 family is then subjected to strict local `g4`
and the degree-three Macaulay bounds in every vertex link.  If the forced
shadow has size `s`, the schedule requires exactly `b-s` additional A5s,
where `b=9+6a`; their total vertex incidence is `5(b-s)`, which must lie
between the sum of the fourteen exact residual lower and upper bounds.

## Positive branches `a=1,...,4`

The active-row/isolate-column quotient minimizes over all row permutations,
the reflection of `P5`, and the full `S9` action on isolates.  The exact
necessary-gate census is:

| `a` | raw A4 families | interval-feasible | strict/Macaulay viable | full graph orbits |
|---:|---:|---:|---:|---:|
| 1 | 663 | 663 | 495 | 6 |
| 2 | 219,453 | 210,177 | 83,946 | 49 |
| 3 | 48,352,811 | 42,603,687 | 1,415,316 | 143 |
| 4 | 7,978,213,815 | 6,218,797,221 | 49,476 | 20 |

These branches do **not** close at this gate.  Among the twenty `a=4`
orbits, three require one additional A5, fourteen require two, and three
require three.  This is the cheapest next all-completion branch.  For `a=3`,
the corresponding orbit counts by required additions `1,2,3,4` are
`55,40,37,11`; for `a=2`, additions `3,4,5,6` occur in `14,14,20,1` orbits;
for `a=1`, additions `5,6,7` occur in `1,1,4` orbits.

## Lossless exclusion of `a=5`

The `a=4` column traversal finds exactly 62,234 full-graph orbits of
interval-feasible four-families.  This quotient is lossless: its
orbit--stabilizer accounting recovers all `6,218,797,221` labelled families.

Every interval-feasible five-family has five interval-feasible four-row
subfamilies, because all local constraints are upper bounds.  Move any one of
those four-row subfamilies to its stored full-graph representative, then add
the moved fifth row.  Consequently, extending every representative by every
one of the 663 admissible rows covers every full-graph orbit of five-families.
The exact extension census is

```text
62,234 * (663 - 4) = 41,012,206 extension paths
                    34,615,149 interval-feasible paths
                       158,736 paths with forced shadow <= 39
                             0 strict/Macaulay survivors.
```

The last 158,736 paths split by forced-shadow size as

```text
37:24, 38:4,896, 39:153,816.
```

Every one fails the exact strict/aggregate local gate.  Hence no family of
five missing four-faces can occur in a compact completion of record 31.

## Structurally independent replay

The independent verifier does not extend the production four-row orbit list.
It directly enumerates all five active `P5` row-mask multisets and all
nine-column binary isolate multisets.  Sorting active rows and isolate
columns is a complete quotient by row order and `S9`; equal active rows may
produce harmless duplicate paths, but cannot omit a literal family.

Its independent census is

```text
6,188 active-row multisets
1,000 pass constraints supported entirely on P5
15,889,463 terminal column paths
12,797,343 interval-feasible paths
31,078 paths with forced shadow <= 39
0 strict/Macaulay survivors.
```

The path counts intentionally differ from the production quotient, while the
terminal zero agrees.

## Next exact completion gate

The most economical continuation is schedule `a=4`:

1. For each of its twenty A4 orbit representatives, enumerate exactly the
   required one, two, or three A5 additions outside the forced shadow.
2. Quotient additions only by the stabilizer of the fixed A4 family; this is
   lossless and avoids both duplicate literal families and an invalid use of
   the full graph group after a representative has been fixed.
3. Enforce every residual vertex-link lower/upper Macaulay bound during the
   traversal, followed by no-transversal-of-order-at-most-four, exact
   order-five blocker count `114`, ridge degree two, connected dual graph,
   complete face vector, and link Euler/mod-two homology audits.
4. Run the unrestricted compact label gate on every surviving nerve: all
   missing five-faces must be Lanner, every actual rank-nine vertex must be
   elliptic, and exact principal-block inertia/interlacing gates must leave no
   residual label state.

Only after `a=4` is closed should the same stabilizer-addition engine advance
through the one-addition `a=3` subbranch and then the remaining `a=3,2,1`
families.  Schedule `a=0` requires a separate lossless nine-A5 quotient.  The
existing selected `a=0` nerve is already compact-label impossible, but that
one selected no-go is not an all-completion theorem.

## Frozen artifacts

- `scripts/explore_h9_14_record31_positive_orbits.cpp`
- `scripts/independent_h9_14_record31_a5_audit.cpp`
- `scripts/certify_h9_14_record31_positive_a5.py`
- `data/enumerations/h9_14_record31_positive_a5.json`
- `tests/test_h9_14_record31_positive_a5.py`

The frozen certificate digest is

```text
05a19d64d4784848ce5c59a1a7581e6bf59068828cd7fd92db6032f1a3176725
```

Replay with

```bash
python3 scripts/certify_h9_14_record31_positive_a5.py
python3 scripts/certify_h9_14_record31_positive_a5.py --fresh
PYTHONPATH=.:src pytest -q tests/test_h9_14_record31_positive_a5.py
```
