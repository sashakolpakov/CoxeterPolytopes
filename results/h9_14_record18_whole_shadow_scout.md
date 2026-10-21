# Record 18 whole-shadow scout: all positive-A4 schedules

## Result and exact scope

This production scout gives a lossless compact-only census through the first
completion gate for every record-18 branch with at least one missing
four-face.  It leaves exactly

```text
a=1:  3 forced-shadow orbits
a=2: 14 forced-shadow orbits
a=3: 40 forced-shadow orbits
a=4: 13 forced-shadow orbits
```

under the complete missing-edge graph action.  A subsequent exact addition
gate already closes all thirteen `a=4` orbits before Coxeter labels.

This is not a whole-shadow theorem.  The `a=1,2,3` additions and the `a=0`
branch remain to be completed, and the production census has not yet been
reproduced by an independent implementation.  The abstract objects here are
only necessary incidence objects: no PL-sphere or ordinary-polytopality claim
is made.  Nothing here constructs or excludes a finite-volume noncompact
chamber; ideal vertices, parabolic diagrams, and Vinberg's endpoint criterion
are not modeled.

Record 18 has graph6 word

```text
M???????????o?K?_
```

and missing edges

```text
(0,12), (1,12), (2,13), (3,13), (12,13).
```

The active graph is the double star with centres `12,13`, two leaves on each
side, and automorphism group of order eight.  Vertices `4,...,11` are eight
isolates, so the full graph action is the active order-eight group times
`S8`.

## Necessary prune 1: all Dehn--Sommerville schedules

Let `a` be the number of minimal missing A4s and let `b` be the number of all
missing independent A5s, including the upper shadow forced by the A4s.  The
complete schedules are

| `a` | `b` | rank-nine facets | sum of the fourteen vertex-link `g4` values |
|---:|---:|---:|---:|
| 0 | 8 | 110 | 42 |
| 1 | 14 | 108 | 36 |
| 2 | 20 | 106 | 30 |
| 3 | 26 | 104 | 24 |
| 4 | 32 | 102 | 18 |
| 5 | 38 | 100 | 12 |
| 6 | 44 | 98 | 6 |

For a compact Coxeter nerve every vertex link is a simplicial sphere and has
strict `g4 >= 1`.  Its summed value is therefore at least fourteen.  This
necessary compact condition excludes `a=5,6` and no lower branch.  The step
does not apply to ideal links in a noncompact chamber.

## Necessary prune 2: lossless local A4 intervals

The five missing edges admit 735 independent four-subsets.  The exact
vertex-, edge-, and triple-link inequalities give 410 incidence constraints,
with upper-cap histogram

```text
0:112, 1:36, 2:132, 3:90, 4:32, 5:8.
```

Every lower endpoint is zero.  Consequently, an A4 violating an individual
cap-zero condition can never be repaired by adding other A4s.  Removing those
240 A4s is therefore lossless, leaving an exact 495-element ground set.

For families of size two, three, and four the program imposes respectively
all cap-one, cap-two, and cap-three constraints.  Higher caps cannot be
violated by four rows.  Thus the pair graph, triple signatures, and quadruple
signatures are an exact implementation of every local interval; they are not
a heuristic compatibility filter.

## Necessary prune 3: the forced-shadow gate

For an A4 family `Q`, let `S(Q)` be its full independent-A5 upper shadow.  At
vertex `v`, write `q_v` and `s_v` for the incidences of `Q` and `S(Q)`.  The
two local quantities before adding the remaining minimal A5s are

```text
previous_v = base_g3(v) - q_v,
before_v   = base_g4(v) + 6 q_v - s_v.
```

Compactness requires both quantities to be positive.  If `c_v` additional
A5s contain `v`, strict local `g4` and the degree-three Macaulay inequality
give the necessary interval

```text
max(0, before_v - previous_v^{<3>}) <= c_v <= before_v - 1.
```

Since every added A5 has five vertices, summing these lower and upper bounds
must bracket five times the scheduled number of additions.  This is the
`strict/aggregate` gate below.  Passing it does not prove that additions
exist; it only reduces each branch to finitely many exact completion targets.

## Exact positive-A4 census

| `a` | raw A4 families | interval-feasible | strict/aggregate viable | full graph orbits |
|---:|---:|---:|---:|---:|
| 1 | 495 | 495 | 145 | 3 |
| 2 | 122,265 | 120,969 | 9,144 | 14 |
| 3 | 20,092,215 | 19,449,199 | 188,944 | 40 |
| 4 | 2,471,342,445 | 2,310,730,275 | 14,812 | 13 |

The exact viable shadow histograms are

```text
a=1: shadow 8 / additions 6 :       145
a=2: shadow 15 / additions 5:       128
     shadow 16 / additions 4:     9,016
a=3: shadow 22 / additions 4:       112
     shadow 23 / additions 3:    14,000
     shadow 24 / additions 2:   174,832
a=4: shadow 28 / additions 4:        28
     shadow 30 / additions 2:    14,784.
```

All seventy orbit keys, exact orbit sizes, and literal representatives are in
`data/enumerations/h9_14_record18_positive_A4_scout.json`.

## Why the `a=4` quotient is lossless

Directly scanning 2.47 billion labelled quadruples would obscure the proof.
Instead every A4 row is split into its subset of the six active vertices and
its isolate incidence.  The four rows are unordered and the eight isolates
are encoded by their complete four-bit row-incidence columns.

The program enumerates all 17,550 multisets of four independent active rows
and every exact multiset of eight isolate columns having the required row
sums.  The transcript is

```text
exact column-multiset paths                    2,228,923
duplicates from permuting equal active rows     448,772
paths with repeated full A4 rows                  90,666
distinct S8 family orbits                      1,689,485
interval-feasible S8 orbits                      280,869
strict/aggregate viable S8 orbits                     54
```

For each column multiset, its labelled `S8` weight is computed by
orbit--stabilizer.  The stabilizer is replayed exactly: for every row
permutation preserving the active rows, the transformed column multiset is
compared literally, and each matching column class contributes the factorial
of its multiplicity.  These weights sum to 2,310,730,275 on the interval
orbits and to 14,812 on the viable orbits.  Canonicalization under all eight
active graph automorphisms then gives thirteen full orbits; their archived
sizes sum again to 14,812 and each divides `8*8!`.

This proves that the quotient loses neither nonsingleton active rows nor
families with repeated active parts.

## The complete `a=4` addition gate

Twelve of the thirteen orbits have forced shadow 30 and require two added
A5s.  Literal enumeration gives four candidate/cap-feasible classes:

| number of orbit representatives | candidates | raw pairs | cap-feasible pairs | no-small-transversal pairs |
|---:|---:|---:|---:|---:|
| 3 | 252 | 31,626 | 3,501 | 0 |
| 2 | 246 | 30,135 | 3,475 | 0 |
| 5 | 126 | 7,875 | 1,375 | 0 |
| 2 | 56 | 1,540 | 540 | 0 |

The exceptional size-28 orbit is

```text
(0,2,10,11), (1,2,10,11), (0,3,10,11), (1,3,10,11).
```

Its forced shadow has size 28, so four A5s must be selected from 764
candidates.  The derived addition intervals are `[0,2]` at vertices
`0,...,9`, `[1,2]` at vertices `10,11`, and `[0,0]` at the centres `12,13`.
Thus restricting candidates to `W`, rejecting a third incidence with the
exact `twice` mask, and requiring `10,11` at the leaves of the search are
consequences of the local inequalities, not an ansatz.  The raw space has
14,084,657,251 unordered quadruples.  The base
clutter has exactly 47 transversals of order at most four.  A cap-aware exact
cover search visits 1,124,699 memoized states and finds 3,960 addition
families that destroy all 47 small transversals.  Direct blocker dualization
gives facet-count histogram

```text
82:160, 83:1920, 85:240, 87:120,
88:900, 92:260, 94:360.
```

No family has the scheduled 102 facets.  Therefore the entire `a=4` branch
is empty at a necessary incidence gate, before topology or Coxeter labels.

## Remaining gates

1. Enumerate the additional A5s on the 40 `a=3` representatives, beginning
   with the literal two-addition cases, and retain exact orbit keys through
   every blocker/ridge decision.
2. Repeat for the fourteen `a=2` and three `a=1` representatives using the
   fixed-A4 stabilizer column quotient where raw combinations are too large.
3. Enumerate `a=0` losslessly from all eight A5 rows.  With
   `L={0,1}`, `R={2,3}`, and `P={12,13}`, the three minimal graph covers are
   `P`, `L+{13}`, and `R+{12}`.  The corresponding `Y/Z/X` blocker sections
   must all be generated; no balanced-cut or singleton-active-row ansatz is
   permitted.  Quotient by the full active order-eight group and `S8`.
4. Only topology-valid targets proceed to complete Lanner placement, every
   actual rank-nine vertex test, and exact Gram-index obstruction.
5. Require a separately written replay before promoting any compact theorem.

## Artifacts and reproduction

```text
scripts/explore_h9_14_record18_positive_orbits.cpp
scripts/certify_h9_14_record18_positive_A4_scout.py
data/enumerations/h9_14_record18_positive_A4_scout.json
scripts/explore_h9_14_record18_a4_exception.cpp
```

The warning-clean production replay took 21.1 seconds for all four positive
branches on the development machine; the exceptional addition gate took 5.6
seconds.  Run

```console
PYTHONPATH=src:. python3 scripts/certify_h9_14_record18_positive_A4_scout.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record18_positive_A4_scout.py --fresh
clang++ -std=c++20 -O3 -Wall -Wextra -pedantic \
  scripts/explore_h9_14_record18_a4_exception.cpp -o /tmp/record18-a4-exception
/tmp/record18-a4-exception
```
