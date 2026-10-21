# Record 32: whole-shadow compact-Coxeter no-go

## Theorem and exact scope

Under the standard compact `H9/14` missing-face reductions, **no incidence
completion of record 32 is the nerve of a compact hyperbolic Coxeter
9-polytope with fourteen facets**.  This excludes the whole record-32
missing-edge graph shadow, not merely its previously selected convex nerve
and selected `a=0` abstract sphere candidate.

This is a compact theorem only.  It uses strict ellipticity of every actual
rank-nine vertex diagram and compact Lanner missing faces.  It does not
enumerate parabolic supports, ideal vertices, cusp extensions, or Vinberg's
finite-volume criterion.  It therefore neither constructs nor excludes a
finite-volume noncompact chamber.

The abstract incidence targets below are used only as necessary objects.
Except for the separately established selected `a=2` nerve, this theorem
does not assert ordinary polytopality, PL-sphere status, simple connectivity,
or integral homology.  Ridge-invalid targets require no stronger topology
inference.

Record 32 has graph6 string

```text
M?????????O?O?g??
```

and missing edges

```text
(0,11), (1,12), (0,13), (2,13).
```

The active graph is `P4 disjoint union K2` on
`(0,1,2,11,12,13)`; vertices `3,...,10` are eight isolates.  The full graph
action is `Aut(P4 disjoint union K2) x S8`, with active factor of order four.

## Complete Dehn--Sommerville schedules

Let `a` be the number of missing A4s and `b` the total number of missing
independent five-subsets, including the upper shadow forced by the A4s.  The
six exact schedules are

| `a` | `b` | rank-nine facets | sum of fourteen vertex-link `g4` values |
|---:|---:|---:|---:|
| 0 | 6 | 114 | 36 |
| 1 | 12 | 112 | 30 |
| 2 | 18 | 110 | 24 |
| 3 | 24 | 108 | 18 |
| 4 | 30 | 106 | 12 |
| 5 | 36 | 104 | 6 |

Every compact vertex link has strict `g4 >= 1`, so the total is at least
fourteen.  This excludes `a=4,5` and leaves `a=0,1,2,3` for exact census.  No
claim about ideal links is inferred.

## Positive missing-four branches

There are 762 independent A4s.  Exact vertex-, edge-, and triple-link
intervals leave 602.  Every lower interval endpoint is zero, so an A4
rejected individually cannot be repaired by adding another A4; the
602-element ground set is lossless.

| `a` | raw families | interval-feasible | strict/Macaulay viable | full graph-action orbits |
|---:|---:|---:|---:|---:|
| 1 | 602 | 602 | 322 | 4 |
| 2 | 180,901 | 172,039 | 2,156 | 6 |
| 3 | 36,180,200 | 30,954,448 | 56 | 1 |

A no-project C++ replay independently derives all 419 local constraints from
the four missing edges, including 130 cap-two constraints, and reproduces
the complete gate histograms and orbit-size histograms in all three rows.

### The `a=3` branch

The 56 viable triples form one orbit represented by

```text
(1,2,3,4), (2,3,4,11), (3,4,11,12).
```

Its forced A5 shadow has size 22 and it needs two additional A5s.  There are
56 candidates and 1,540 pairs; 540 satisfy all exact vertex caps, and every
one has a transversal of order at most four.  Thus `a=3` is empty before
labels.

### The `a=2` branch

The six strict-shadow orbit representatives have exact
`(shadow, additions, candidates, raw, cap-feasible, no-small, target)` rows

```text
(18,0,  0,         1,         1,   0,0)
(17,1,  1,         1,         1,   1,1)
(15,3,462,16,328,620, 1,794,605, 840,0)
(16,2, 56,     1,540,       540,   0,0)
(16,2, 56,     1,540,       945,   0,0)
(16,2, 56,     1,540,     1,540,   0,0).
```

The sole target is exactly the older selected record-32 nerve

```text
A4 = (2,3,4,5), (3,4,5,11)
A5 = (6,7,8,9,10).
```

Its previously frozen unrestricted Lanner/actual-vertex obstruction closes
the entire `a=2` branch.  The independent replay reproduces all six addition
rows and literal equality with that selected nerve.

### The `a=1` branch

For a fixed A4, the isolates split into those inside and outside it.  The
lossless census admits every independent active subset in every additional
A5 row and represents both isolate classes by their complete row-incidence
column multisets.  This quotients only the two genuine isolate symmetric
groups in the fixed-A4 stabilizer; no singleton-active-row ansatz is used.

| A4 representative | column paths | distinct-row paths | no-small paths | literal targets | target orbits | ridge-valid orbits |
|:---|---:|---:|---:|---:|---:|---:|
| `(3,4,5,6)` | 1,161 | 1,102 | 9 | 1 | 1 | 1 |
| `(2,3,4,5)` | 43,229 | 42,865 | 271 | 3 | 3 | 1 |
| `(1,2,3,4)` | 2,323,448 | 2,278,508 | 12,873 | 0 | 0 | 0 |
| `(2,3,4,11)` | 6,068,459 | 5,934,396 | 19,836 | 8 | 6 | 0 |

The four A4s lie in different full graph-action orbits.  Hence there are
exactly ten scheduled `a=1` target orbits, and exactly two have blocker
histogram `{5:112}` and ridge histogram `{2:504}`.  The independent
implementation reproduces every path count, complete facet histogram,
target orbit, blocker, and ridge decision.

Complete Lanner joins leave 10,368 and 864 sparse locally elliptic states on
the two ridge targets.  Positive sparse counts are not hyperbolic candidates
and are not obstructions by themselves.  In each target, choose the displayed
disjoint order-four and order-five Lanner supports.  If all twenty cross
pairs are right-angled, their principal Gram is an orthogonal sum with
negative index two, contradicting ambient compact-H9 inertia `(9,1,4)`.
Thus at least one cross edge is non-right.  Lower it to `m=3` and every other
unforced finite edge to `m=2`; Perron--Frobenius monotonicity makes this a
necessary actual-vertex test.

| target | sparse states | cross pairs | single-`m=3` tests | passing |
|:---|---:|---:|---:|---:|
| first | 10,368 | 20 | 207,360 | 0 |
| second | 864 | 20 | 17,280 | 0 |

The first-bad restriction histograms are `{7:206208,8:1152}` and
`{7:17280}`.  A no-project implementation independently reproduces the
Lanner states and every cross test.  The twenty-pair argument uses no
symmetry assumption: edges to facets outside the chosen principal pair
cannot alter that principal submatrix.  Hence the zero tests exclude
arbitrary stronger and additional finite labels.

## Lossless `a=0` census

Here the six missing A5s are minimal.  All 24 independent active masks are
allowed, with size histogram

```text
0:1, 1:6, 2:11, 3:6.
```

The six rows are unordered.  Each isolate is represented by its exact six-bit
row-incidence column, and the eight isolate columns form a complete multiset.
This is a lossless `S8` quotient: every labelled family can first have its
rows sorted and then its isolates renamed to put the columns in nondecreasing
order.  Equal active rows may overcount a final literal family, but cannot
omit one.

The six minimal covers of the active missing-edge graph all have order three.
Every forbidden transversal of order at most four is one such cover, possibly
extended by one active vertex or one isolate.  The active-row and allowed-
column tests exhaust exactly those alternatives.  A full literal postcheck
on every nonrepeated path records zero discrepancies.

The production stream gives

```text
22,513 active-row multisets after caps
14,731 after the active-extension test
1,976 distinct column-distribution problems
179,785,655 exact column-multiset paths
3,513,588 paths with repeated A5 rows
176,272,067 blocker-valid paths
1,263 paths with the scheduled 114 facets
1,025 production literal representatives.
```

A separately written no-project C++ implementation uses numeric active-mask
order, explicit column-multiplicity recursion, `std::bitset` facet words, and
the literal small-transversal postcheck.  It reproduces every path count and
emits 969 literal representatives under its deliberately different choice of
isolate representative.  A third no-project comparator exhausts the genuine
`Aut(P4 disjoint union K2) x S8 x S6` action.  Both lists reduce to exactly
292 orbits, with zero set differences and common orbit-key SHA-256

```text
97f630bd464ec01d77087019c5dd02b549f0cbaebd466d5e621ae065797dc848.
```

The unequal literal representative counts are normalization-dependent; the
identical exact orbit set is the invariant completeness statement.

The production target archive SHA-256 is

```text
ae84f9c8dde0f66cd4d8fc8c0b9dee2afab985c7d484f2e052709ad3c3734938.
```

Independent literal blocker dualization leaves 941 families with histogram
`{5:114}`.  The other 84 have additional minimal order-six transversals.
Only indices

```text
6, 27, 43, 57, 132, 139, 436, 439, 481, 574, 626, 629
```

have ridge histogram `{2:513}`.  The other pure targets fail the ridge law;
no PL or full-homology conclusion is used.

For all twelve ridge candidates, the complete 420-placement A5-Lanner join
against all 114 actual rank-nine vertices has an uncapped zero terminal.
The independent implementation deliberately uses a different support order
and reaches the same twelve zero terminals.  Perron--Frobenius monotonicity
then excludes arbitrary additional non-right finite labels.

## Independent replay and conclusion

The hostile replay imports no record-32 producer.  It independently covers

- all six schedules and the strict compact cut;
- the complete 602-A4 positive family census and all orbit histograms;
- the distinct `a=0` path/target generator and exact 292-orbit comparison;
- literal blocker and ridge audits for every `a=0` and `a=1` target;
- all order-four and order-five Lanner placements and actual vertices;
- both unrestricted `a=1` cross gates; and
- the `a=2` selected-nerve chain and empty `a=3` addition branch.

Combining these independently replayed branches proves the whole record-32
graph shadow impossible on the compact Coxeter track.  No finite-volume
noncompact conclusion follows.

## Frozen proof objects and reproduction

The promoted compact-only certificate is
[`h9_14_record32_all_completions.json`](../data/enumerations/h9_14_record32_all_completions.json),
verified by
[`certify_h9_14_record32_all_completions.py`](../scripts/certify_h9_14_record32_all_completions.py).
Its file SHA-256 is
`600d107e9739085897349cf905d411df1af260afb75facb589c0edadf74a52b0`
and its internal canonical certificate SHA-256 is
`486b063921841de9aa3c4611ea171801c3368c1c1c7c5cdcef938fabedfb69d9`.
The chained hostile archive SHA-256 is
`381c97ac36f85d146f96601e19de0200efe68cd48bb7ad98728035f4cc0a0417`.
Its historical production predecessor deliberately retains its pre-replay
`theorem=false` flags; the promoted object chains to it without rewriting
that history.

Run the frozen aggregate verifier with

```console
PYTHONPATH=src:. python3 scripts/certify_h9_14_record32_all_completions.py
```

The hostile archive/source check and fresh positive-orbit replay are

```console
python3 scripts/certify_h9_14_record32_whole_shadow_hostile.py
python3 scripts/certify_h9_14_record32_whole_shadow_hostile.py --positive
python3 scripts/certify_h9_14_record32_whole_shadow_hostile.py --downstream
```

The `--downstream` mode freshly runs the positive-orbit, additions, and
blocker/ridge/Lanner/cross programs and compares their complete canonical
outputs. Passing `--full` additionally rebuilds the 179,785,655-path `a=0`
stream and compares both target lists under the complete graph/row action.  The focused
regression suite is

```console
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_h9_14_record32_whole_shadow_scout.py \
  tests/test_h9_14_record32_all_completions.py
```

The production and replay sources are listed in
[`h9_14_record32_whole_shadow_candidate.md`](h9_14_record32_whole_shadow_candidate.md)
and the independent transcript is in
[`h9_14_record32_whole_shadow_hostile_audit.md`](h9_14_record32_whole_shadow_hostile_audit.md).
