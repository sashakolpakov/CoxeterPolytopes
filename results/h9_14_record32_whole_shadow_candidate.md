# Record 32 whole-shadow compact no-go: production candidate

## Status and scope

The production census now excludes every record-32 completion on the compact
Coxeter track.  **This is not yet promoted to a theorem:** the new lossless
completion census and the new all-target label passes still require a
structurally independent replay.  Until that replay agrees, the result is a
candidate proof package and the older selected-nerve theorems remain the
formally established claims.

The candidate conclusion is compact-only.  It uses strict ellipticity of
rank-nine actual-vertex diagrams and compact Lanner missing faces.  It does
not enumerate parabolic supports, ideal vertices, cusp extensions, or the
Vinberg finite-volume criterion.  Consequently it neither constructs nor
excludes a finite-volume noncompact chamber.

Record 32 has graph6 string

```text
M?????????O?O?g??
```

and missing edges

```text
(0,11), (1,12), (0,13), (2,13).
```

The active missing-edge graph is `P4 disjoint union K2` on
`(0,1,2,11,12,13)`; vertices `3,...,10` are eight isolates.  Its active
automorphism group has order four, and the full graph action is
`Aut(P4 disjoint union K2) x S8`.

## Exact compact schedules

Let `a` be the number of missing A4s and let `b` be the total number of
graph-independent missing five-subsets, including every A5 in an A4 upper
shadow.  Direct Dehn--Sommerville replay and a double count of vertex-link
faces give

| `a` | `b` | rank-nine facets | sum of vertex-link `g4` |
|---:|---:|---:|---:|
| 0 | 6 | 114 | 36 |
| 1 | 12 | 112 | 30 |
| 2 | 18 | 110 | 24 |
| 3 | 24 | 108 | 18 |
| 4 | 30 | 106 | 12 |
| 5 | 36 | 104 | 6 |

Every compact vertex link has strict `g4 >= 1`, so the total is at least
fourteen.  This excludes `a=4,5`; it makes no statement about ideal links.

## Positive missing-four reduction

There are 762 independent A4s and 602 pass the exact singleton vertex-,
edge-, and triple-link intervals.  Every lower interval endpoint is zero.
Thus an A4 rejected individually cannot be repaired by adding other A4s, and
the 602-element ground set is lossless.

| `a` | raw A4 families | interval-feasible | strict/Macaulay viable | graph-action orbits |
|---:|---:|---:|---:|---:|
| 1 | 602 | 602 | 322 | 4 |
| 2 | 180,901 | 172,039 | 2,156 | 6 |
| 3 | 36,180,200 | 30,954,448 | 56 | 1 |

For `a=3`, the exact interval generator encodes cap-zero constraints on
singletons, cap-one constraints on pairs, and cap-two constraints on triples.
The 56 viable families form one orbit, represented by

```text
(1,2,3,4), (2,3,4,11), (3,4,11,12).
```

Its forced A5 shadow has size 22.  Of the 1,540 pairs of its 56 possible
additional A5s, 540 satisfy the exact vertex caps and all 540 have a
transversal of order at most four.  Hence `a=3` is empty before labels.

### Complete `a=2` additions

The six viable A4-pair orbits have the following exact addition censuses.

| orbit | forced shadow | added A5s | candidates | raw choices | cap-feasible | no small transversal | 110-facet targets |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 18 | 0 | 0 | 1 | 1 | 0 | 0 |
| 1 | 17 | 1 | 1 | 1 | 1 | 1 | 1 |
| 2 | 15 | 3 | 462 | 16,328,620 | 1,794,605 | 840 | 0 |
| 3 | 16 | 2 | 56 | 1,540 | 540 | 0 | 0 |
| 4 | 16 | 2 | 56 | 1,540 | 945 | 0 | 0 |
| 5 | 16 | 2 | 56 | 1,540 | 1,540 | 0 | 0 |

The sole target is

```text
A4 = (2,3,4,5), (3,4,5,11)
A5 = (6,7,8,9,10).
```

It is literally the previously certified ordinary convex record-32 nerve.
Its already frozen unrestricted compact-label obstruction therefore closes
the `a=2` target, subject only to independent replay of the new six-orbit
incidence reduction.

## Lossless `a=1` census

For a fixed A4, the eight graph isolates split into those inside and outside
that A4.  The production enumeration allows every independent active subset
in every additional-A5 row and represents each isolate class by the complete
multiset of its row-incidence columns.  This is a lossless quotient by the
two genuine isolate symmetric groups in the fixed-A4 stabilizer.  It is not a
singleton-active-row normalization.

| orbit representative | shadow / additions / candidates | raw additions | column paths | distinct-row paths | no-small paths | canonical 112-facet families | stabilizer orbits | ridge-valid orbits |
|:---|:---:|---:|---:|---:|---:|---:|---:|---:|
| `(3,4,5,6)` | 10 / 2 / 664 | 220,116 | 1,161 | 1,102 | 9 | 1 | 1 | 1 |
| `(2,3,4,5)` | 9 / 3 / 672 | 50,351,840 | 43,229 | 42,865 | 271 | 3 | 3 | 1 |
| `(1,2,3,4)` | 8 / 4 / 672 | 8,421,345,240 | 2,323,448 | 2,278,508 | 12,873 | 0 | 0 | 0 |
| `(2,3,4,11)` | 8 / 4 / 1,252 | 101,887,955,625 | 6,068,459 | 5,934,396 | 19,836 | 8 | 6 | 0 |

Because the four distinguished A4s lie in different full graph-action
orbits, these rows give exactly ten scheduled `a=1` incidence orbits, of
which exactly two pass the ridge law.

Every scheduled target has 112 order-five minimal transversals.  The
ridge-invalid targets are rejected only as pseudomanifolds; no PL-sphere,
ordinary-polytopality, or full-homology conclusion is inferred.  Exactly two
fixed-A4 target orbits have ridge histogram `{2:504}`:

```text
Q=(3,4,5,6)
R=(2,7,8,9,10)
B=(7,8,9,10,11)
```

and

```text
Q=(2,3,4,5)
A=(6,7,8,9,10)
R=(3,4,5,6,11)
B=(7,8,9,10,13).
```

Complete Lanner joins and every actual rank-nine restriction leave 10,368
and 864 sparse locally elliptic states, respectively.  These positive sparse
counts are not hyperbolic candidates and are not themselves obstructions.

For each target, the displayed `Q` and the second displayed disjoint support
(`R` in the first target, `A` in the second) are disjoint Lanners.  If every
cross pair between them is right-angled, their principal Gram is an
orthogonal sum with negative index two, contradicting ambient compact-H9
inertia `(9,1,4)`.  Thus at least one of the twenty cross pairs must be
non-right.  Retaining the forced Lanner state, lowering that edge to label
three, and lowering every other unforced finite edge to label two is
necessary by Perron--Frobenius monotonicity.  The complete tests are

| target | sparse states | cross pairs | single-`m=3` tests | passing |
|:---|---:|---:|---:|---:|
| first | 10,368 | 20 | 207,360 | 0 |
| second | 864 | 20 | 17,280 | 0 |

For the first target the first failed actual-vertex restriction histogram is
`{7:206208, 8:1152}`; for the second it is `{7:17280}`.  The twenty-pair
scope does not use an unproved symmetry: arbitrary edges to facets outside
the chosen disjoint principal pair cannot alter that principal submatrix.
Hence either its cross block is zero and already has index two, or it contains
a tested non-right edge.  The zero single-edge counts therefore close
arbitrary stronger and additional finite labels.

## Lossless `a=0` census

Here all six missing A5s are minimal.  The active graph admits all 24
independent masks, with size histogram

```text
0:1, 1:6, 2:11, 3:6.
```

All are available in all six unordered rows.  The strict A5 incidence caps
are

```text
(1,2,3,5,5,5,5,5,5,5,5,3,2,1).
```

The six inclusion-minimal covers of the active missing-edge graph all have
order three.  Any forbidden transversal of order at most four is therefore
one of these covers, possibly extended by one active vertex or one isolate.
The active-row filter handles the first possibility and the allowed-column
filter handles the second, so this prefilter is logically equivalent to the
literal small-transversal postcheck for the encoded paths.  A complete second
stream with that literal postcheck enabled records zero failures among all
nonrepeated paths.

The exact row/column stream gives

```text
22,513 active-row multisets after caps
14,731 after the exact active-cover filter
5,389 active configurations with a positive column distribution
1,976 distinct column-distribution problems after memoization
179,785,655 exact S8 column-multiset paths
3,513,588 paths with repeated A5 rows
176,272,067 blocker-valid paths
1,263 paths with the scheduled 114 facets
1,025 exact literal unordered A5 families.
```

The exact target archive has SHA-256
`ae84f9c8dde0f66cd4d8fc8c0b9dee2afab985c7d484f2e052709ad3c3734938`.

A separately written no-project C++ stream uses numeric active-mask order and
explicit column-multiplicity recursion.  It independently reproduces all path
counts and emits 969 literal representatives under its deliberately different
choice of isolate representative.  A no-project comparator then exhausts the
true `Aut(P4 disjoint union K2) x S8 x S6` action: the production 1,025 and
independent 969 lists both reduce to exactly 292 orbits, their orbit-key sets
are identical, and the common orbit-key SHA-256 is
`97f630bd464ec01d77087019c5dd02b549f0cbaebd466d5e621ae065797dc848`.
The unequal literal representative counts are normalization-dependent; exact
equality of the genuine orbit sets is the completeness invariant.

Direct blocker dualization leaves 941 families whose 114 minimal
transversals are all order five.  Only the following twelve are ridge
pseudomanifolds with histogram `{2:513}`:

```text
6, 27, 43, 57, 132, 139, 436, 439, 481, 574, 626, 629.
```

The other 84 families have an impure blocker, and the other 929 pure families
fail the ridge law.  No full topology claim is needed for these necessary
rejections.

Every one of the twelve ridge candidates has zero terminal state in the
uncapped complete A5-Lanner/actual-vertex join.  The successive
`support:passing-state` transcripts are

| family | transcript |
|---:|:---|
| 6 | `0:420, 1:96, 2:0` |
| 27 | `0:420, 1:96, 2:0` |
| 43 | `0:420, 1:96, 3:144, 2:2,304, 4:0` |
| 57 | `0:420, 1:96, 5:30, 2:72, 3:144, 4:0` |
| 132 | `0:420, 1:96, 3:0` |
| 139 | `0:420, 1:96, 2:1,152, 3:0` |
| 436 | `0:420, 1:96, 5:0` |
| 439 | `0:420, 1:96, 5:144, 2:2,304, 3:3,456, 4:0` |
| 481 | `0:420, 1:96, 5:0` |
| 574 | `0:420, 3:96, 5:0` |
| 626 | `0:420, 5:96, 2:1,152, 1:2,304, 3:0` |
| 629 | `0:420, 4:96, 5:0` |

Because the terminals are zero, assigning arbitrary additional non-right
finite labels cannot repair a failed elliptic actual-vertex restriction.

## Candidate conclusion and remaining promotion step

The production chain covers `a=0,...,5`: strict link sums remove `a=4,5`;
`a=3` is empty; `a=2` reduces to the existing selected-nerve no-go; the two
ridge-valid `a=1` targets fail unrestricted cross-edge gates; and all twelve
ridge-valid `a=0` targets fail unrestricted sparse gates.  If a structurally
independent implementation reproduces the lossless orbit reductions, target
families, and label transcripts, this proves that the whole record-32 graph
shadow is compact-Coxeter impossible.

Until that replay is green, the formal status remains **candidate proof / no
theorem promotion**.

Production files are

- `scripts/explore_h9_14_record32_all_completions.py`;
- `scripts/explore_h9_14_record32_a0.cpp`;
- `scripts/explore_h9_14_record32_a0_target_gate.py`;
- `scripts/explore_h9_14_record32_a1_completions.py`;
- `scripts/explore_h9_14_record32_a1_cross_gate.py`;
- `scripts/explore_h9_14_record32_a2_additions.py`;
- `scripts/explore_h9_14_record32_a3_additions.py`;
- `scripts/independent_h9_14_record32_a0_replay.cpp`;
- `scripts/compare_independent_h9_14_record32_a0_targets.py`;
- `scripts/certify_h9_14_record32_all_completions_candidate.py`;
- `data/enumerations/h9_14_record32_a0_all_targets_scout.json`;
- `data/enumerations/h9_14_record32_whole_shadow_candidate.json`; and
- `tests/test_h9_14_record32_whole_shadow_scout.py`.

Representative reproduction commands are

```sh
PYTHONPATH=src:. python3 scripts/explore_h9_14_record32_all_completions.py --positive
c++ -std=c++20 -O3 -pthread scripts/explore_h9_14_record32_a0.cpp -o /tmp/h9_14_record32_a0
/tmp/h9_14_record32_a0 --count-only
PYTHONPATH=src:. python3 scripts/explore_h9_14_record32_a0_target_gate.py data/enumerations/h9_14_record32_a0_all_targets_scout.json --summary-only
c++ -std=c++20 -O3 -pthread scripts/independent_h9_14_record32_a0_replay.cpp -o /tmp/independent_record32_a0
/tmp/independent_record32_a0 --target-output /tmp/independent_record32_a0_targets.txt
python3 scripts/compare_independent_h9_14_record32_a0_targets.py data/enumerations/h9_14_record32_a0_all_targets_scout.json /tmp/independent_record32_a0_targets.txt
PYTHONPATH=src:. python3 scripts/explore_h9_14_record32_a1_completions.py
PYTHONPATH=src:. python3 scripts/explore_h9_14_record32_a1_cross_gate.py
PYTHONPATH=src:. python3 scripts/explore_h9_14_record32_a2_additions.py
PYTHONPATH=src:. python3 scripts/explore_h9_14_record32_a3_additions.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record32_all_completions_candidate.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_record32_whole_shadow_scout.py
```
