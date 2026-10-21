# Record 32 whole-shadow hostile audit

## Verdict

Every record-32 production gate has been independently reproduced on the
compact-Coxeter track.  In particular, a separate positive-orbit program
exhausts all 602 admissible A4s through triples, and a structurally different
`a=0` generator gives exactly the same full-action target-orbit set.  The
certified branch conclusion is

```text
a=0: 1,025 archived target families -> 941 pure blockers
     -> 12 ridge pseudomanifolds -> 0 compact label states
a=1: four strict-shadow A4 orbits -> 1,3,0,6 target orbits
     -> exactly two ridge pseudomanifolds -> 0 unrestricted states
a=2: six strict-shadow A4-pair orbits -> the one older selected nerve
a=3: one strict-shadow orbit -> no no-small-transversal completion
a=4,5: excluded by the strict compact vertex-link g4 sum
```

This promotes only a **whole-graph-shadow compact-Coxeter no-go**.  It does
not assert that the abstract incidence targets are PL spheres or ordinary
polytopes, and it does not exclude a finite-volume noncompact chamber.

## Independent implementation boundary

All hostile Python programs import only the standard library, and the two
hostile C++ programs include no project header or generated producer table.
The finite-diagram recognizer in the downstream program is a separately
implemented exact `A/B/D/E/F/H/I2` classifier; it regenerates the five
order-five and nine order-four classical Lanner types and all 420/108 labelled
placements.  Actual vertices are reconstructed directly as the nine-subsets
avoiding the displayed minimal nonfaces.

The independent sources and SHA-256 digests are

```text
2b0d11f800c88f2d6bb588163b1a7fcf601b1f2c7c48f664432990447bcee843  scripts/independent_h9_14_record32_a0_replay.cpp
2fd380cd8449fcbfbfd8b0461b5c60374aecddd80dab8176a5d8c860590089e2  scripts/compare_independent_h9_14_record32_a0_targets.py
b0576bc3bbdeb603a721d1a8906c3cdb27b65ce0abb0a54ca768855bad81aa3d  scripts/independent_h9_14_record32_positive_orbits.cpp
7702ace985a5c1b5010a7de9f56e611c5a783afce12504eaad85f9d7baf7b544  scripts/independent_h9_14_record32_additions_audit.py
70dc0c02182e218d5c1485536c15f34665186d868e1413f2da0b3f7f06dbf25a  scripts/independent_h9_14_record32_downstream_audit.py
```

The additions and downstream programs' complete canonical JSON outputs have
respective SHA-256 digests
`90ac595c643fa7d0958a6824cf81488824e969d57bc0359bd5578753ba03a95c`
and `0d600a61d85c32d82de35cf38d090494a03c080d9c3ce8f0a65e266dad9c27f1`.
The hostile certifier reruns and compares those complete outputs, not just the
headline counts.

The target archive has SHA-256

```text
ae84f9c8dde0f66cd4d8fc8c0b9dee2afab985c7d484f2e052709ad3c3734938
```

## Full positive-orbit replay

The independent C++ program derives all 419 local constraints directly from
the four missing edges.  Every lower endpoint is zero.  It finds 762
graph-independent A4s, 602 individually interval-feasible A4s, and the exact
full-action branch censuses

| `a` | raw | interval-feasible | strict/aggregate viable | orbits | orbit-size histogram |
|---:|---:|---:|---:|---:|:---|
| 1 | 602 | 602 | 322 | 4 | `{28:1,70:1,112:2}` |
| 2 | 180,901 | 172,039 | 2,156 | 6 | `{56:2,112:1,420:1,672:1,840:1}` |
| 3 | 36,180,200 | 30,954,448 | 56 | 1 | `{56:1}` |

It emits every canonical orbit key and size.  The separate additions program
independently computes the keys and literal orbit--stabilizer sizes of the
four, six, and one representatives; the key/size maps are exactly equal.
Their sizes sum to `322`, `2,156`, and `56`, so the representatives used by
the downstream additions gates exhaust the viable ground sets.

The complete strict/aggregate shadow histograms, including every rejected
shadow/addition bucket, are frozen in
`data/enumerations/h9_14_record32_whole_shadow_hostile_audit.json` and are
checked by the certifier against a fresh C++ run.

## Losslessness of the `a=0` quotient

The active missing-edge graph is `P4 disjoint union K2` on
`(0,1,2,11,12,13)`.  Literal enumeration gives all 24 independent active
masks, with size histogram `1,6,11,6` in orders zero through three.  All 24
occur both before and after the active-cover filter.

The six A5 rows are first made unordered by sorting their active parts.  For
that fixed row order, an isolate is exactly a six-bit incidence column, and
the eight isolates are made unordered by taking a nondecreasing column
multiset.  Conversely, every labelled family can first have its rows sorted
and then its isolates renamed to put those columns in nondecreasing order.
This proves surjectivity of the normalization.  Equal active rows can leave
more than one path for a final family, but the literal family set removes
that overcount; it cannot cause an omission.

The six inclusion-minimal graph covers all have order three.  A transversal
of order at most four is therefore a minimal graph cover, possibly extended
by one vertex.  The active-row test checks every active extension.  An
isolate extension exists exactly when one isolate column contains all row
positions missed by that cover, which is exactly the forbidden-column test.
The production literal postcheck also reports zero discrepancy.

The independent staged counts are

```text
24 active masks, histogram {0:1,1:6,2:11,3:6}
22,513 active-row multisets after caps
14,731 after the active extension test
1,976 distinct allowed-column/row-sum problems
```

The production standalone stream reports

```text
179,785,655 column-multiset paths
3,513,588 repeated-row paths
176,272,067 blocker-valid paths
1,263 scheduled-target paths
1,025 literal target families
```

The structurally different hostile C++ generator uses numeric active-mask
order, explicit isolate-column multiplicity recursion, and a literal
small-transversal postcheck.  It independently reports

```text
24 active masks
22,513 active-row multisets after caps
14,731 after the cover filter
179,785,655 column-multiset paths
3,513,588 repeated-row paths
0 literal-small-transversal discrepancies
176,272,067 blocker-valid paths
1,263 scheduled-target paths
969 literal representatives under its different normalization
```

The no-project comparator canonicalizes the production 1,025 and hostile 969
lists under the full `Aut(P4 disjoint union K2) x S8 x S6` action.  Both give
exactly 292 orbits, with zero orbit on either side only, and common orbit-key
digest

```text
97f630bd464ec01d77087019c5dd02b549f0cbaebd466d5e621ae065797dc848
```

Thus the unequal literal representative counts reflect only the deliberately
different quotient sections; the genuine orbit sets agree exactly.

## Hostile blocker and ridge replay

A literal scan of all `2^14-1` nonempty subsets was used for each archived
family.  It recomputes inclusion-minimality after deleting each chosen
vertex; it does not use the producer's five-cover word shortcut.

Exactly 941 families have blocker histogram `{5:114}`.  The other 84 have
the same 114 order-five transversals and additional minimal order-six
transversals, with multiplicities distributed as follows:

```text
{6:1}:3, {6:2}:30, {6:3}:8, {6:4}:20, {6:8}:19, {6:16}:4.
```

Direct facet-to-ridge incidence leaves exactly indices

```text
6, 27, 43, 57, 132, 139, 436, 439, 481, 574, 626, 629
```

with ridge histogram `{2:513}`.  This independently reproduces
`1,025 -> 941 -> 12`.

For those twelve families, a complete 420-placement order-five Lanner join
was replayed with the opposite support tie-breaking from production.  All
twelve terminal state sets are empty.  The different ordering gives a
different intermediate transcript on several rows, but the same uncapped
zero conclusion.  This is the needed unrestricted obstruction: increasing
any unassigned nonnegative cosine adjacency cannot turn a non-positive
actual-vertex Gram restriction positive definite.

## `a=1` replay

The fixed-A4 hostile census independently derives all graph-independent A5
rows, the forced shadow, vertex caps, all independent active row parts, and
both isolate-column multisets.  Its exact transcripts are

| A4 representative | active configurations | column paths | distinct-row paths | no-small paths | literal targets | target orbits |
|:---|---:|---:|---:|---:|---:|---:|
| `(3,4,5,6)` | 58 | 1,161 | 1,102 | 9 | 1 | 1 |
| `(2,3,4,5)` | 70 | 43,229 | 42,865 | 271 | 3 | 3 |
| `(1,2,3,4)` | 155 | 2,323,448 | 2,278,508 | 12,873 | 0 | 0 |
| `(2,3,4,11)` | 410 | 6,068,459 | 5,934,396 | 19,836 | 8 | 6 |

The full facet-count histograms also agree.  Literal blocker/ridge replay
leaves only the two production targets, each with 112 five-transversals and
ridge histogram `{2:504}`.

The independent Lanner join gives 10,368 and 864 sparse states.  For the
chosen disjoint order-four/order-five Lanner supports there are twenty cross
pairs and no dotted pair.  All `207,360` and `17,280` weakened single-`m=3`
tests fail; the first-bad restriction histograms are respectively
`{7:206208,8:1152}` and `{7:17280}`.

The unrestricted argument is valid.  If the cross block is zero, the
principal Gram is the orthogonal sum of two Lanner Grams and has negative
index two.  (Each classical representative has negative determinant and all
proper principal subdiagrams are elliptic.)  Otherwise choose a non-right
cross edge, lower it to `m=3`, and lower every unforced edge to `m=2`.
Perron--Frobenius monotonicity makes that tested adjacency entrywise no
larger than the original, so its failed elliptic restriction cannot be
repaired by stronger or additional labels.

## `a=2` and `a=3`

The six `a=2` representatives were replayed from direct independent-set,
link-g/Macaulay, small-transversal, and five-cover computations.  Their
`(shadow, additions, candidates, raw, cap-feasible, no-small, target)` rows
are

```text
(18,0,  0,        1,        1,   0,0)
(17,1,  1,        1,        1,   1,1)
(15,3,462,16,328,620,1,794,605, 840,0)
(16,2, 56,    1,540,      540,   0,0)
(16,2, 56,    1,540,      945,   0,0)
(16,2, 56,    1,540,    1,540,   0,0)
```

The sole target is exactly the older certified record-32 nerve with A4s
`(2,3,4,5)`, `(3,4,5,11)` and A5 `(6,7,8,9,10)`.  The whole-shadow theorem
therefore chains to its already frozen unrestricted compact-label no-go; it
does not create a new polytopality claim.

For `a=3`, the cap-zero/singleton, cap-one/pair, and cap-two/triple generator
is lossless because every local interval has lower endpoint zero and a
three-member family can violate only an upper cap zero, one, or two.  The
sole viable orbit representative has forced shadow 22 and needs two A5s.
The independent addition replay gives 56 candidates, 1,540 raw pairs, 540
cap-feasible pairs, and zero pairs without a transversal of order at most
four.

## Promoted theorem scope

The safe promoted statement is:

> Every necessary compact completion of H9/14 record 32 in the complete
> Dehn--Sommerville schedules fails either strict compact link positivity,
> blocker purity, the ridge-degree law, or the unrestricted
> Lanner/actual-vertex Coxeter gate.  Hence record 32 is impossible as a
> compact hyperbolic Coxeter graph shadow.

The statement must retain all of the following qualifiers:

- compact only;
- abstract target families are necessary incidence objects, not asserted PL
  spheres or ordinary polytopes;
- ridge-invalid targets need no stronger topology inference;
- the unique `a=2` target uses the older selected-nerve theorem;
- no finite-volume noncompact conclusion follows.

## Reproduction

The quick frozen-source and semantic check is

```sh
python3 scripts/certify_h9_14_record32_whole_shadow_hostile.py
```

The fresh positive-orbit, additions, blocker/ridge/Lanner, and a1 cross replay
is

```sh
python3 scripts/certify_h9_14_record32_whole_shadow_hostile.py --downstream
```

Adding `--full` also reruns the 179,785,655-path independent `a=0` generator
and the exhaustive 292-orbit comparator.  Every mode retains the compact-only
scope and explicitly reports that it makes no finite-volume noncompact
inference.
