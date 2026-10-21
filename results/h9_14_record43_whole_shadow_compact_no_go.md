# Record 43: whole-shadow compact-Coxeter no-go

## Theorem and scope

Under the standard compact `H9/14` missing-face reductions, **record 43 has
exactly two residual incidence orbits, and neither can be a compact
hyperbolic Coxeter nerve**.  Thus the entire record-43 missing-edge graph
shadow is excluded on the compact track.

This is not a finite-volume noncompact theorem.  No ideal-vertex,
maximal-parabolic, or Vinberg endpoint-extension census is performed.  The
argument uses strict ellipticity at every actual vertex and exact inertia of
principal Gram blocks.  Its conclusion is therefore deliberately
`compact-only`; it says nothing about whether a noncompact chamber has finite
or infinite volume.

## All four schedules

Write `a` for the number of missing four-subsets and `b` for the total number
of graph-independent five-subsets that are nonfaces, including forced upper
shadows of missing A4s.  Dehn--Sommerville gives

| `a` | `b` | facets | sum of the fourteen vertex-link `g4` values |
|---:|---:|---:|---:|
| 0 | 4 | 90 | 24 |
| 1 | 10 | 88 | 18 |
| 2 | 16 | 86 | 12 |
| 3 | 22 | 84 | 6 |

Every vertex link has strict `g4 >= 1`, so the sum is at least fourteen.
This excludes `a=2,3` and leaves `a=0,1`.

## Exact exclusion of `a=1`

There are 690 graph-independent A4s.  Exact vertex-, edge-, and triple-link
intervals leave 378.  The forced-five-shadow and aggregate Macaulay gate has
the complete histogram

| strict and aggregate feasible | forced A5 shadow | additional minimal A5s | A4 count |
|---|---:|---:|---:|
| no | 8 | 2 | 168 |
| no | 9 | 1 | 112 |
| no | 10 | 0 | 70 |
| yes | 8 | 2 | 28 |

The 28 survivors are exactly

```text
{2,11,i,j},  with {i,j} contained in {3,4,5,6,7,8,9,10}.
```

They form one orbit under the full `S8` isolate action.  For the
representative `(2,3,4,11)`, the verifier constructs all 1,084 possible
additional minimal A5s.  Of the `C(1084,2) = 586,986` pairs, exactly 540
satisfy every vertex-link incidence interval.  Every one of those 540 pairs
has a transversal of order four.  Such a transversal would be a forbidden
minimal nonface of order at most four, so the entire `a=1` branch is empty.

## Exhaustive `a=0` census

The active vertices are `0,1,2,11,12,13`; the other eight vertices are graph
isolates.  The four minimal graph covers are

```text
(0,1,2), (0,1,13), (0,12,13), (11,12,13).
```

Four missing A5s are encoded as four active rows plus a multiset of eight
four-bit isolate-incidence columns.  Exact integer recursion enumerates every
column multiset, rejects repeated rows and small transversals, and counts all
five-element blocker members.  The production order gives

```text
1,073 enumeration paths
  887 exact four-mask clutters
  394 canonical S4-row x S8-isolate states.
```

The number 394 is explicitly not called a full graph-isomorphism count: the
main census deliberately does not quotient active-graph automorphisms.  For
the scheduled 90-facet rows, however, the verifier separately enumerates the
complete active missing-graph automorphism group.  It has order two, and the
two target states remain inequivalent after that action, all row
permutations, and all isolate permutations.  Hence there are exactly two
residual full-graph incidence orbits.

The first A5 family is the frozen selected record-43 nerve:

```text
(3,4,5,6,7)
(1,2,8,9,10)
(2,8,9,10,11)
(8,9,10,11,12).
```

The second is

```text
(1,3,4,5,6)
(2,7,8,9,10)
(3,4,5,6,11)
(7,8,9,10,12).
```

Both blockers have exactly ninety minimal transversals, all of order five.
Both generated complexes have face vector

```text
(14,85,300,690,1088,1190,880,405,90),
```

all 405 ridges of degree two, connected dual graph, spherical Euler
characteristic in every link, and mod-two homology-sphere links.  The first
facet-complement / packed-face-table digests are

```text
e93bd0a0353dc03b18ba499454c08f7f267637a18ccb0e953d22d8faec66effd
32cb562b346841aa4c4a850a911a04da96d907cf6503ed6947f4960bd0212b51
```

and the corresponding second-state digests are

```text
7afc91cca45c5d5796350082aa06d63e0df593734b2fb3fa200c615b6a20faa9
13c95733594fcb5f91ef74a7c261cbc75683f42283884696e574616107dd8b2b
```

## First-orbit Coxeter obstruction

The frozen selected-nerve theorem has now been re-audited for scope and treats
the first orbit with three complete disjoint-A5 support gates.  The
normalizer of every fixed left support has the full order `120` induced
action.  Each gate has 52,500 label-three trials;
labels four and five have no survivors, all 1,890,000 raw two-edge trials are
covered, and the two exceptional weighted-block types have exact inertia
`(8,2,0)`.  This is an unrestricted obstruction: entries incident to nodes
outside a fixed ten-node principal block cannot change that block or its
inertia, so Cauchy interlacing preserves the contradiction.  Thus the first
orbit has no compact Coxeter labeling.

## Second-orbit Coxeter obstruction

The second orbit has the disjoint A5 supports

```text
L = (1,3,4,5,6),   R = (2,7,8,9,10).
```

Vertex `1` is distinguished inside `L`, and vertex `2` is distinguished
inside `R`.  The four remaining vertices on each side are interchangeable.
The verifier explicitly checks all `24 x 24 = 576` elements of

```text
S4({3,4,5,6}) x S4({7,8,9,10})
```

against the entire ten-member minimal-nonface clutter.  On each support, all
420 labelled A5 Lanner placements split into twenty orbits:

```text
5 orbits of size 12 + 15 orbits of size 24 = 420.
```

Consequently the exact representative census has `20 x 20 x 25 = 10,000`
label-three cross-edge trials.  Fifty-five orbit blocks survive local
ellipticity.  Replacing their cross edge by label four or five leaves zero
survivors.

The two-edge statement uses monotonicity carefully.  If a locally elliptic
state had two or more non-right cross edges, deleting all but two and lowering
those two labels to three could only improve positive definiteness of every
finite vertex block.  Thus both chosen edges must first occur among the
single-edge label-three survivors for the same pair of support states.  The
complete representative search covers

```text
20 x 20 x C(25,2) = 120,000
```

raw pairs; after the single-edge gate there are zero residual pairs.  Hence
at most one cross edge can occur.

With no cross edge, the two Lanner blocks form an orthogonal principal sum
with negative index two.  With one cross edge, the 55 exceptional orbit
blocks fall into two weighted-graph classes, of multiplicities 25 and 30.
Exact leading principal minors give sign sequence `++++-----+` and inertia
`(8,2,0)` for both.  Cauchy interlacing contradicts the single negative
direction of an ambient hyperbolic Gram matrix.  The second orbit is therefore
also impossible in the compact case.

### Withdrawn shortcut

A preliminary exploratory calculation used five canonical left Lanner types
as though the full `S5` of the left support acted on the nerve.  That shortcut
was invalid because vertex `1` is distinguished by the missing-edge graph.
It is not used in the theorem, certificate, or tests.  The final proof uses
the verified `S4 x S4` subgroup and the complete twenty-orbit partition on
both supports described above.

## Certificates and reproduction

The theorem is implemented by

- `src/coxeter_polytopes/combinatorics/h9_14_record43_all_completions.py`;
- `src/coxeter_polytopes/combinatorics/h9_14_record43_second_nerve_coxeter_no_go.py`;
- `scripts/certify_h9_14_record43_all_completions.py`;
- `scripts/independent_h9_14_record43_all_completions_audit.py`;
- `data/enumerations/h9_14_record43_all_completions.json`; and
- `tests/test_h9_14_record43_all_completions.py`.

Reproduce it with

```sh
PYTHONPATH=src:. python3 scripts/certify_h9_14_record43_all_completions.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record43_all_completions_audit.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_record43_all_completions.py
```

The hostile replay imports no project module.  It reconstructs the schedules,
the `a=1` and `a=0` censuses, both blockers and face tables, the active-graph
target inequivalence, the distinguished-support Lanner orbit partition, all
finite-label gates, and the exact inertia witnesses.

The frozen certificate digest is
`ca83707489c65f23711849417a28bf548cb5d456a80b9f075e0e642cec785530`;
the independent hostile-summary digest is
`3193e0570170f2e4896b27f5404bee1d23f2240c47dd6ed62f42830254223d61`.
