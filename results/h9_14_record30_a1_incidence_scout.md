# Record 30, `a=1`: exploratory incidence census and compact-label gate (superseded)

This exploratory report has been superseded by the independently replayed
compact theorem and frozen certificate in
[`h9_14_record30_a1_compact_schedule_no_go.md`](h9_14_record30_a1_compact_schedule_no_go.md).
The historical scout below is retained to record the discovery-stage scope.

Status: an exact nauty-quotiented discovery search leaves three necessary-
incidence target orbits.  A separate no-nauty program expands them to all 168
labelled families and replays their literal blockers, face tables, ridges, and
necessary links.  Complete order-four/order-five Lanner joins kill all three
displayed targets.  Search exhaustion has not yet been reproduced by a second
enumerator, so this is not promoted to a whole-`a=1` theorem.

## Corrected schedule bookkeeping

The earlier proposed shadow labels had the first and third representatives
reversed.  Literal independent five-set supersets give

| fixed missing `A4` | forced nonminimal missing five-sets | additional minimal `A5`s |
|---|---:|---:|
| `{3,4,5,6}` | 10 | 4 |
| `{0,3,4,5}` | 9 | 5 |
| `{0,1,3,4}` | 8 | 6 |

Thus every branch has 14 missing independent five-sets, not fourteen minimal
`A5` nonfaces.  The fixed `A4` makes its five-set supersets nonminimal.  The
scheduled face vector is

```text
(14,88,328,805,1358,1582,1232,585,130).
```

## A strict-link exclusion before search

The two-active representative `{0,1,3,4}` is already impossible.  In the
link of the face `{0,1}`, the fixed structure gives

```text
f-prefix = (1,10,43,104)
h        = (1,3,4,4)
g        = (2,1,0)       # g1,g2,g3
```

This is the boundary of a putative simplicial seven-polytope.  Section 3 of
[`h9_14_local_link_constraints.md`](h9_14_local_link_constraints.md) applies:
Murai--Nevo equality at `g3=0` forces a link missing face of order at least
`D-r+2=6`, whereas every link minimal nonface inherits the Coxeter/Lanner
order bound at most five.  Therefore `g3` must be positive.  Additional
minimal `A5`s containing `{0,1}` only decrease this entry and cannot repair
the contradiction.

The older positive-A4 scout did not see this because its post-shadow strict
gate audited vertex links and their aggregate incidence only; it did not
recompute residual edge-link `g3`.

## Quotient census

The remaining branches use the complete fixed-structure stabilizer, a raw
canonical colored-incidence graph rather than a digest, literal transversals
of orders at most four, monotone order-five blocker intersections, and all
430 vertex/edge/triple link incidence caps.

| fixed `A4` | stabilizer | raw addition choices | nodes | scheduled blocker orbits | final necessary-incidence orbits | labelled expansions |
|---|---:|---:|---:|---:|---:|---:|
| `{3,4,5,6}` | 27,648 | 69,512,943,340 | 697,960 | 5 | 1 | 48 |
| `{0,3,4,5}` | 5,760 | 3,741,224,212,670 | 2,402,884 | 64 | 2 | 120 |
| `{0,1,3,4}` | 5,760 | excluded by strict link | 0 | 0 | 0 | 0 |

For the isolate branch, all five scheduled terminals have a literal pure
130-blocker, but four fail the ridge-degree law.  For the one-active branch,
the terminal blocker histogram is `130:59, 131:5`; among the 59 pure
130-blockers, 37 fail the exact face table and 20 fail the ridge gate, leaving
two.

The three surviving addition representatives are

```text
A4 {3,4,5,6}:
  {1,7,8,9,10} {4,7,8,9,10} {7,8,9,10,11} {0,3,5,6,12}

A4 {0,3,4,5}, target 1:
  {3,4,5,7,8} {6,7,8,9,10} {1,6,9,10,11}
  {6,8,9,10,11} {3,4,5,7,12}

A4 {0,3,4,5}, target 2:
  {3,4,5,7,8} {6,7,8,9,10} {2,6,9,10,11}
  {6,9,10,11,12} {1,3,4,5,13}
```

Every survivor has exactly 130 actual rank-nine vertices, the scheduled face
table, 585 degree-two ridges, connected dual graph, all local strict
`g`/Macaulay inequalities, and the correct Euler characteristic in all 5,993
face links.  These are necessary tests, not proofs of PL-sphericity or convex
polytopality.

## No-nauty literal replay

An implementation with no nauty dependency explicitly enumerates
`48 * 8! = 1,935,360` graph actions, obtains stabilizers

```text
27648, 5760, 5760,
```

and target orbit sizes

```text
48, 80, 40.
```

All 168 labelled expansions independently reproduce the literal pure
130-blocker, exact face table, and connected degree-two ridge table.  The
three representatives independently pass all local `g`/Macaulay and
every-face Euler tests.  This validates the emitted target orbits; it does not
independently establish that the quotient search emitted every target.

## Compact Coxeter-label gate

The standalone label replay regenerates all 108 labelled order-four and all
420 labelled order-five Lanner placements and tests them on the literal 130
rank-nine vertices of each target.  The passing-state transcripts are

```text
target 0: 108 -> 84 -> 0
target 1: 108 -> 84 -> 36 -> 108 -> 0
target 2: 108 -> 84 -> 0
```

Thus each displayed target is compact-Coxeter impossible.  Unspecified
meeting pairs are assigned the minimal label `2`; once an actual-vertex
principal diagram is nonfinite, increasing an unspecified Coxeter label
enlarges the corresponding nonnegative cosine adjacency matrix and cannot
restore positive definiteness.  No connected-cross or rank-ten Gram search
is needed because no finite-label state survives.

## Artifacts and exact scope

- `scripts/explore_h9_14_record30_a1_nauty.cpp`
- `scripts/independent_h9_14_record30_a1_targets.cpp`
- `scripts/independent_h9_14_record30_a1_labels.py`
- `data/enumerations/h9_14_record30_a1_scout.json`

The discovery-stage result was compact-only.  Its then-missing distinct replay
of search exhaustion has now been supplied by
`scripts/independent_h9_14_record30_a1_exhaustive.cpp`; see the superseding
theorem report linked at the top.  The resulting theorem still says nothing
about ideal vertices or noncompact finite volume.
