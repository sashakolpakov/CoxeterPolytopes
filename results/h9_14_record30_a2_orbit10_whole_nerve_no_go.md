# Whole-nerve compact no-go for Record 30 `a=2` orbit 10

> **Later schedule closure.** This whole-nerve theorem is now a subcase of the
> [complete Record-30 compact `a=2` schedule no-go](h9_14_record30_a2_compact_schedule_no_go.md),
> which exhausts all 21 incidence orbits and 77 targets.  The subsequent
> [complete `a=1` theorem](h9_14_record30_a1_compact_schedule_no_go.md) leaves
> only compact schedule `a=0` open; finite-volume noncompact cases remain
> untouched.

## Theorem

The Record 30, `a=2`, positive-`A4` orbit-10 nerve admits no compact
Coxeter realization in \(\mathbb H^9\).

This is a whole-nerve result: all finite labels allowed by the target
\(\{2,3,4,5\}\) alphabet are exhausted, and the three dotted entries remain
arbitrary throughout the decisive rank test.  By itself this certificate is
not a finite-volume noncompact theorem and does not close the rest of the
Record 30 `a=2` schedule; the later theorem linked above supplies that compact
schedule closure.

## Nerve

The inclusion-minimal nonfaces are

```text
dotted: (0,11), (1,12), (2,13)
A4:     (0,3,4,5), (6,7,8,11)
A5:     (3,4,5,9,10), (6,7,8,9,10).
```

There are 128 actual rank-nine vertices.  Exhaustive placement of all 108
order-four and all 420 order-five classical Lanner diagrams, with overlap and
actual-vertex checks after every join, gives

```text
420 -> 84 -> 216 -> 648
```

sparse states.  Their canonical digest is

```text
611d1a2403863268b5410a79b034d729ee9c3fcbc63f197bcc561b6a80111b34.
```

The independent replay uses the opposite join order, its own Lanner table
expansion, and a standalone finite Coxeter forest recognizer.  It recovers
the same 648-state set and digest.

## All residual internal labels are forced right

Inside the ten-node main support component there are nineteen meeting pairs
not already fixed by a Lanner support:

```text
(0,6),(0,7),(0,8),(0,9),(0,10),
(3,6),(3,7),(3,8),(3,11),
(4,6),(4,7),(4,8),(4,11),
(5,6),(5,7),(5,8),(5,11),(9,11),(10,11).
```

For every sparse state, each pair is independently assigned each label
\(3,4,5\), and every affected actual vertex is retested.  All

\[
648\cdot19\cdot3=36,936
\]

trials fail.  This single-addition test is unrestricted by monotonicity: in
normalized form an elliptic Coxeter Gram is \(I-A\), with \(A\) nonnegative;
increasing a right-angle entry to a non-right entry cannot lower the Perron
root and restore positive definiteness.  Consequently all nineteen residual
internal pairs are forced to have label 2 in any compact realization.

## Complete cross-component CSP

The minimal-nonface hypergraph has three components

```text
M = {0,3,4,5,6,7,8,9,10,11},
P = {1,12},
Q = {2,13}.
```

If their finite non-right connection graph is disconnected, the Gram matrix
is block diagonal into at least two blocks with negative directions.  The
main block contains a Lanner principal submatrix; each two-node dotted block
also has one negative direction because its dotted entry has absolute value
greater than the diagonal.  Such a Gram has negative index at least two, not
the required one.  Hence a candidate must connect \(M,P,Q\) and needs at
least two cross edges.

There are 44 cross-component meeting pairs.  The exact single-edge audit has
\(648\cdot44\cdot3=85,536\) literal trials.  Its passing set is identical for
all 648 sparse states and consists of only twenty labelled choices:

- the twelve edges of \(K_{2,2,2}\) on
  \(\{0,11\}|\{1,12\}|\{2,13\}\);
- all eight choices incident to the main dotted pair have label 3; and
- the four edges between the other two dotted pairs may have label 3, 4, or
  5.

No facet 3 through 10 occurs in a locally admissible cross edge.

Starting with every connected two-edge spanning tree and adding every
compatible labelled edge gives the complete per-state shells

```text
non-right cross edges:       2    3    4    5
locally finite assignments: 56  128   38    0.
```

The literal totals over all sparse states are therefore

```text
k=2: 36,288
k=3: 82,944
k=4: 24,624.
```

The empty five-edge shell also excludes every larger shell: any connected
assignment with at least five edges contains a connected five-edge subset,
and local ellipticity is downward closed on deleting non-right edges.

The full nerve automorphism group has order 1,152 and partitions the 648
sparse states into six orbits of sizes

```text
72, 144, 144, 72, 144, 72.
```

The production and hostile computations agree on the literal shell counts
and all three literal assignment-set digests.

## Exact Gram-rank obstruction

For each assignment in each of the six sparse-state orbits, examine the eight
maximal principal blocks that avoid all three dotted entries.  Scale twice
the normalized Gram once more, so

```text
diagonal 4; m=2: 0; m=3: -2;
m=4: -2*sqrt(2); m=5: -(1+sqrt(5)).
```

Every entry lies in \(\mathbb Z[\sqrt2,\sqrt5]\).  Production specializes

\[
\sqrt2\mapsto439,\qquad\sqrt5\mapsto244
\quad\text{in }\mathbb F_{1009}.
\]

Every one of the \(6(56+128+38)=1,332\) quotient assignments has a nonzero
dotted-free \(11\times11\) principal determinant.  Symmetry transports these
witnesses to all 143,856 literal assignments.  Thus the full Gram rank is at
least 11 for every choice of the dotted entries, while a chamber in
\(\mathbb H^9\) has Gram rank at most 10.

The independent program repeats the determinant calculation with

\[
\sqrt2\mapsto406,\qquad\sqrt5\mapsto474
\quad\text{in }\mathbb F_{1129}
\]

and again obtains no zero determinant.

## Scope

This certificate proves compact-Coxeter impossibility for this entire
orbit-10 nerve.  By itself it does not prove:

- impossibility for the other Record 30 `a=2` targets (now supplied by the
  later schedule theorem);
- impossibility for all of Record 30 or all fourteen-facet nerves; or
- impossibility of a finite-volume noncompact chamber with different ideal
  face data.

No new hyperbolic polytope is produced.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/explore_h9_14_record30_a2_orbit10_cross_supersets.py --verify
PYTHONPATH=src:. python3 scripts/independent_h9_14_record30_a2_orbit10_whole_nerve.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_record30_a2_orbit10_whole_nerve.py
```

The frozen machine certificate is
`data/enumerations/h9_14_record30_a2_orbit10_whole_nerve_no_go.json`.
