# Record 30 `a=2` orbit 10: the complete two-cross-edge shell

## Result

The complete minimal cross-component finite-label shell of the surviving
Record 30, `a=2`, orbit-10 nerve is impossible in \(\mathbb H^9\).

There are 36,288 labelled assignments in this shell which pass every local
Lanner condition and all 128 actual rank-nine vertex tests.  Every one of
them nevertheless has a fixed, dotted-variable-free \(11\times11\) principal
Gram minor which is nonzero.  Its full Gram matrix therefore has rank at
least 11, whereas wall normals in \(\mathbb R^{9,1}\) have Gram rank at most
10.

This is a scoped exact no-go theorem, not a new hyperbolic example.  It closes
**exactly two** cross-component non-right meeting edges; assignments with
three or more such edges remain open.

## Nerve and sparse states

The inclusion-minimal nonfaces are

```text
dotted pairs: (0,11), (1,12), (2,13)
A4 supports:  (0,3,4,5), (6,7,8,11)
A5 supports:  (3,4,5,9,10), (6,7,8,9,10)
```

Their support hypergraph has three components:

```text
M = {0,3,4,5,6,7,8,9,10,11},
P = {1,12},
Q = {2,13}.
```

The production join places every one of the 420 order-five and 108
order-four classical Lanner diagrams, checks overlap consistency, and tests
the induced diagrams on the actual vertices after each join.  Its state
counts are

```text
420 -> 84 -> 216 -> 648.
```

Thus 648 sparse label states survive when every as-yet unassigned meeting
pair is labelled 2.  Their canonical digest is

```text
611d1a2403863268b5410a79b034d729ee9c3fcbc63f197bcc561b6a80111b34
```

The independent program reconstructs the same 648 states using the opposite
support-join order and a standalone finite Coxeter forest recognizer.

## Complete two-edge census

There are 20 meeting pairs between \(M,P\), 20 between \(M,Q\), and four
between \(P,Q\).  Two edges connect all three components precisely when they
belong to two different component-pair classes.  Hence the literal placement
count is

\[
 20\cdot20+20\cdot4+20\cdot4=560.
\]

Each edge independently receives a label in \(\{3,4,5\}\), so the complete
literal census contains

\[
 648\cdot560\cdot3^2=3,265,920
\]

assignments.  A monotonic single-edge pass reduces this to 72,576 interacting
two-edge trials without discarding a possible survivor.  All actual vertices
containing both new edges are then retested exactly.  The outcome is 36,288
local survivors—56 for every sparse state—with label histogram

```text
(3,3): 15,552
(3,4): 10,368
(3,5): 10,368
```

and canonical survivor digest

```text
151a6305767cb86e227efb0b1be7559c519c1f1841897d181e42599c586a6806
```

In particular, the previously displayed state-0 assignment with new
`m=3` edges `(0,1)` and `(2,11)` occurs in this complete census; it is no
longer being used as evidence for the rest of the shell.

## Exact global rank obstruction

Use twice the normalized Coxeter Gram matrix and then multiply once more by
two.  The resulting finite entries are

```text
diagonal 4,  m=2: 0,  m=3: -2,
m=4: -2*sqrt(2),  m=5: -(1+sqrt(5)).
```

They lie in \(\mathbb Z[\sqrt2,\sqrt5]\).  For every local survivor, the
principal block on facets

```text
3,4,5,6,7,8,9,10,11,12,13
```

contains no dotted entry: facets 0, 1, and 2 remove one endpoint of each
dotted pair.  The production verifier specializes

\[
 \sqrt2\mapsto439,\qquad \sqrt5\mapsto244
 \quad\text{in }\mathbb F_{1009};
\]

indeed \(439^2=2\) and \(244^2=5\pmod {1009}\), and 1009 is prime.  The
36,288 determinants fall into fifteen residue classes, none zero.  A nonzero
specialization is an exact proof that the original algebraic-integer
determinant is nonzero.  It follows that every full Gram completion, for all
values of the three dotted entries, has rank at least 11.

The hostile replay uses a separately implemented determinant eliminator and
the different specialization

\[
 \sqrt2\mapsto406,\qquad \sqrt5\mapsto474
 \quad\text{in }\mathbb F_{1129}.
\]

It independently obtains fifteen nonzero residue classes and agrees on all
census counts and both state-set digests.

## Exact scope

The theorem excludes every compact Coxeter labeling of this nerve having
exactly two cross-component finite non-right meeting edges, in the complete
target label alphabet \(\{3,4,5\}\).  It does **not** exclude:

- assignments with three or more cross-component non-right edges;
- the whole orbit-10 nerve;
- the whole Record 30 `a=2` schedule; or
- a finite-volume noncompact chamber with different ideal-face data.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/explore_h9_14_record30_a2_orbit10_two_cross.py --verify
PYTHONPATH=src:. python3 scripts/independent_h9_14_record30_a2_orbit10_two_cross.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_record30_a2_orbit10_two_cross.py
```

The frozen machine-readable certificate is
`data/enumerations/h9_14_record30_a2_orbit10_two_cross.json`.
