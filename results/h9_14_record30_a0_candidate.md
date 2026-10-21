# Record 30, `a=0`: an abstract multiwedge sphere and its compact-label no-go

Status: one explicit `a=0` completion has been independently reconstructed
and passes every implemented abstract-sphere incidence test.  It is an exact
multiwedge of a ten-vertex simplicial 4-sphere.  Nevertheless, an exhaustive
order-five-Lanner join proves that this completion has no compact Coxeter
labeling in the full reduced alphabet `{2,3,4,5}`.

This is a theorem about this completion only.  The full `a=0` orbit census
is still running, convex realizability of the base sphere is not proved, and
finite-volume noncompact chambers are not considered.

## Explicit minimal nonfaces

Besides the missing edges

```text
{0,11}, {1,12}, {2,13},
```

the completion has exactly these eight minimal missing five-faces:

```text
{0,1,4,6,7}       {0,4,6,7,13}
{1,4,5,6,7}       {4,5,6,7,8}
{2,3,9,10,12}     {3,5,8,9,10}
{3,8,9,10,11}     {3,9,10,11,12}.
```

The standalone replay scans all `2^14` subsets.  The minimal-transversal
clutter is pure `5^132`; its complements generate

```text
f = (14,88,328,806,1364,1596,1248,594,132).
```

All 594 ridge groups have degree two and the dual graph is connected.  The
recovered minimal nonfaces are exactly the three edges and eight five-sets.
The replay also checks all 430 vertex/edge/triple local `g` bases, 874
Macaulay transitions, 6,039 Euler face links, and 6,171 mod-two homology
face links.  Thus this is a strong abstract homology-sphere candidate, not a
certificate of PL or convex polytopality.

## Exact multiwedge description

The two clone classes are

```text
{3,9,10},       {4,6,7}.
```

Collapsing each class to one vertex produces a ten-vertex simplicial
4-sphere with

```text
f_base = (10,42,88,90,36).
```

Expanding the two base vertices with multiplicity three recovers every
minimal nonface above exactly.  This supplies a compact algebro-combinatorial
description of the candidate.  Whether the base sphere is the boundary of a
convex 5-polytope remains open in this certificate.

## Unrestricted compact Coxeter-label obstruction

The replay independently generates all 420 labelled placements of the five
classical order-five Lanner types.  Joining supports in the production order
gives locally elliptic state counts

```text
420 -> 96 -> 6 -> 0.
```

A deliberately different reverse order gives

```text
420 -> 96 -> 0.
```

At each stage, every restriction is the literal intersection with one of the
132 actual rank-nine vertices.  Unspecified meeting pairs are assigned the
minimal label `2`.  Once such a principal diagram is nonelliptic, increasing
an unspecified label to `3`, `4`, or `5` only increases the associated
nonnegative adjacency matrix; Perron--Frobenius monotonicity cannot restore
positive definiteness.  The zero terminal count is therefore an unrestricted
compact-label no-go, not merely a sparse-label failure.

## Artifacts

- `scripts/independent_h9_14_record30_a0_candidate_audit.py`
- `scripts/certify_h9_14_record30_a0_candidate.py`
- `data/enumerations/h9_14_record30_a0_candidate.json`
- `tests/test_h9_14_record30_a0_candidate.py`

Frozen certificate SHA-256:

```text
e40b22cc094901627c31186405faf6bd238dfd0633440677e788d02b5003b5fd
```
