# Record 30, `a=1`: compact schedule no-go

## Theorem certified

There is no compact hyperbolic Coxeter polytope in dimension nine with
fourteen facets, Record-30 face vector

```text
(14,88,328,805,1358,1582,1232,585,130),
```

and exactly one minimal missing four-set (`a=1`).  This is a compact-only
statement.  It neither excludes noncompact finite-volume chambers nor
produces a new hyperbolic Coxeter polytope.

The frozen certificate is
[`h9_14_record30_a1_schedule_no_go.json`](../data/enumerations/h9_14_record30_a1_schedule_no_go.json).

## Correct schedule and the three fixed-`A4` orbits

The missing-edge graph is `3K2`.  Its full action has exactly three orbits on
locally admissible missing four-sets.  A missing `A4` makes each independent
five-set containing it nonminimal, so the scheduled count fourteen is a
count of *all missing independent five-sets*, not fourteen minimal `A5`s.

| fixed missing `A4` | forced nonminimal five-set shadow | minimal `A5`s still to add |
|---|---:|---:|
| `{3,4,5,6}` | 10 | 4 |
| `{0,3,4,5}` | 9 | 5 |
| `{0,1,3,4}` | 8 | 6 |

## One orbit dies analytically

For the two-active representative `{0,1,3,4}`, the residual link of the face
`{0,1}` has

```text
f-prefix = (1,10,43,104)
h        = (1,3,4,4)
g        = (2,1,0)       # g1,g2,g3
```

This would be the boundary of a simplicial seven-polytope.  The strict local
link theorem in [`h9_14_local_link_constraints.md`](h9_14_local_link_constraints.md)
forces `g3>0`: equality `g3=0` would, by the Murai--Nevo equality theorem,
force a missing face of order at least six, while every link minimal nonface
inherits the Coxeter/Lanner order bound at most five.  Further minimal `A5`s
containing `{0,1}` only lower `g3`, so this branch cannot be repaired.

## Two independent exhaustive incidence searches

The other two branches were traversed twice.

The discovery implementation uses nauty on a raw colored incidence graph.
The independent implementation has no nauty or SAT/SMT dependency.  It
canonizes each partial family as follows:

1. enumerate the `3!` orders of the three missing-edge rows;
2. retain the unique missing-`A4` row as its own color;
3. partition added `A5` rows by their exact one-row canonical form and
   enumerate every row order inside each invariant class; and
4. sort the fourteen column signatures.

This is an exact colored-incidence canonical form.  Any isomorphism must send
an added row to a row with the same one-row form, and all exchanges inside
each such class are enumerated.  Conversely, equality of the sorted column
signatures reconstructs a color-preserving row/column isomorphism.

Both traversals use literal order-at-most-four transversals, monotone
order-five cover intersections, all 430 vertex/edge/triple link caps, and at
terminals the complete `2^14` blocker, exact face table, ridge-degree and dual
connectivity equations, local strict-`g`/Macaulay inequalities, and every
face-link Euler equation.

Their complete transcripts agree:

| fixed `A4` | canonical states by depth | terminals | terminal disposition | incidence targets |
|---|---|---:|---|---:|
| `{3,4,5,6}` | `1,11,256,9229,309175` | 5 | five 130-blockers; four ridge failures | 1 |
| `{0,3,4,5}` | `1,9,153,4123,159275,1399930` | 64 | `130:59, 131:5`; then 37 face and 20 ridge failures | 2 |
| `{0,1,3,4}` | strict-link theorem | 0 | analytically impossible | 0 |

The three target orbit sizes under the fixed-`A4` stabilizers are

```text
48, 80, 40,
```

for 168 labelled families in total.  The certificate expands the production
and independent representatives and compares the literal Python sets before
computing any digest.  They are exactly equal.  The combined literal-set
SHA-256

```text
628d18323217fba80cea165f79e2ff047f9bc9ed06c06aa62265d89c5a43fa07
```

is only a compact reproducibility checksum of that prior exact equality; it
is not the mathematical meaning of “hash” and is not used as a substitute
for comparing the families.

## Complete compact Coxeter-label obstruction

For this precise compact dimension-nine/fourteen-facet problem, Burcroff's
Lemma 10.10 excludes multi-multiple edges globally.  Hence the finite meeting
alphabet is exactly `{2,3,4,5}`; there is no unbounded `I2(m)` tail hidden from
the enumeration.

The standalone label replay regenerates all 108 labelled order-four Lanner
placements and all 420 labelled order-five Lanner placements.  On the exact
130 rank-nine vertices of the three incidence targets, its state counts are

```text
target 0: 108 -> 84 -> 0
target 1: 108 -> 84 -> 36 -> 108 -> 0
target 2: 108 -> 84 -> 0
```

Thus no target admits even a finite labeling on every actual vertex.  Every
unspecified meeting label is initially set to `2`.  Increasing such a label
increases the corresponding entry of the nonnegative cosine adjacency
matrix, so Perron--Frobenius monotonicity cannot restore positive definiteness
once an actual-vertex principal diagram is nonfinite.  Hence the sparse gate
covers all larger labels as well.

The three fixed-`A4` orbits are exhaustive, one is analytically impossible,
the two incidence searches agree exactly on all survivors of the other two,
and every survivor is label-impossible.  This proves the stated compact
`a=1` no-go.

## Reproduction and scope

Relevant sources are:

- `scripts/explore_h9_14_record30_a1_nauty.cpp`
- `scripts/independent_h9_14_record30_a1_exhaustive.cpp`
- `scripts/independent_h9_14_record30_a1_labels.py`
- `scripts/certify_h9_14_record30_a1_schedule_no_go.py`
- `tests/test_h9_14_record30_a1_schedule_no_go.py`

Verify the frozen archive with

```text
python3 scripts/certify_h9_14_record30_a1_schedule_no_go.py
python3 -m pytest -q tests/test_h9_14_record30_a1_schedule_no_go.py
```

Use `--recompute --write` on the certifier to compile and rerun both exhaustive
incidence programs and the label replay.  The theorem does not address ideal
vertices or finite-volume noncompact polytopes, and it does not by itself
exclude the other Record-30 values of `a`.
