# Record 18: independent high-addition compact audit

Status: the two `a=2/add5` branches are incidence-impossible.  Of the three
`a=1/add6` branches, two are incidence-impossible and the third has exactly
eight quotient targets of total literal orbit weight 1,080; all eight targets
are topology-valid but admit no compact Coxeter labeling.  This is
compact-only and makes no finite-volume noncompact claim.

## Solver-independent incidence search

The standalone C++ replay uses no SAT/SMT backend and imports no producer
code.  It quotients every partial row-incidence matrix by row permutations,
the exact stabilizer of the prescribed missing four-faces, and all
permutations of the residual isolate columns.  Its canonical state is stored
without truncation in 128 bits; the `a=1` depth-six normal form needs 72 bits.

Mandatory sets drive the recursion:

- an uncovered transversal of order at most four must be killed by a future
  missing five-face;
- if more than the scheduled number of facets survive, a future row must
  delete one;
- a current ridge degree one or greater than two must decrease;
- local upper incidence bounds and an optimistic facet-deletion bound prune
  only impossible continuations.

Every accepted terminal is dualized literally.

The two `a=2/add5` representatives reproduce the base counts

```
candidates 777, small transversals 56,
base five-transversals 233, nonempty ridge equations 630
```

and have zero ridge-scheduled terminals.  Their search sizes are 413,064 and
87,819 states.

The three `a=1/add6` representatives reproduce

```
candidates 784, small transversals 64,
base five-transversals 254, nonempty ridge equations 665.
```

Representatives `15` and `23` have zero terminals.  Representative `53`
has exactly eight canonical targets with orbit sizes

```
20, 40, 60, 120, 180, 180, 240, 240,
```

whose sum is 1,080.

## Independent topology and Coxeter replay

Each of the eight `a=1` targets has

```
f = (14, 86, 310, 734, 1198, 1358, 1036, 486, 108),
blocker sizes = {5: 108},
ridge degrees = {2: 486}.
```

The dual graph is connected, the reconstructed minimal nonfaces are exactly
the five missing edges, one missing four-face, and six added missing
five-faces, and every face link has the exact mod-two homology of a sphere.
The strict local `g` and Macaulay gates also pass.

The independent label replay regenerates all 108 placements of the nine
rank-four Lanner types and all 420 placements of the five rank-five Lanner
types, then checks every induced restriction at all 108 actual vertices.
Every class terminates with zero states.  The stage sequences, paired with
literal orbit size, are

```
20  : 108 -> 0
40  : 108 -> 0
60  : 108 -> 84 -> 8 -> 0
120 : 108 -> 84 -> 8 -> 0
180 : 108 -> 84 -> 8 -> 0
180 : 108 -> 84 -> 24 -> 16 -> 72 -> 0
240 : 108 -> 84 -> 8 -> 0
240 : 108 -> 84 -> 24 -> 16 -> 0.
```

Thus none of the 1,080 candidate incidence objects is a compact hyperbolic
Coxeter polytope.

## Artifacts

- `scripts/independent_h9_14_record18_a2_add5_audit.cpp`
- `scripts/certify_independent_h9_14_record18_high_addition_audit.py`
- `scripts/independent_h9_14_record18_a1_downstream_audit.py`
- `data/enumerations/h9_14_record18_high_addition_independent_audit.json`
- `data/enumerations/h9_14_record18_a1_downstream_independent_audit.json`
- `tests/test_h9_14_record18_high_addition_independent_audit.py`

The combined focused independent suite passes: `6 passed`.
