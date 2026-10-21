# Record 5: an exact rational ordinary-polytopal completion

## Result

**TOOT.** The selected H9/14 record-5 nerve is the boundary of a rational
simplicial convex 9-polytope on fourteen vertices. Dually, it is the nerve
of an ordinary rational simple 9-polytope with fourteen facets.

This is an ordinary-polytopality theorem, not a hyperbolic Coxeter example.
No Coxeter labels or Lorentzian Gram matrix are supplied, and neither
compactness nor finite hyperbolic volume follows.

## Construction

Start with the cyclic 5-polytope `C(8,5)` in moment-curve order

```text
0, 1, 2, 3, 4, 5, 6, 7.
```

Perform two boundary stellar subdivisions:

1. subdivide the edge `{0,2}` with new vertex `8`;
2. subdivide the edge `{1,3}` with new vertex `9`.

The facet counts are `20 -> 28 -> 34`. Then apply four simplicial wedges,
replacing

```text
0 -> {0,10},  2 -> {2,11},  4 -> {4,12},  8 -> {8,13},
```

and relabel construction vertices by

```text
0:3, 1:0, 2:4, 3:12, 4:5, 5:1, 6:2, 7:6,
8:7, 9:13, 10:8, 11:9, 12:10, 13:11.
```

The expanded minimal-nonface clutter is exactly the record-5 clutter: the
three missing edges

```text
{0,12}, {1,13}, {2,13},
```

the three missing four-faces

```text
{0,2,5,10}, {1,6,7,11}, {3,4,8,9},
```

and the four missing five-faces

```text
{1,3,6,8,12}  {2,4,5,9,10}
{2,5,7,10,11} {5,6,7,10,11}.
```

## Exact audits

The first and second stellar vertices use rational radial-beyond parameters
`1/16` and `1/32`. Exhaustive exact rational determinants recover precisely
the 34 constructed simplicial base facets. The base face vector is

```text
(10,41,84,85,34).
```

After the four wedges, complete blocker enumeration gives 136 minimal
transversals, all of order five. The full selected complex has

```text
f = (14,88,329,813,1385,1631,1282,612,136).
```

This is schedule index 3, with three missing four-faces. All 136 facet
complements lie in the frozen record-5 ridge-support core; every one of the
612 ridges has degree two, the dual graph is connected, and the generic
strict-link, Macaulay, Euler, and mod-two link audits pass. A separate
standard-library verifier reconstructs the rational realization, clutter,
blocker, ridges, and face vector without importing the proof module.

Rational convex polytopality is preserved by stellar subdivision and
simplicial wedge. Polarizing gives a rational simple 9-polytope; its rational
normal fan supplies a projective Q-factorial toric 9-fold whose toric
boundary-stratum dual complex is this nerve. No smoothness assertion is
made.

## Frozen artifacts

- `src/coxeter_polytopes/combinatorics/h9_14_record5.py`
- `scripts/certify_h9_14_record5.py`
- `scripts/independent_h9_14_record5_audit.py`
- `data/enumerations/h9_14_record5_polytopality.json`
- `tests/test_h9_14_record5.py`

The proof-source and certificate SHA-256 digests are respectively

```text
37a9c0b6139bd03f1dddb0f2e6ba0298ef4a0c6566ecb4718d62940bb551e598
6cdd1b7951093bb6765a3d77376e3bbad50d9f470d1752d0abc80eb061cdb015
```

Replay with

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_record5.py
python3 scripts/independent_h9_14_record5_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record5.py
```

The selected-nerve Coxeter label/Gram gate is a separate metric calculation.
Until that is run, neither compact nor finite-volume noncompact hyperbolic
realizability of this selected nerve has been decided.
