# Record 1: an exact rational ordinary-polytopal completion

## Result

**TOOT.** The selected H9/14 record-1 nerve is the boundary of a rational
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

1. subdivide the triangle `{0,1,2}` with new vertex `8`;
2. subdivide the triangle `{0,1,3}` with new vertex `9`.

The facet counts are `20 -> 30 -> 36`. Then apply four simplicial wedges,
replacing

```text
1 -> {1,10},  3 -> {3,11},  4 -> {4,12},  8 -> {8,13},
```

and relabel construction vertices by

```text
0:3, 1:4, 2:0, 3:5, 4:6, 5:1, 6:2, 7:7,
8:8, 9:13, 10:9, 11:10, 12:11, 13:12.
```

The expanded minimal-nonface clutter is exactly the record-1 clutter: the
three missing edges

```text
{0,13}, {1,13}, {2,13},
```

the four missing four-faces

```text
{0,1,3,7} {0,2,6,11} {0,3,4,9} {1,7,8,12},
```

and the ten missing five-faces

```text
{0,3,6,7,11}  {1,3,5,7,10}  {1,4,5,9,10}
{1,5,8,10,12} {2,4,5,9,10}  {2,4,6,9,11}
{2,5,8,10,12} {2,6,8,11,12} {3,4,5,9,10}
{6,7,8,11,12}.
```

## Exact audits

Both stellar vertices use rational radial-beyond parameter `1/16`.
Exhaustive exact rational determinants recover precisely the 36 constructed
simplicial base facets. The base face vector is

```text
(10,42,88,90,36).
```

After the four wedges, complete blocker enumeration gives 152 minimal
transversals, all of order five. The full selected complex has

```text
f = (14,88,331,831,1451,1757,1414,684,152).
```

This is schedule index 4, with four missing four-faces. All 152 facet
complements lie in the frozen record-1 ridge-support core; every one of the
684 ridges has degree two, the dual graph is connected, and the generic
strict-link, Macaulay, Euler, and mod-two link audits pass. A separate
standard-library verifier reconstructs the rational realization, clutter,
blocker, ridges, and face vector without importing the proof module.

Rational convex polytopality is preserved by stellar subdivision and
simplicial wedge. Polarizing gives a rational simple 9-polytope; its rational
normal fan supplies a projective Q-factorial toric 9-fold whose toric
boundary-stratum dual complex is this nerve. No smoothness assertion is
made.

## Frozen artifacts

- `src/coxeter_polytopes/combinatorics/h9_14_record1.py`
- `scripts/certify_h9_14_record1.py`
- `scripts/independent_h9_14_record1_audit.py`
- `data/enumerations/h9_14_record1_polytopality.json`
- `tests/test_h9_14_record1.py`

The proof-source and certificate SHA-256 digests are respectively

```text
d11365f9794a1b2fd450e08f60a5f334b4203b2e48a69ed44578555e3de9a6a3
5684d466a7bca45828da56aae3bd1ea91c5390e08f966c98bccf71b427f2b976
```

Replay with

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_record1.py
python3 scripts/independent_h9_14_record1_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record1.py
```

The selected-nerve Coxeter label/Gram gate is a separate metric calculation.
Until that is run, neither compact nor finite-volume noncompact hyperbolic
realizability of this selected nerve has been decided.
