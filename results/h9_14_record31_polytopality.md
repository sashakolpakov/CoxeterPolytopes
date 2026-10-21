# Record 31: an exact rational ordinary-polytopal completion

## Result

**TOOT.** The selected H9/14 record-31 nerve is the boundary of a rational
simplicial convex 9-polytope on fourteen vertices.  Dually, it is the nerve
of an ordinary rational simple 9-polytope with fourteen facets.

This is an ordinary-polytopality theorem, not a hyperbolic Coxeter example.
No Coxeter labels or Lorentzian Gram matrix are supplied, and neither
compactness nor finite hyperbolic volume follows.

## Construction

Start with the cyclic 5-polytope `C(8,5)`, whose moment-curve labels occur in
the order

```text
2, 14, 15, 16, 17, 11, 13, 12.
```

Perform two boundary stellar subdivisions:

1. subdivide the triangle `{2,15,17}` with new vertex `1`;
2. subdivide the triangle `{2,14,16}` with new vertex `0`.

The facet counts are `20 -> 26 -> 34`.  Then apply four simplicial wedges,
replacing

```text
14 -> {3,4},  15 -> {5,6},  16 -> {7,8},  17 -> {9,10}.
```

The expanded minimal-nonface clutter consists exactly of the record-31
missing-edge path

```text
11--0--13--1--12
```

and the nine missing five-faces

```text
{0,5,6,9,10}       {2,3,4,7,8}       {2,5,6,9,10}
{2,5,6,11,12}      {2,7,8,11,12}     {3,4,7,8,11}
{3,4,7,8,13}       {3,4,9,10,13}     {5,6,9,10,13}.
```

## Exact audits

Both stellar vertices use the rational radial-beyond formula

```text
y + (y-p)/32,
```

where `y` is the carrier-face barycentre and `p` is the current vertex
barycentre.  Exhaustive exact rational determinants recover precisely the 34
constructed simplicial base facets.  The base face vector is

```text
(10,41,84,85,34).
```

After the four wedges, complete blocker enumeration gives 122 minimal
transversals, all of order five.  The full selected complex has

```text
f = (14,87,319,771,1287,1491,1158,549,122).
```

All 122 facet complements lie in the frozen record-31 ridge-support core;
every ridge has degree two, the dual graph is connected, and the generic
strict-link, Macaulay, Euler, and mod-two link audits pass.  A separate
standard-library verifier reconstructs the rational realization, clutter,
blocker, ridges, and face vector without importing the proof module.

Rational convex polytopality is preserved by stellar subdivision and
simplicial wedge.  Polarizing gives a rational simple 9-polytope; its rational
normal fan supplies a projective Q-factorial toric 9-fold whose toric
boundary-stratum dual complex is this nerve.  No smoothness assertion is
made.

## Frozen artifacts

- `src/coxeter_polytopes/combinatorics/h9_14_record31.py`
- `scripts/certify_h9_14_record31.py`
- `scripts/independent_h9_14_record31_audit.py`
- `data/enumerations/h9_14_record31_polytopality.json`
- `tests/test_h9_14_record31.py`

The proof-source and certificate SHA-256 digests are respectively

```text
ad30ec2b567c865165d805f1f2375211cf03aec766b574e8fe2fcf6cb987f890
12b01aec7f3de0bb4e71040291753e8362787b9c7063704877bd24693884a799
```

Replay with

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_record31.py
python3 scripts/independent_h9_14_record31_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record31.py
```

The selected-nerve compact-Coxeter label/Gram gate is a separate metric
calculation. It has now excluded this selected nerve in the compact case;
see `results/h9_14_record31_coxeter_label_gate.md`. That later theorem still
does not address finite-volume noncompact realizations or other completions
of graph-shadow 31.
