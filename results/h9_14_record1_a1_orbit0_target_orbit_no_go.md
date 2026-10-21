# Record 1, `a=1`, source-orbit 0: a 60,480-target compact-label no-go

**Status: exact bounded theorem, 2026-08-07.**  This certificate closes one
complete labelled incidence orbit for the first positive-`A4` source of
Record 1.  It does not close the source, the `a=1` schedule, or the whole
shadow.

## The correction

An earlier exploratory command used `--max-models 10` and its stdout was not
archived.  A report subsequently described nine observed pure targets as an
exact census.  That was incorrect: without a terminal `UNSAT`, a model limit
is a cutoff, not a completeness certificate.

There is also a decisive internal contradiction to the old count.  The first
frozen target has 60,480 distinct images under the stabilizer of its fixed
source.  Consequently there are at least 60,480 labelled pure targets, before
considering any other target orbit.

## Exact orbit theorem

The Record-1 missing graph has edges

```text
{0,13}, {1,13}, {2,13},
```

and source orbit 0 is represented by the four-nonface `{0,3,4,5}`.  Graph
degree and source membership partition the vertices as

```text
{13} | {0} | {1,2} | {3,4,5} | {6,7,8,9,10,11,12}.
```

It follows directly that the full fixed-source stabilizer is

```text
S2 x S3 x S7,
```

of order `2! 3! 7! = 60,480`.  Exhaustive streaming over these 60,480 group
elements finds exactly one element fixing the representative nineteen-row
target.  Its target stabilizer is therefore trivial, and orbit--stabilizer
gives a labelled target orbit of size 60,480.  Production and independently
written replays have the same action-image stream hash

```text
3f36c8a0e2cc7501a3b65c17f25c9e73ee09b843df2714a85f9ffe0be6fbecee.
```

The representative itself independently passes the complete literal
incidence tests:

- 19 additional minimal five-nonfaces;
- 158 minimal transversals, all of order five;
- ridge-degree histogram `{0:2292, 2:711}`;
- hence every active six-set ridge has exactly two incident facets.

## Compact-label obstruction

The production Lanner join and a hostile replay using independently tabulated
classical order-four and order-five Lanner diagrams both reconstruct 158
actual rank-nine vertices.  Their complete survivor counts are

```text
108 -> 84 -> 8 -> 0,
```

with state hashes

```text
857c6268544a6e22216a4068a1642a2e48730ff5288ae810cdaa398a9517f730
de2dcf81429b02fe4970c82cc168428caddf019335da066bed4c228047e8b2e5
a3cd5a07b59bf0c8e0b299bd334049f235dc243abded9be5f31d1ad6f011f270
4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945
```

Thus the representative admits no compact Coxeter labelling.  The property is
equivariant under facet relabelling, so every one of its 60,480 labelled
images is excluded as well.

## Exact scope

What is proved:

- one full 60,480-member labelled incidence orbit is literal and pure;
- every target in that orbit has no compact Coxeter labelling.

What is not proved:

- completeness of the target-orbit quotient for source orbit 0;
- absence of a second incidence orbit for that source;
- exclusion of source orbits 1 and 2;
- a Record-1 `a=1` schedule theorem or whole-shadow theorem;
- any finite-volume noncompact conclusion.

The next terminating step is a fresh symmetry-quotiented incidence run to
terminal `UNSAT`, followed by a separately derived quotient replay.  The
enumerator must stream targets/checkpoints so an interruption cannot again be
mistaken for completeness.

## Frozen artifacts and replay

- `data/enumerations/h9_14_record1_a1_orbit0_target_orbit_no_go.json`
- `scripts/certify_h9_14_record1_a1_orbit0_target_orbit.py`
- `scripts/independent_h9_14_record1_a1_orbit0_target_orbit_audit.py`
- `scripts/independent_h9_14_record1_a1_target0_label_audit.py`
- `tests/test_h9_14_record1_a1_orbit0_target_orbit_no_go.py`

```text
PYTHONPATH=src:. python3 \
  scripts/certify_h9_14_record1_a1_orbit0_target_orbit.py
PYTHONPATH=src:. python3 \
  scripts/independent_h9_14_record1_a1_orbit0_target_orbit_audit.py
PYTHONPATH=src:. pytest -q \
  tests/test_h9_14_record1_a1_orbit0_target_orbit_no_go.py
```

Frozen certificate digest:

```text
199a9cbfa94d32cf96371d8dec1dd1411ba6d2dbd1afc5d0d7e4fc2e93fb4eae
```
