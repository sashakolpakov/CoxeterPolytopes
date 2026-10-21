# Record 14: complete compact whole-shadow no-go

## Theorem

No compact hyperbolic Coxeter 9-polytope with fourteen facets has the
Record-14 missing-edge graph

```text
M???????????o?w??
```

All seven Dehn--Sommerville schedules `a=0,1,...,6` are now closed.  This is
a whole graph-shadow theorem, not a statement about one selected ordinary
polytope.

The conclusion is compact-only.  It does not exclude a finite-volume
noncompact chamber with ideal faces, and it uses neither endpoint counts nor
parabolic counts as a volume surrogate.

## The final `a=0` incidence census

Here there is no positive `A4` row, so the source branch is unique.  Its exact
schedule requires eight minimal independent `A5` rows and 110 facet
complements.  The source symmetry is the full group

```text
S9({3,...,11}) x S2({0,1}),  order 725,760.
```

Production Z3 and an independently written CP-SAT model each impose:

- exactly eight selected `A5` rows and all local incidence intervals;
- destruction of every base transversal of order at most four;
- an iff variable for every surviving order-five transversal;
- exactly 110 surviving facets; and
- degree zero or two for every six-set.

Both use the same lossless type of generator-lex slice but independently
derive all variables and clauses.  Every finite group orbit meets the slice,
because its global lexicographic minimum is no larger than every generator
image.  Both solvers terminate and agree literally on all 52 pure targets,
not merely on their count.

## Exact full-source quotient

Each target is an unordered eight-row incidence matrix.  The quotient keeps
the five fixed columns `0,1,2,12,13`, tries the exact twin action on `0,1`,
permutes rows only within equal fixed-column codes, and sorts the remaining
nine isolate-column words.  This is an exact canonical form for the full
`S9 x S2` action.  The 52 generator-slice targets split into exactly three
orbits of sizes

```text
14, 20, 18.
```

The groups cover target indices `0,...,51` once each.  No representative is
discarded by a guessed or partial symmetry.

A hostile downstream replay imports neither incidence solver, the production
canonicalizer, nor the project Coxeter classifier.  It independently recovers
the same three classes by colored bipartite graph isomorphism, regenerates the
classical rank-five Lanner table, and repeats the actual-vertex label gate.

## Unrestricted compact-label extinction

For one representative of each full target orbit, all eight minimal
five-faces are assigned complete labelled order-five Lanner diagrams.  There
are 420 placements for the first support.  After every join, every actual
rank-nine vertex is tested for ellipticity.  The exact state sequences are

```text
420 -> 96 -> 144 -> 16 -> 0
420 -> 96 ->   8 -> 0
420 -> 96 ->   8 -> 0
```

Each candidate has 110 actual rank-nine vertices.  The final zero is
unrestricted over the compact finite alphabet `2,3,4,5`: unassigned meeting
pairs begin at label two, and increasing an entry increases the nonnegative
cosine adjacency matrix entrywise.  Perron--Frobenius monotonicity therefore
cannot repair a failed elliptic principal restriction.  Since no label state
survives, no exact Lorentzian Gram candidate reaches the Gram stage.

## Chaining every schedule

The separately frozen certificates give:

| schedules | terminal obstruction |
|---|---|
| `a=5,6` | strict summed vertex-link `g4` |
| `a=3,4` | exhaustive incidence/blocker UNSAT |
| `a=2` | 36 source orbits, 24 literal targets, 6 target orbits, zero labels |
| `a=1` | 5 source orbits, 19 generator-slice targets, 10 target orbits, zero labels |
| `a=0` | 1 source, 52 generator-slice targets, 3 target orbits, zero labels |

Thus the completed schedule list is empty.  Record 14 moves from the open
compact ledger to the wholly excluded ledger.

## Frozen artifacts

- `scripts/certify_h9_14_record14_a0_whole_shadow_no_go.py`
- `scripts/explore_h9_14_record14_a0_symmetry_z3.py`
- `scripts/independent_h9_14_record14_a0_audit.py`
- `scripts/explore_h9_14_record14_a0_downstream.py`
- `scripts/independent_h9_14_record14_a0_downstream_audit.py`
- `data/enumerations/h9_14_record14_a0_whole_shadow_no_go.json`
- `data/enumerations/h9_14_record14_a0_downstream_independent_audit.json`
- `tests/test_h9_14_record14_a0_whole_shadow_no_go.py`

The frozen internal certificate SHA-256 is
`75fa29552f12864ca668e82550a0d63c9b643d5b114f85a15fd641c4b523cd0e`.
The hostile downstream replay summary SHA-256 is
`e871f807999f07e3e89a213669d34fab8ff717abec74c16fe9a96a61337966c3`.
Reproduce the full certificate with

```console
PYTHONPATH=src:. caffeinate -i python3 \
  scripts/certify_h9_14_record14_a0_whole_shadow_no_go.py --fresh \
  --write data/enumerations/h9_14_record14_a0_whole_shadow_no_go.json
PYTHONPATH=src:. python3 \
  scripts/certify_h9_14_record14_a0_whole_shadow_no_go.py
PYTHONPATH=src:. python3 \
  scripts/independent_h9_14_record14_a0_downstream_audit.py
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_h9_14_record14_a0_whole_shadow_no_go.py
```
