# Record 42: whole-shadow compact-Coxeter no-go

## Theorem and exact scope

Under the standard compact `H9/14` missing-face reductions, **no incidence
completion of record 42 is the nerve of a compact hyperbolic Coxeter
9-polytope with fourteen facets**.  This excludes the whole record-42
missing-edge graph shadow, not merely one selected completion.

This is a compact theorem only.  It uses strict ellipticity at actual
vertices and compact Lanner missing faces.  It does not enumerate parabolic
subdiagrams, ideal endpoints, or Vinberg cusp extensions, so it neither
excludes nor constructs a finite-volume noncompact chamber.

Record 42 has graph6 string

```text
M?????????O?o?O@?
```

and missing edges

```text
(0,11), (0,12), (1,12), (1,13), (11,13).
```

They form a `C5` on `(0,11,13,1,12)`; vertices `2,...,10` are isolates.

## Dehn--Sommerville schedules

Write `a` for the number of missing A4s and `b` for the total number of
missing independent five-subsets, including the upper shadows forced by the
A4s.  The five exact schedules are

| `a` | `b` | facets | sum of the fourteen vertex-link `g4` values |
|---:|---:|---:|---:|
| 0 | 5 | 102 | 30 |
| 1 | 11 | 100 | 24 |
| 2 | 17 | 98 | 18 |
| 3 | 23 | 96 | 12 |
| 4 | 29 | 94 | 6 |

Every compact vertex link has strict `g4 >= 1`; hence the sum must be at
least fourteen.  This excludes `a=3,4` and leaves `a=0,1,2` for exact
enumeration.

## The positive branches

There are 546 individually link-feasible A4s.

For `a=2`, all `C(546,2) = 148,785` unordered pairs are considered.  Exact
family-local intervals leave 131,355.  Their forced A5-shadow histogram is

| shadow size | pair count |
|---:|---:|
| 15 | 420 |
| 16 | 70,140 |
| 17 | 2,520 |
| 18 | 50,400 |
| 19 | 1,260 |
| 20 | 6,615 |

There are 58,275 immediate schedule overruns, and the strict local/Macaulay
gate leaves zero pairs.  Thus `a=2` is empty.

For `a=1`, the forced-shadow histogram is `{8:420, 10:126}`.  The strict
gate leaves precisely the 126 choices of four of the nine isolates.  For
each, the last minimal A5 is forced to be the complementary five isolates.
Consequently all survivors form one `S9` orbit, represented by

```text
A4 = (2,3,4,5),    A5 = (6,7,8,9,10).
```

This is exactly the frozen `C5 x Delta3 x Delta4` product nerve.  Its
re-audited selected-support certificate uses the full normalizing support
symmetries of orders `24` and `120`, covering all 108 labelled A4 Lanner
placements and all 420 labelled A5 placements.  It tests 900 single
label-three cross edges, leaving two; labels four and five leave none.  Of
8,550 raw two-edge trials, only one monotonicity residual remains and it also
fails.  The two exceptional blocks are symmetry-equivalent and have exact
inertia `(7,2,0)`.

This is an unrestricted obstruction on the selected nine-node principal
block.  Labels incident to facets outside that block do not change the block
or its inertia, so Cauchy interlacing prevents any outside dotted component
or additional Lanner support from removing its two negative directions.
Hence the unique `a=1` orbit is not a compact Coxeter nerve.

## Corrected exhaustive `a=0` census

Five missing A5s are represented as five rows.  Strict singleton-link
intervals put every active `C5` vertex in at most one row.  The active part of
a row must be an independent subset of the `C5`, and therefore has size at
most two.  There are exactly eleven possible active row masks:

```text
1 empty + 5 singletons + 5 nonadjacent pairs.
```

Taking every unordered multiset of five such masks and imposing pairwise
disjointness gives exactly **82 active-row configurations**.  This is a
lossless quotient by row permutations.  Pair- and triple-link intervals are
deliberately not used here, so the enumeration is a proof-safe superset.

For each active configuration, the nine isolate columns are enumerated as a
complete `S9` multiset with the exact row sums.  The five minimal covers of
the missing-edge `C5` give a cover prefilter.  A direct postcheck records zero
small-transversal failures, so that prefilter is exact.  The census is

```text
82 active-row configurations
1,098,455 graph-cover-prefilter paths
28,805 paths with repeated A5 rows
0 small-transversal failures after the prefilter
1,069,650 blocker-valid paths
35 paths with the scheduled 102 facets
15 exact unordered A5 families.
```

The exact target-family digest is
`27412b42ae9cbe2c09e915c70b2986f3d5ca4b90a11b758d0df674727534d9d5`.

### Withdrawal of the 32-configuration draft

An earlier exploratory draft enumerated only 32 active configurations by
allowing every A5 row to contain either zero or one active vertex.  That was
not complete: an independent A5 row can contain two nonadjacent vertices of
the active `C5`.  The 32-configuration draft and every conclusion based on
its purported completeness are withdrawn.  The theorem, certificate, hostile
replay, and tests use the corrected eleven-mask / 82-configuration census
above.

## Topology diagnostic and orbit

Every one of the fifteen scheduled target families has 102 minimal
transversals, all of order five.  Five families, with indices
`(0,1,2,3,6)`, have ridge-degree histogram `{2:459}` and pass the full
connected mod-two homology-sphere audit.  The other ten have histogram
`{1:8,2:443,3:8}` and already fail the pseudomanifold ridge condition.

The five topology-valid families form exactly one orbit under

```text
Aut(C5) x S9 x S5 = D5 x S9 x S5.
```

The verifier explicitly uses all ten `C5` automorphisms and all 120 row
permutations; modulo `S9`, equality is exactly equality of isolate-column
multisets.  The orbit-key digest is
`05c1d904e9aab3257c140300f94b7a3d3d13029736f70ab730db0b286cd0b9e2`.
A representative is

```text
(0,3,4,5,6)
(1,3,4,5,6)
(2,3,4,5,6)
(2,7,8,9,10)
(7,8,9,10,12).
```

The compact theorem does not need to discard the ten ridge-invalid families
first: the unrestricted label gate below kills all fifteen exact scheduled
targets.

## Complete compact Coxeter-label obstruction

Each minimal A5 must carry one of the five classical order-five Lanner
types.  The program regenerates all 420 labelled placements.  The finite
meeting-label alphabet is the complete `{2,3,4,5}` alphabet for compact
`H9/14`, by Burcroff's Lemma 10.10.

Supports are joined sequentially.  At every stage, every actual rank-nine
vertex restriction contained in the current union must be finite.  Pairs not
yet assigned by a Lanner support are provisionally right-angled.  Rejecting a
nonelliptic provisional restriction is safe: adding a non-right edge or
increasing its Coxeter weight cannot restore ellipticity, by
Perron--Frobenius monotonicity for Coxeter Gram matrices.

For representative family zero, the exact stage counts are

| stage | raw joins | compatible joins | passing states |
|---:|---:|---:|---:|
| 0 | 420 | 420 | 420 |
| 1 | 176,400 | 1,332 | 96 |
| 2 | 40,320 | 384 | 0 |

The raw-join column and the `compatible -> passing` columns are distinct
counts; they must not be read as one linear survivor chain.  The same complete
search ends at zero for each of the fifteen families.  The production label
transcript digest is
`5c9bab8a166c58bdbaa1e3179d900a42059ba53ff2f3aa15dbdb349e882d3554`;
the independent standard-library implementation has transcript digest
`0e726439a3f85977f4456b3184e74106050804c0f2e106456aad115df47b298e`.

Combining the schedule exclusions, the empty `a=2` branch, the single
excluded product orbit at `a=1`, and the complete fifteen-family label no-go
at `a=0` proves the whole-shadow compact theorem.

## Certificates and reproduction

The theorem is implemented by

- `src/coxeter_polytopes/combinatorics/h9_14_record42_all_completions.py`;
- `scripts/certify_h9_14_record42_all_completions.py`;
- `scripts/independent_h9_14_record42_all_completions_audit.py`;
- `data/enumerations/h9_14_record42_all_completions.json`; and
- `tests/test_h9_14_record42_all_completions.py`.

Reproduce it with

```sh
PYTHONPATH=src:. python3 scripts/certify_h9_14_record42_all_completions.py --verify-only
PYTHONPATH=src:. python3 scripts/certify_h9_14_record42_all_completions.py --verify-only --recompute
python3 scripts/independent_h9_14_record42_all_completions_audit.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_record42_all_completions.py
```

The hostile replay imports no project module.  It independently reconstructs
the schedules, `a=1,2` incidence gates, corrected 82-pattern `a=0` census,
exact topology orbit, all 420 labelled Lanner placements, and all fifteen
sequential label searches.  Its selected-product chain requires the
re-audited aggregate schema, the explicit support-symmetry orders `[24,120]`,
the unrestricted Record-42 flag, and the fixed-principal-block lemma.

The frozen hostile-summary digest is
`55b0206407df3aa6c3d9c14a76070db8a03f165a06798b670b148eb9ce5e30e6`.
The certificate digest is
`4590a0f0365dd26abd416cf4f73b909028a2863913f7c6d179a2a3e248122e9d`;
the archive file digest is
`c9beee7a93965774bff273b89c95430a2b9504b9d094b4730f78033aa659e74c`.
