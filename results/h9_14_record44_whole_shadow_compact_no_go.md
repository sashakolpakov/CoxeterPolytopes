# Record 44: whole-shadow compact-Coxeter no-go

## Theorem and scope

Under the standard compact `H9/14` missing-face reductions, **record 44 has
exactly one residual incidence orbit, and that orbit is not a compact
hyperbolic Coxeter nerve**.  Thus the entire record-44 missing-edge graph
shadow is excluded on the compact track.

This is not a finite-volume noncompact theorem.  The proof uses strict
ellipticity at every actual vertex.  It does not enumerate ideal vertices or
maximal parabolic diagrams and does not apply Vinberg's finite-volume
criterion.  It is exhaustive under the compact-Coxeter missing-face bound;
it does not purport to classify arbitrary ordinary complexes with minimal
nonfaces of order greater than five.

## All four schedules

Write `a` for the number of missing four-subsets and `b` for the total number
of graph-independent five-subsets that are nonfaces.  When `a>0`, this
includes the forced upper shadow of every missing A4 and must not be confused
with the number of additional minimal A5s.  The complete
Dehn--Sommerville table is

| `a` | `b` | facets | sum of the fourteen vertex-link `g4` values |
|---:|---:|---:|---:|
| 0 | 4 | 90 | 24 |
| 1 | 10 | 88 | 18 |
| 2 | 16 | 86 | 12 |
| 3 | 22 | 84 | 6 |

Every vertex link of a compact Coxeter nerve satisfies strict `g4 >= 1`, so
the total must be at least fourteen.  This excludes `a=2,3` immediately and
leaves only `a=0,1`.

## Exact exclusion of `a=1`

The missing-edge graph has 690 independent four-subsets.  The exact
vertex-, edge-, and triple-link interval constraints leave 378.  For each of
these the verifier constructs every forced five-subset containing the A4 and
applies strict vertex-link positivity before adding the remaining minimal
A5s, together with the aggregate Macaulay incidence bound.  The complete
histogram is

| forced five-shadow | additional A5s | A4s | strict feasible |
|---:|---:|---:|---:|
| 8 | 2 | 252 | 0 |
| 10 | 0 | 126 | 0 |

Hence the `a=1` branch is empty.

## Exhaustive `a=0` census

The five active vertices are

```text
0, 1, 11, 12, 13
```

and vertices `2,...,10` are nine interchangeable graph isolates.  There are
four inclusion-minimal graph covers, all of order three:

```text
(0,1,11), (0,1,13), (0,12,13), (11,12,13).
```

There must be exactly four missing A5s.  The census is lossless:

1. an unordered A5 family is first split into its four active-vertex rows;
2. the strict vertex-link caps leave six possible active row masks and 12
   active-intersection multisets;
3. for fixed rows, an assignment of the nine isolates is a multiset of
   four-bit incidence columns, so its `S9` orbit is determined exactly by the
   sixteen column multiplicities;
4. every resulting column-distribution problem is solved by exact integer
   recursion; and
5. any family with a transversal of order at most four is discarded, while
   every surviving order-five graph cover meeting all four A5s is counted.

The initial recursion does not claim to be the final quotient: equal active
rows can create duplicate enumeration paths.  The output is therefore
deduplicated first as exact four-mask clutters and then canonically over all
24 row permutations and the `S9` isolate-column action.  No active-graph
automorphism quotient is used, so 502 is not claimed to be the number of
full graph-isomorphism orbits.  This deliberately coarser quotient can only
retain additional representatives; it cannot merge distinct active
configurations accidentally.  Since the 90-facet target is already unique
as a labelled enumeration path, no further quotient is needed for the theorem.

The census rejects 79 column assignments because at least two A5 rows
coincide, and retains 2,782 enumeration paths with four distinct rows.  These
give 2,214 exact four-mask clutters and 502 canonical
`S4 x S9` subgroup-quotient states.  None of the 2,782 paths has a missed
small transversal.  The required facet count 90 occurs exactly once at the
path, exact-mask, and canonical-quotient levels.  Its four A5s are

```text
(2,3,4,5,6)
(1,7,8,9,10)
(7,8,9,10,11)
(7,8,9,10,12).
```

They are literally the frozen `TARGET_MISSING_FIVES` of the previously
selected record-44 nerve.  Direct blocker dualization gives exactly ninety
minimal transversals, all of order five.  The resulting face table has all
405 ridges of degree two, connected dual graph, the scheduled face vector,
and spherical Euler and mod-two homology for every face link.  Its frozen
facet-complement and packed-face-table digests are

```text
a706233135e61171db1da19172f4cb5234191f0a0a703a3388b4a4d9a9877917
449e514bd9f090a3c52dba38db833014e696c7d9cb20753a6b8f4f7390cb6ad6
```

The topology audit is not being used to infer uniqueness: uniqueness already
holds in a deliberately coarse necessary-condition census.  The full audit
verifies that the sole survivor is the previously frozen globally audited
incidence completion.  Ordinary polytopality is supplied separately by the
multiwedge theorem, not by the selection archive alone.

## Chained unrestricted Coxeter obstruction

The re-audited seven-selected-nerve theorem is unrestricted in the finite
Coxeter labels for this exact nerve.  For record 44 it verifies that each
fixed left support has the full order `120` normalizing action and uses three
complete support gates, each with 52,500 one-edge trials.  Labels 4 and 5 have no
survivors; Perron--Frobenius monotonicity covers all larger labels.  All
1,890,000 raw two-edge trials are covered, and the only exceptional block
types have exact inertia `(8,2,0)`, contradicting the single negative
direction of a hyperbolic Gram matrix.  Entries involving nodes outside the
fixed ten-node principal block do not change that block or its inertia, so
the obstruction is genuinely unrestricted by Cauchy interlacing.  Therefore
this unique incidence orbit has no compact Coxeter labeling.

The label certificate correctly describes itself as
"selected-nerve-only."  The new exhaustive uniqueness theorem is precisely
what upgrades that fixed-nerve conclusion to the whole record-44 graph
shadow.

## Certificates and reproduction

The theorem is implemented by

- `src/coxeter_polytopes/combinatorics/h9_14_record44_all_completions.py`;
- `scripts/certify_h9_14_record44_all_completions.py`;
- `scripts/independent_h9_14_record44_all_completions_audit.py`;
- `data/enumerations/h9_14_record44_all_completions.json`; and
- `tests/test_h9_14_record44_all_completions.py`.

Reproduce it with

```sh
PYTHONPATH=src:. python3 scripts/certify_h9_14_record44_all_completions.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record44_all_completions_audit.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_record44_all_completions.py
```

The independent replay uses only the Python standard library and does not
import the production combinatorics package.  It reconstructs the schedule,
both surviving branch gates, the unique blocker and face table, and the
record-44 projection of the frozen unrestricted label certificate.

The frozen certificate digest is
`bf58d2c32c5791c61fadf9b8b9b386f98e2e1e65a429bfd9bbfb32d134a97e90`;
the independent hostile-summary digest is
`caf58222d6ee59e8a658a441b0991641118a8f6e87b2750b065c4b46cc972ec5`.
