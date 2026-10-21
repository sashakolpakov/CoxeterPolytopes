# Record 47: whole-shadow compact-Coxeter no-go

## Theorem and scope

Under the standard compact `H9/14` missing-face reductions, **record 47 has
exactly one residual incidence orbit, and that orbit is not a compact
hyperbolic Coxeter nerve**.  Thus the entire record-47 missing-edge graph
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
vertex-, edge-, and triple-link interval constraints leave 322.  For each of
these the verifier constructs every forced five-subset containing the A4 and
applies strict vertex-link positivity before adding the remaining minimal
A5s, together with the aggregate Macaulay incidence bound.  The complete
histogram is

| forced five-shadow | additional A5s | A4s | strict feasible |
|---:|---:|---:|---:|
| 8 | 2 | 84 | 0 |
| 9 | 1 | 168 | 0 |
| 10 | 0 | 70 | 0 |

Hence the `a=1` branch is empty.

## Exhaustive `a=0` census

The six active vertices are

```text
0, 1, 2, 11, 12, 13
```

and vertices `3,...,10` are eight interchangeable graph isolates.  There are
four inclusion-minimal graph covers, all of order three:

```text
(0,2,12), (0,1,13), (0,12,13), (11,12,13).
```

There must be exactly four missing A5s.  The census is lossless:

1. an unordered A5 family is first split into its four active-vertex rows;
2. the strict vertex-link caps leave eight possible active row masks and 78
   active-intersection multisets;
3. for fixed rows, an assignment of the eight isolates is a multiset of
   four-bit incidence columns, so its `S8` orbit is determined exactly by the
   sixteen column multiplicities;
4. every resulting column-distribution problem is solved by exact integer
   recursion; and
5. any family with a transversal of order at most four is discarded, while
   every surviving order-five graph cover meeting all four A5s is counted.

The initial recursion does not claim to be the final quotient: equal active
rows can create duplicate enumeration paths.  The output is therefore
deduplicated first as exact four-mask clutters and then canonically over all
24 row permutations and the `S8` isolate-column action.  No active-graph
automorphism quotient is used, so the final count can only retain additional
representatives, never identify distinct active configurations accidentally.

The census rejects 18 column assignments because at least two A5 rows
coincide, and retains 744 enumeration paths with four distinct rows.  These
give 618 exact four-mask clutters and 247 canonical
`S4 x S8` states.  None of the 744 paths has a missed small transversal.  The
required facet count 90 occurs exactly once at the path, exact-mask, and
canonical-quotient levels.  Its four A5s are

```text
(3,4,5,6,7)
(1,2,8,9,10)
(1,8,9,10,11)
(2,8,9,10,11).
```

They are literally the frozen `TARGET_MISSING_FIVES` of the previously
selected record-47 nerve.  Direct blocker dualization gives exactly ninety
minimal transversals, all of order five.  The resulting face table has all
405 ridges of degree two, connected dual graph, the scheduled face vector,
and spherical Euler and mod-two homology for every face link.  Its frozen
facet-complement and packed-face-table digests are

```text
8306d9c6c39a90aae519bd5ae587743b3737371f9efbe01ebd3757105b0b56ed
4d07f78cba8551c9d0fe9fc4d6041f3df8ae35c5017e9942e10005686c628aba
```

The topology audit is not being used to infer uniqueness: uniqueness already
holds in a deliberately coarse necessary-condition census.  The full audit
verifies that the sole survivor is the previously frozen globally audited
incidence completion.  Ordinary polytopality is supplied separately by the
multiwedge theorem, not by the selection archive alone.

## Chained unrestricted Coxeter obstruction

The re-audited seven-selected-nerve theorem is unrestricted in the finite
Coxeter labels for this exact nerve.  For record 47 it verifies that each
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
what upgrades that fixed-nerve conclusion to the whole record-47 graph
shadow.

## Certificates and reproduction

The theorem is implemented by

- `src/coxeter_polytopes/combinatorics/h9_14_record47_all_completions.py`;
- `scripts/certify_h9_14_record47_all_completions.py`;
- `scripts/independent_h9_14_record47_all_completions_audit.py`;
- `data/enumerations/h9_14_record47_all_completions.json`; and
- `tests/test_h9_14_record47_all_completions.py`.

Reproduce it with

```sh
PYTHONPATH=src:. python3 scripts/certify_h9_14_record47_all_completions.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record47_all_completions_audit.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_record47_all_completions.py
```

The independent replay uses only the Python standard library and does not
import the production combinatorics package.  It reconstructs the schedule,
both surviving branch gates, the unique blocker and face table, and the
record-47 projection of the frozen unrestricted label certificate.

The frozen certificate digest is
`f9a3dd87824110249815aca2abaae8348aba0bfd8c439732aafe141494a7e700`;
the independent hostile-summary digest is
`67dff778d408b831c2dd499a634caa7c3af816bcdecfe59e01ffbd9dc44ac6fa`.
