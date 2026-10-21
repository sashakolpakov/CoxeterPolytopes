# H9/14 Record 30: whole-shadow compact Coxeter no-go

## Theorem

There is no compact hyperbolic Coxeter polytope in dimension nine with
fourteen facets and Record-30 missing-edge graph

```text
3K2 + 8K1, with missing edges (0,11), (1,12), (2,13).
```

This is a whole-graph-shadow theorem: it covers every compact
Dehn--Sommerville schedule and every residual necessary-incidence completion,
not merely the displayed ordinary product nerve.

It is also deliberately **compact-only**.  No ideal-vertex, affine-cusp,
maximal-parabolic, or Vinberg finite-volume endpoint audit is part of the
proof.  Consequently the theorem does not exclude a finite-volume noncompact
Record-30 chamber.  It constructs no Lorentzian Gram realization and gives no
new hyperbolic Coxeter example.

The frozen aggregate is
[`h9_14_record30_whole_shadow_compact_no_go.json`](../data/enumerations/h9_14_record30_whole_shadow_compact_no_go.json),
verified by
[`certify_h9_14_record30_whole_shadow_compact_no_go.py`](../scripts/certify_h9_14_record30_whole_shadow_compact_no_go.py).

## Why the eight schedules are exhaustive

Let `a` be the number of missing independent four-sets.  For the fixed
missing-edge graph, Dehn--Sommerville gives

```text
b = 8 + 6a,       number of facets = 132 - 2a,
global g4 = 8-a,  sum of the fourteen vertex-link g4 values = 6(8-a).
```

Every minimal nonface of a compact Coxeter nerve is a Lanner diagram of order
at most five.  The Murai--Nevo equality theorem therefore makes the relevant
global and link `g` entries strict; the full argument is frozen in
[`h9_14_local_link_constraints.md`](h9_14_local_link_constraints.md).  Since
global `g4>0`, every possible nonnegative integer `a` satisfies `a<=7`.
Moreover every one of the fourteen vertex-link entries is at least one, so
their sum is at least fourteen; this already excludes `a=6,7`.

The complete numerical table is

| `a` | `b` | facets | global `g4` | sum link `g4` | disposition |
|---:|---:|---:|---:|---:|---|
| 0 | 8 | 132 | 8 | 48 | dual CNF/DRAT exhaustion plus label no-go |
| 1 | 14 | 130 | 7 | 42 | strict edge-link branch plus two exact incidence replays and label no-go |
| 2 | 20 | 128 | 6 | 36 | 77-orbit incidence exhaustion plus complete label/Gram no-go |
| 3 | 26 | 126 | 5 | 30 | exact facet-count incidence contradiction |
| 4 | 32 | 124 | 4 | 24 | exact facet-count incidence contradiction |
| 5 | 38 | 122 | 3 | 18 | unique positive orbit has a forced edge-link with top `g=0` |
| 6 | 44 | 120 | 2 | 12 | strict summed vertex-link bound |
| 7 | 50 | 118 | 1 | 6 | strict summed vertex-link bound |
| `>=8` | — | — | `<=0` | — | strict global `g4` bound |

Thus the branch certificates close `a=0,...,7`, and the analytic strict bound
closes every larger nonnegative value.  There is no omitted schedule.

## Exact branch certificates

### `a=0`

The production connected-ridge census has eight incidence orbits, and the
complete finite-label consumer leaves zero states on each.  Completeness is
certified by a structurally different dual Boolean encoding.  It chooses the
eight missing independent five-sets directly and imposes the blocker
biconditional, exact 132-facet count, local caps, absence of smaller
transversals, and every ridge degree `0` or `2`.  These are only necessary
conditions, so this CNF describes a superset of the desired targets.

After fixing a lossless type-three seed and blocking the complete 63,960
seed forms of the eight known target orbits, the formula has 226,358 variables
and 984,935 clauses.  Kissat reports `UNSATISFIABLE`; the archived DRAT proof
is independently `VERIFIED` by `drat-trim` in 34,109,825 resolution steps.
The checker was rerun independently in 488.306 seconds with the exact pinned
binary.  The hardened deep verifier rejects a solver or checker whose binary
SHA-256 differs from the archived provenance.

A separate no-project-import hostile replay reconstructs
`Aut(3K2) x S8`, all eight full family orbits, their stabilizers, all 63,960
seed forms, and the exact set of final CNF blocking clauses.  The combined
hostile audit reports `PASS` with no remaining fatal or material gap.  See
[`the a=0 theorem`](h9_14_record30_a0_dual_exhaustion.md) and its
[`final hostile audit`](h9_14_record30_a0_dual_hostile_audit.md).

### `a=1`

There are three fixed-`A4` orbits.  The two-active orbit violates strict
edge-link `g3>0`.  On the other two branches, a nauty colored-incidence
quotient and a no-nauty row/column canonical traversal have identical state
transcripts and literally equal target sets: three target orbits of sizes
`48,80,40`, or 168 labelled families.  Complete order-four/order-five Lanner
joins on all 130 actual rank-nine vertices leave zero finite-label states.

### `a=2`

Two independent quotients exhaust all 21 positive-`A4` orbits and obtain
exactly 77 pure 128-facet target orbits.  Complete Lanner/actual-vertex joins
kill 74.  The selected product target has its frozen principal-block no-go.
For each of the other two targets, every one of 1,332 symmetry-quotient
cross-shell assignments has an exactly positive, nonzero determinant on an
eleven-wall dotted-free principal Gram block.  Rank at least eleven
contradicts Gram rank at most ten in `H9`.

### `a=3,4`

For `a=3`, independent meet-in-the-middle and pivoted exact-cover searches
agree across all 32 positive-`A4` orbits and 1,896,496 no-small-transversal
families.  None has the scheduled 126 order-five blockers.

For `a=4`, six of seven positive orbits retain a small transversal.  The last
has 840 admissible triples, whose order-five blocker counts are only
`92,97,102,107`, never the required 124.  A literal 16,328,620-triple search
and an independent pivoted replay agree.

### `a=5,6,7`

The `a=5` positive quotient has one full orbit after representing 672 labelled
families.  Its forced edge link has

```text
h=(1,3,4,4),  g=(2,1,0),
```

and no later missing five-set can make the final entry positive.  For `a=6,7`,
the summed vertex-link `g4` values are respectively 12 and 6, both below the
strict minimum fourteen.

## Frozen composition and checksum roles

The aggregate verifier first runs every component's semantic verifier.  Only
after those checks succeed does it compare the following archive-file and
internal-certificate digests.  They are reproducibility checksums, not a
substitute for the mathematical and combinatorial verification above.

| schedules | component archive SHA-256 | internal certificate SHA-256 |
|---|---|---|
| `a=0` | `7934ffec51f8bc3c5f209d016f75c0aa6feac34343ad34a58f16a0d4660a64e5` | `436b9bc834529f7d2d94754e43c3f785bb937ec0f7157570ad9a06d56858ab2d` |
| `a=1` | `65152e8f4b687f487852932a5dae3896d189a56d1ed99c6764f7389aaea8d341` | `932d69e5646106785171562174e5423329ac09137e0224ac65e6ce436187067d` |
| `a=2` | `06759a184708d9c2b8a4ab15fc772fc2a4b7e66ac275d2ff537bebf92cb9ebe6` | `3ac5e4ea4ddd624aef17dec45f64d403165b736bab015d088e2b8d78f43a5adf` |
| `a=3` | `2e959be2b47424ebe6b9780848202fd800055dd8f093ea3cf675336947abcb5c` | `e3afaf33526c003307d797bd0789bd2b4275d3d74c83f0a8d9792bb1ad0b9d54` |
| `a=4` | `4c23f1a87a35f51c9c277e08a18e7f89ca02eb51757c244d59270d005ecc6c46` | `12918a797758ec9b14ed02657c8a98060c77b54abaa93d0e627dc31d0efb3bb5` |
| `a=5,6,7` | `841ebf50748747a8dc96fb4ce337570f29f1f9455f2206609e18d31187eb3771` | `02a2238457b33d0cba59bb9fb85579a5e27c4520ea86f17a4ce0614fcabb111e` |

The aggregate certificate's internal SHA-256 is
`4215b56a4fe830825da8a8689395ae33f2fba105b83bbbc468bb6a7f054f170f`;
the archive-file SHA-256 is
`17a28b69ffc46e9712386971ed5e3c19edf6b952e0474272a6b00c51ce4a2b09`.
Both are checked by the certifier.

## Frontier consequence

Record 30 was ordinary-polytopal at the incidence level, so this theorem does
not change the incidence partition: it remains

```text
123 incidence-impossible + 19 ordinary-polytopal.
```

It does promote one formerly open compact graph shadow to a whole-shadow
no-go.  The compact ledger changes from `132 + 10` to

```text
133 wholly excluded + 9 open,
open = {0,1,3,4,5,7,9,14,31}.
```

## Reproduction

Run the frozen aggregate and focused tests with

```bash
PYTHONPATH=src:. python3 scripts/certify_h9_14_record30_whole_shadow_compact_no_go.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_record30_whole_shadow_compact_no_go.py
```

The expensive `a=0` proof can be replayed separately using the exact binaries
whose hashes are stored in the theorem object:

```bash
python3 scripts/certify_h9_14_record30_a0_dual_exhaustion.py \
  --deep --solver PATH_TO_KISSAT_4_0_4 \
  --drat-trim PATH_TO_PINNED_DRAT_TRIM
```

This aggregation completes one more primal/nonexistence step.  It does not
move the counterexample boundary and warrants no new-example `TOOT`.
