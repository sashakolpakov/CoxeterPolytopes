# Exact partition of the frozen distance-nine replay

Source-only successor:
`scripts/check_selective_distance9_batches_20261020.py`, SHA256
`d9ab46e4404ed8dc02b0166f220cb0e63bc62eb42beaf68a5a7823fa65f84745`.
This file and contract are **UNEXECUTED**; no import, syntax-check job,
scientific replay, or new timing claim has been made. Root owns all runs.

The immutable original source is
`scripts/fresh_selective_distance9_dimension108_20261017.py`, SHA256
`f5cd7d0181cf9d5f96b10b6d4d6e7654caa3369b18ffa7473a16033f3124c81a`.
Its unchanged 22,358,471-byte certificate is
`results/fresh_selective_distance9_dimension108_20261017.json`, SHA256
`e7c6d1a1e56d8788e7d7a7d02bd5dd71837e3a3aef4ae768c1efc0a187bcacd9`.
The original monolithic replay timed out at the parent cap. This successor
changes its execution partition, not its mathematical domain or evidence.

The initial source-only revision was
`66dd350e1271e6f24aec231abe7917c6b4dbd33ce22c37b48472ff03abb12427`.
Root preserved it and made two small driver corrections: an explicit
`if not __debug__: raise RuntimeError(...)` guard remains effective under
optimization, and the root JSON endpoint is found structurally before
requiring only whitespace through EOF. The latter accepts root arrays with
trailing whitespace and rejects trailing nonwhitespace. The updated pinned
projection helper is
`be8a42d1d3687204af4d559c6bd8a28779dbfac68aeccdd95d41a7172645d555`.
These corrections do not change any mathematical candidate or output.

## Invocation and finite coverage contract

One invocation verifies exactly a half-open interval `[start,stop)` of one
phase. `stop` defaults to `start+1`; empty or out-of-range intervals fail.
The required complete integer domains are:

| Phase | Domain | Ownership |
|---|---:|---|
| `metadata` | `[0,1)` | Full schema, scalar constants, aggregate counters, exact survivor/case bijection and order, profile domain, boundary and compact-facet polynomials |
| `domains` | `[0,1)` | Complete 384 contexts/1262 ports, all4442 join keys, every1332 multiple-edge key, sorted85 inertia-type keys |
| `overlap-domain` | `[0,384)` | Every raw partial matching and its injectivity/contact-consistency pruning, by entire left context |
| `disjoint` | `[0,384)` | Join obstructions and **all** multi-edge subsets, by entire left context |
| `inertia` | `[0,85)` | Exact two-negative-direction congruence, by recomputed sorted type key |
| `overlaps` | `[0,384)` | Inverse/actual-edge rejection or survival for every consistent matching, by entire left context |
| `isolation` | `[0,949)` | Full extension/corner proof, by contiguous survivor index |
| `profiles` | `[0,8)` | Complete reduced partition profiles at link dimensions108 through115 |

A complete replay means successful intervals covering every integer of
every domain. Intervals may be subdivided without changing any test.
Overlap between successful intervals is harmless but does not replace any
missing index. The printed receipt includes its phase, bounds, complete
domain map, own source hash and immutable certificate hash. No new scientific
certificate is written. The wrapper may retain its ordinary compact ledger.

For example, root may run one interval as follows:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py \
  --parent-lock --ledger results/distance9_disjoint_000_008_replay_20261020.json -- \
  scripts/check_selective_distance9_batches_20261020.py --verify \
  --phase disjoint --start 0 --stop 8
```

This is an example partition, not a measured safe chunk size. Start with
small intervals and choose larger ones only from actual capped timings.
`metadata` and `domains` are bounded finite-domain/aggregation checks, not
the full obstruction computations. All other expensive domains can be
reduced to single ownership units if necessary. A timeout or memory failure
is a failed/incomplete interval; it never suppresses candidates or a phase.

Every invocation requires the existing parent-held exclusive proof lock,
checks its live wrapper parent, and has a45-second child alarm. The unchanged
wrapper enforces59 seconds,1GiB process-tree RSS, one numerical thread, and
the user resource gates. The script is replay-only and rejects disabled
assertions. No detached workers or parallel proof processes are introduced.

## Why the candidate partitions are exhaustive

The frozen `contexts_and_ports` function is reused unchanged. It derives
every pair of finite black deletions from every original census diagram,
recomputes distances and inverse-diagonal tests, and checks384 contexts,
1262 ports and378 eligible ports. No saved context list defines this domain.

The join-domain routine takes the **entire ordered Cartesian product** of
eligible ports and applies precisely the original distance and white Schur
tests. It reconstructs4442 candidates,1103 white rejections, all3339 other
joins, the complete context-pair multiplicity distribution and85 Gram keys.
All multi-edge subsets are enumerated inside the whole key
`(left diagram,left black pair,right diagram,right black pair)`.
The left context owns every such subset. In particular, no left-port split
can lose a subset using two different ports. `domains` compares the whole
ordered candidate-key and subset-key arrays; the obstruction phases compare
every complete output record in their corresponding domains.

The overlap generator is a localized transcription of the frozen outer
loops: only the left-context interval is restricted. Every right context,
eligible node pair, matching size, injectivity test and shared-contact test
is unchanged. It does not start from saved overlaps. An independent
combinations count gives the global312316 unpruned matchings. The domain
phases compare every consistent matching key against the saved20718 records.
`metadata` checks all saved owners are actual contexts and that global
left-context order is nondecreasing. Hence the interval comparisons prove
the complete global list and its original order, not just membership of
selected saved rows. The obstruction phase then recomputes each original
inverse/PD decision and every witness field.

Isolation consumes the949 surviving records established by those complete
overlap phases. The metadata phase verifies an exact index-by-index
bijection between these records and the isolation cases, including left and
right marks, matching, type and core key. The54 core-key set is checked
exactly. Each isolation interval rebuilds its native unions, all extra-cross
tests, rationalization, kernel or inverse, exact inertia, every preferred
corner candidate and the frozen tie-breaking rule. It compares **all** case
fields and each encountered whole core proof. Core inverses are cached only
within an interval, keyed by the whole native core; black-context variation
cannot change that matrix. This replaces repeated identical inversions,
not a mathematical test or candidate. Every core occurs in at least one
of the949 bijected cases.

## Equivalence to all original assertions and outputs

The only localized loop transcriptions are the original `disjoint`,
`overlaps`, and `isolation` loops. They retain the same loop order, bond
comparison, positive-vector witness, preferred-corner choice and output
construction. Matrix calculations reuse the frozen functions, including
`rationalize`, `z_obstruction`, `exact_inertia`, `edge_obstruction`,
`extra_cross`, `corner`, `ordinary_corner` and `relevant_corner`.
Single-join inertia is assigned its own complete85-key phase. The main
frozen `run()` is never invoked.

The original scalar functions `profiles` and `budget` are reused. The
metadata phase independently reconstructs the complete756 profile keys
and their order; each profile phase recomputes all coefficients, moments,
curvature, costs, strict margins, finite-vertex slack and dimension minimum.
The negative boundary profile and both compact-facet polynomials are
recomputed too. Aggregate checks retain20718 consistent overlaps,18721
inverse rejections,1048 PD rejections,949 survivors,54 cores,60831 shortcuts,
888 extra-bond PD rejections,69 relevant-component cases, order cap17,
756 scalar rows and minimum margin45873/2476099. Global counters are tied
to the exact case records, which are themselves recomputed by the required
intervals; they are not accepted solely from an outer status string.

Every top-level and nested output family is covered: constants/schema,
sources, ports, joins, inertia, multi-edge obstructions, overlap counts and
records, isolation cores/cases/scalars, reduced profiles, dimension minima,
boundary and compact-facet polynomials. JSON normalization matches the
frozen run (tuples become arrays; integer object keys become strings).
Arrays retain their original order; dictionary key order is immaterial,
just as in the original `json.loads` equality.

## Input fidelity, freshness and allocation

Every invocation checks the original full23-path outer source domain and
all seven inherited input statuses and `sources` maps **before** importing
the frozen mathematical definitions. Shared input paths are hashed once
per invocation and incompatible inherited hashes fail. Every collected
input is streamed and rehashed before success. The worker, original source,
certificate, existing wrapper and reusable projection helper are also
checked before/after. Thus separating original inherited checks into an
earlier successful invocation cannot leave later computations unbound.

The existing `JsonProjection` is reused for small object-path projections.
The new array reader traverses the read-only mmap and decodes one item at
a time, with a1MiB item cap and a32MiB total-certificate cap. The full21MiB
certificate is never decoded as one object. A selected ownership interval
is retained only as needed for exact equality; each scan traverses the
complete array and returns its length, and the domain/metadata phases check
all global lengths. All allocation bounds fail closed. Source/script
imports are definition-only; guarded mathematical producers are not run
on import. The pure-stdlib projection/driver retains the original exact
SymPy arithmetic through those frozen functions.

No frozen source, certificate, manuscript, catalogue or status file was
changed. Runtime/RSS feasibility of particular interval widths remains for
root's guarded executions; no successful replay is claimed by this source
preparation.
