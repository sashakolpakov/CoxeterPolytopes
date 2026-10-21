# Distance-nine replay: source-directed batching advice

Source-only analysis of
`scripts/check_selective_distance9_batches_20261020.py`,
SHA256 `126d5fe3588cfbb56009bff1f748dbab2154582e3df814feff50b32003cf023b`.
No worker, import, fixture or scientific job was run for this note; no source
was changed. Root owns every resource-gated serial pilot.

## What is a unit, and what is repeated

The complete domains have \(1+1+384+384+85+384+949+8=2196\) integer units.
They are **coverage units, not 2196 required processes**. Any contiguous
half-open interval is accepted. A receipt for a larger interval covers all
its units; overlapping receipts are harmless but do not cover missing units.

Every invocation repeats input pin checks, imports only the frozen
definitions, and rehashes all collected inputs before success. Every phase
except `profiles` also recomputes the full census: 384
contexts, 1262 ports, 378 eligible ports and their deletion inverses.
No phase interval restricts this setup.

| Phase | Current indivisible coverage unit | Slice-independent and variable work |
|---|---|---|
| metadata | the whole [0,1) | Entire schema/counter/bijection/profile-domain checks; several full certificate traversals. Cannot be made smaller by changing start/stop. |
| domains | the whole [0,1) | All \(378^2=142884\) ordered port pairs, distance-pruned to4442 joins; all1332 multi-edge keys and85 type keys. No inertia witnesses are recomputed here. |
| overlap-domain | one entire left context | All right contexts and every partial matching size for that left context; all injectivity/contact-consistency decisions. The whole phase has312316 raw matchings and20718 consistent ones. |
| disjoint | one entire left context | **Every slice first recomputes the complete join domain**, even a context owning no joins. It then recomputes its single-join and whole context-pair multi-edge witnesses. |
| inertia | one sorted underlying Gram type | **Every slice first recomputes the complete join domain.** Each selected type then gets its exact congruence, at matrix order at most20. |
| overlaps | one entire left context | The same raw matching generator as overlap-domain, then all inverse/PD/survival decisions for its consistent matchings. |
| isolation | one complete survivor case | Every slice first scans all20718 overlap rows to select949 survivors and scans all949 saved cases. Selected cases rebuild all extra-cross tests and all four marked corner candidates. Core cache is local to this invocation. |
| profiles | one whole link dimension | Skips census setup. It still performs all pin/freshness work, then recomputes every scalar profile at the selected dimension(s). |

`Arrays.select` decodes every row of its selected array before
filtering; it is memory bounded but not random-access by owner. Its path
locator and small `saved(...)` projections also traverse skipped
frozen subtrees. Therefore a one-unit disjoint/overlap/isolation invocation
does **not** have one-unit input-reading cost. New processes also lose the
deletion/inverse caches. Repeating hundreds of tiny invocations pays these
fixed costs repeatedly.

## Ownership and the largest atomic work

Single joins and multiple-edge subsets must stay under the entire
`(left diagram,left black pair)` context. The right context
and **all its ports** remain unrestricted. In particular, splitting by a
left port would lose subsets spanning left ports. The `good`
groups have1,2,3 or4 cross edges; every subset of size2 through the full group
is checked. Across all groups there are1332 such subsets.

Overlap units likewise retain all right contexts and all partial matchings
for one left context. A single match can identify several ports; selecting
one matching edge is not a valid ownership partition. The largest match
size is five. Native unions have order at most19, and the exact rejected
edge blocks are smaller. Most of the20718 rows stop at the inverse test
(18721); 1048 stop at the actual-edge PD test and949 survive.

An isolation unit is substantially richer than a stored preferred-corner
record: it recomputes extra-cross tests and *every* marked corner candidate,
then the original tie-breaking rule. With two black marks on each side
there are four candidate corners. At most nine unshared nodes remain per
side, hence at most81 candidate extra-cross edges before existing-edge
skips. The case may need rationalization, nullspace/inverse, exact inertia,
corner block decompositions, and residual/relevant-component tests. It
cannot be split further by current start/stop options. An actual single-case
timeout would therefore require a separately reviewed finer partition,
not skipping that case or raising a cap.

There are54 native core keys for949 isolation cases. Larger chunks can
reuse a key across different black contexts because the whole Gram matrix
depends only on the two diagrams and the identification matching. This
cache benefit disappears in a one-case process.

For profiles the complete counts can be derived directly from the printed
generator: at link dimension \(d=108,\ldots,115\), there are \(d-17\)
profiles, namely91 through98, totaling756. Polynomial degree is at most115,
with at most seven factors, and the independent moment uses at most128
factor subsets. This is the strongest candidate for one full [0,8) batch;
it is not eight required processes.

## Economical measured pilot strategy

The existing predecessor metadata pilot took14.295s and133562368 bytes RSS.
The one-line Boolean count fix does not touch that phase, but its old-source
receipt does **not** establish current-source coverage. The old domains
failure at1.954s stopped before the work completed and is not a domains
runtime estimate. The current fixture's resource-gate deferral is not a
successful fixture run. At the time this note was prepared, no further
current-source phase timing was asserted here.

After root's current-source fixture/coverage self-tests and launch gate:

1. Run metadata [0,1) and domains [0,1). These are indivisible; their measured
   fixed costs are useful baselines. Run profiles [0,8) as one capped pilot.
2. Probe the larger-order contexts, rather than only an empty/light early
   context. Census diagram25 is the final diagram; its nine finite deletion
   nodes produce the final36 contexts, [348,384). Diagram24 contributes
   [327,348). Thus [348,349) is a sensible first hot overlap pilot; [348,352)
   is a useful small context block after that passes. This identifies matrix
   order and long distances, not a proven maximal runtime.
3. For overlap-domain, try to coalesce aggressively once the hot block is
   measured: its whole raw domain is only312316 combinations. Use64/128 or
   the full384 contexts if actual timings justify that; do not default to384
   processes. The overlap *classification* phase has heavier per-row work,
   so it needs its own pilot rather than borrowing the domain-only timing.
4. For disjoint, measure a hot whole-context block before using32/64-context
   chunks. The full join reconstruction is paid regardless of chunk width,
   so even light contexts should be bundled. Use the saved row counts only
   to choose a workload-balanced schedule, never to prune the recomputed
   candidate domain.
5. For inertia, probe the final5 sorted types [80,85), which include the
   largest-order end of the keys, then choose16/32 or all85 types from the
   result. A one-type plan repeats the full join reconstruction85 times.
6. For isolation, start with8 contiguous cases, plus a pilot containing a
   saved relevant-component case if the first block contains none. Then grow
   to16/32/64 cases according to actual wall/RSS results. A late block such
   as [941,949) samples the final large-diagram cases. Keep adjacent cases
   together where possible to exploit repeated core keys.

These widths are pilot suggestions, **not measured safe sizes**. Keep normal
successful chunks comfortably below the45s child alarm, for example around
20--30s including pinning/scanning, rather than budgeting to the59s parent
limit. Doubling a short measured chunk is preferable to assuming linear
cost from a single case: owner loads and core-cache reuse are nonuniform.
Compare two widths to distinguish fixed scanning/setup from variable work.
If a proposed interval times out, retain all its indices as incomplete and
bisect it; do not issue coverage credit for partial progress.

Receipts already covering a hot pilot need not be repeated merely to get
regular interval widths; cover the complementary intervals or accept the
small duplication. The final coverage check must still find every integer
in all eight domains under the **current** worker pin. Neither selecting a
pilot from saved data nor completing metadata proves the missing units.

## When source work, rather than more tiny batches, would be warranted

No source change is proposed now. If measured fixed overhead dominates,
the current implementation already permits larger intervals and should use
them first. If a whole left-context overlap unit itself exceeds45s, a future
reviewed partition could refine ownership by *whole left/right context pair*,
preserving all partial matchings inside that pair. A disjoint refinement
would likewise keep every multi-edge group whole. A single isolation case
would need a separate output-field/inner-domain contract. None of those
finer partitions exists in the current worker, and none is needed merely
because2196 coverage units are required.

## Engineering appendix: coordinator source review

This is a metadata/source review, not mathematical replay coverage. The
coordinator's author re-read its complete receipt and cleanup paths against
the separate coverage and bounded-runner contracts; this is not an
outside-author audit. The held coordinator is
`scripts/resume_distance9_replay_20261020.py`, SHA256
`1e3432ec083ce1fd6fce65f8b9ad50ca202e8680935eaa13c8cf6b97ebc4f6d1`.
Its byte-verified predecessor is archived as
`manuscripts/source_archives/resume_distance9_replay_before_frozen_guard_20261020.py.gz`,
decoded SHA256
`372b569f15fc6d57486ef55139d1b0f3d0f54ce4cda9f7bf1457c15fef917fd8`.
The exact diff adds only the streaming frozen-path test, its execution/save
guards, and the flushed per-job outcome. No replay enumeration changed.

The current coverage contract is
`scripts/check_partitioned_replay_coverage_20261020.py`, SHA256
`41f3c82b86ad05aa0877a34cf4718f92d5253e72fe669c083dd161a22c4342c6`.
The coordinator imports only this and the stdlib-only definitions in
`verify_manuscript.py`. It checks their hashes and all of the coverage
checker's current fixed inputs. `analyze` credits only records whose
`parse_record` validates the current worker/source/certificate, exact argv
and printed interval, output hashes, and the wrapper's lock/resource
receipt. Stale, failed and malformed rows remain present but contribute
no bits. Top-level status, plans and historical nested objects add no credit.
For the first gap `[l,h)` in fixed phase order the next job is exactly
`[l,min(h,l+width))`; covered islands are not crossed. Each new successful
record is independently parsed before normal progress continues. Complete
coverage therefore means the union of accepted intervals contains all2196
units, not that2196 processes ran.

Execution refuses a ledger whose repository-relative path occurs in any
`results/*freeze.json`, before entering the job loop and before each save;
planning may read it. This is the existing builder's conservative literal
substring rule, not a new general freeze-manifest format. The separate
coordinator flock serializes coordinators. The old ledger digest is checked
before atomic replacement, preventing silent overwrite by a detected
noncooperating writer. Every saved job retains the actual `run_job` record,
including failure/defer rows, and then prints/flushed its phase, interval,
actual status and elapsed time. Final output separately reports coverage.

For handled SIGINT/TERM/HUP and parent-loss polling, the first interrupt
disarms the poll and suppresses repeated termination signals before raising
`KeyboardInterrupt`. The active `run_job` then executes its existing
`finally`: terminate/wait/kill the wrapper; the wrapper's own `finally`
stops its checker process group and releases the shared proof lock. Only
after this unwind does the coordinator save an interrupted metadata event;
it manufactures no successful receipt for an interrupted job. Final ledger
writes disable the poll and repeated termination signals. There is no
new child lifecycle or background proof process.

No source-level blocker was found for ordinary coordinated use. Explicit
limits remain: this is an unsigned-receipt protocol; the coordinator reads
one supplied jobs ledger rather than discovering other receipt files; all
writers must respect ownership of that ledger; freezes/external writes,
I/O errors and hard SIGKILL/power loss can prevent recording the latest
attempt. It refuses to launch when its conservative worst-output reserve
would exceed the existing8MiB ledger or record caps. It performs no adaptive
splitting, cap increase or completion credit for timed-out partial work.

### Standalone profiles pilot, outside coordinator phase order

Root may run the following **one** bounded pilot when the unchanged live
resource gates admit it. The output ledger must be new; this command was
not run in this review.

```sh
python3 -B scripts/run_bounded_research_check_20261018.py \
  --parent-lock \
  --ledger results/finite57_distance9_profiles_0_8_pilot_20261020.json -- \
  scripts/check_selective_distance9_batches_20261020.py \
  --verify --phase profiles --start 0 --stop 8
```

A successful wrapper writes the complete direct receipt, whose schema is
accepted by coverage41f3. A failed/deferred pilot writes no successful
wrapper ledger and must remain failed/deferred in root's actual execution
record. For final coverage, supply this receipt explicitly with another
`--distance9-ledger` option alongside the runtime ledger, or have root
preserve its actual receipt in the coordinator's jobs ledger. The coordinator
does not automatically see the standalone file and must not skip profiles
merely because this command was proposed.
