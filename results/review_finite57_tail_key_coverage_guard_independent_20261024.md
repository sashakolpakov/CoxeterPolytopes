# Independent source review of the finite57 tail key supplement

24 October2026. Reviewer `/root/audit_finite57_algorithms`; author of the
guard `/root/consolidate_six_scope`. No authorship of the reviewed guard
is attributed to this reviewer. Source-only review: no import, execution,
test, certificate arithmetic, Lean or TeX job ran.

## Current successor verdict

The corrected338-line guard, 16,559 bytes, SHA256
`bdda6b5cd4b4da3894f8981d5c67fd226a595dad18cfb963663f4b48244f3b1c`,
is source-clear at its stated **key-only** contract, pending root's
coordinating acceptance. The bounded-read defect and unit two-contact
scope issue below have been corrected. The observed-input-hash boundary
remains explicit and requires external clearance/receipt binding.
No test or positive execution result follows from this source review.

The original source and contract are preserved byte-exact in
`results/source_archive/finite57_tail_key_coverage_supplement_pre_bounded_read_20261024/`.
The initial independent review is separately preserved, with a verified
gzip byte round trip, by
`results/source_archive/finite57_tail_guard_independent_review_first_version_20261024/manifest.json`;
its SHA256 is `49593e1011f4df3bfc614dd494617f4763abae3704fa3f45e1e91114f7b87b05`.
The final successor review appears after the original findings below.

## Original-version findings (preserved history)

The complete324-line source reviewed is
`scripts/verify_finite57_tail_key_coverage_supplement_20261024.py`,
15,629 bytes, SHA256
`2480998ff71620dd9413abd63d7c3337f341f90c9cbcd40a2d406782e797a137`.
The complete contract is
`results/finite57_tail_key_coverage_supplement_source_plan_20261024.md`.
This review also inspected the legacy unit-charge two-contact producer and
reviewer branches and the classical-coordinate review loop to ascertain
the supplemental binding's actual limits. It does not repeat root's
separate mathematical audit of the unbounded tail.

The literal key sets, type checks, duplicate rejection and stated
cross-bindings are sound at their explicit inherited reference-domain
premises. No numerical or global verification overclaim was found.
One bounded-input implementation correction is recommended before a
positive launch: replace `stat` followed by unbounded `read_bytes` with
an actually bounded read. Two provenance/coverage limitations must remain
explicit, as described below. No mathematical counterexample is alleged.

## Strict independent key domains

`require` uses a real conditional and raises `CoverageError`; no assertion
can disappear under optimized Python. Every row index first checks list
and object shapes, rejects duplicate computed keys, and then compares the
complete actual key set with the expected set. Duplicate JSON object keys
are rejected during parsing as well. There is no silent dictionary overwrite.

The four E7 keys include side, a strict boolean own-omitted state, and both
ordered support tuples. `integer` requires `type(x) is int`, excluding
Python booleans, while `flag` requires `type(x) is bool`. Support entries
are distinct strict integers. Consequently neither boolean/int aliasing,
four repetitions of one valid row, an exchanged side, nor replacement of
one omitted-node state by another can satisfy this set comparison. The
four literal expected tuples agree with the stated two ordered sides and
two own-node states. This is stronger than a four-row count.

The cycle wrappers require exactly the sixteen integer sizes2–17, each
once, in all three primary/review wrapper arrays. Inner support tuples
are nonempty, have distinct nodes and lie inside that cycle's range.
The widened scalar receipts are then required to have exactly the580
cross-bound `(size,left_support,right_support)` keys, once each.

The exceptional descriptors are independently constructed for ranks
1,6,7,8,4,2, using `None` and all r+1 affine-node avoid states. This gives
40 descriptors. All seven relevant primary/review descriptor arrays must
match that complete set without duplicates. The expected singleton,
multiple, subunit and unit singleton key sets are then cross-bound
strictly at each descriptor. The71/16/46 totals are supplemental metadata
checks; they are not presented as proofs of geometric completeness.

The fifteen failure occurrences are indexed by full descriptor plus the
ordered contact-support tuple. Both prescribed E7 same-row occurrences
must be present, and the thirteen remaining floors must be precisely the
set difference. The two rows cannot be removed merely by a count filter.

The three packing budget keys are ordered string pairs. `(1,3/2)` and
`(3/2,1)` are distinct required keys. The five B/C/D analytic-tail keys
include strict boolean modified status and integer start rank; the two
A keys include parity and the complete start-form tuple. The24 classical
coordinate classes preserve the distinction between integer side0/1 and
symbolic string sides `r`/`r-1`. Boolean sides are rejected. The seven
four-thirds cycle-class keys similarly include the complete allowed-distance
tuple. All these arrays require set equality, not containment or counts alone.

## Cross-bound reference assumptions

Inner expected support keys for small cycles and exceptional factors are
taken from named primary or independent-review outputs. The supplement
checks uniqueness and exact agreement; it does not regenerate those
mathematical domains. A coordinated omission from a reference and every
bound output could therefore pass. The contract states this correctly:
legacy domain reconstruction and numerical verification remain premises.

The `two_contact_case` unit-charge check is narrower still. It requires
that exactly descriptor `(E7,None)` has a non-None record, but does not
bind that record's internal `nodes` field to `((0,4),(1,7))`. The original
producer explicitly constructs those contacts, and its independent reviewer
reconstructs the unique cheap support and its polynomial, while checking
the saved record's presence. Thus the inherited mathematical reconstruction
is not missing, but this supplement must not be described as a full
inner-key binding of that single two-contact record. An optional literal
node-key comparison would make the supplement stronger without any new
scientific computation.

Similarly the guard does not validate ignored numerical fields, polynomial
coefficients, KKT slacks, model matrices, packing predecessor arithmetic
or omission choices. Its `PASS_FINITE57_TAIL_KEY_COVERAGE_ONLY` status and
returned scope correctly exclude such claims.

## Input-size correction and execution/provenance boundary

The reviewed `load` first checks `path.stat().st_size<=32MiB`, then calls
`path.read_bytes()`. A file can grow or be replaced between those operations.
The latter read is consequently not bounded by the advertised limit. This
is a resource-hardening defect, not an invalid key comparison. The minimal
fix is to open the file and read at most `MAX_INPUT_BYTES+1` bytes, reject
an over-limit length, then parse/hash precisely those captured bytes. A
post-read length test alone would detect oversize after the allocation,
but would not enforce the intended allocation-read bound.

The source imports only the standard library and launches no child process.
There are no project-worker imports or hidden compilation paths before its
guards. It does not acquire the shared proof lock or monitor resources
itself; the contract requires the coordinating admitted bounded wrapper.
Through the existing public `check` route it is an ordinary script and
must receive the default parent lock, not a self-locking exemption. Direct
invocation is not evidence of compliance with that resource contract.
This review did not run the wrapper or test the lock choice.

The returned `sources` map records SHA256 of the exact consumed JSON byte
buffers and the guard source. It does not compare those input hashes with
an independently fixed expected manifest, and it does not recheck files
after processing. A future receipt therefore needs an outer source/input
identity check against the reviewed frozen bytes before reuse. Alternatively,
the guard can enforce a fixed expected-input manifest. Recording hashes
alone is provenance, not acceptance of arbitrary replacement certificates.
The contract's stated key-only success remains valid for the captured bytes.

All loops are finite iterations over parsed arrays or literal finite sets;
there is no scientific recursion or optimization. Malformed rows fail by
explicit checks or ordinary exceptions and cannot produce a successful
result before all checks have returned. The bounded runner remains necessary
for aggregate JSON memory, parsing depth and runtime control; a per-file
byte limit is not a process heap cap.

## Recommendation

Fix the bounded read, keep the reference-domain and unit two-contact scope
limitations explicit, and bind any future runtime receipt to its exact
seven input hashes and final source hash. Root should review the changed
bytes before admitting the minimal key-only check. This report supplies
no test or execution result and does not authorize a raw proof-process
launch or a numerical/global finite57 PASS.

## Final review of the corrected source and prepared structural driver

The complete old-to-new source diff contains only the following changes:
bounded cap-plus-one reads and explicit actual-length tests for input JSON
and the source hash; exact ordinary-E7 two-contact support recognition
with every other descriptor required null; and an explicit observed-hash
policy. I read every changed line and its surrounding control flow.
The archived original hashes to2480998f as reported. No numerical checker,
mathematical domain, expected descriptor or worker was altered.

The actual-read fix now bounds the byte allocation before parsing, even
if the preexisting stat result becomes stale. The new E7 branch requires
an object with the strict ordered `support_key=((0,4),(1,7))`; all other
descriptors require `None`. Thus the two initial findings have been
resolved in the successor. The unchanged reference-domain limitation is
accurately stated in both code and the revised contract, SHA256
`fe8b26013e77a2920dc7f59778887275668679731d2bb38833259aefabbeacda`.

I also read the complete prepared235-line structural regression driver,
`scripts/test_finite57_tail_key_coverage_supplement_20261024.py`,
11,368 bytes, SHA256
`3e9e47313d14cfb214f6e1fdce2e027188c386ddd13e80c3f6457c6652de08c1`.
It is source-clear for its narrow regression contract, pending root
acceptance and resource admission. It has not been imported or run here.

The driver reads at most four MiB plus one for each file/source and
checks length, pins the seven original JSON inputs, and caps the owned
fixture total at eight MiB. It executes exactly the captured/hashed new
guard bytes into a synthetic module whose name is not `__main__`.
This loads standard-library definitions and no mathematical worker or
subprocess. Its first positive result must bind all seven inputs and the
same guard hash as the captured source. The eleven negative fixtures
must raise the guard's own `CoverageError` containing the specific
expected clause; unrelated exceptions and accepted invalid fixtures
cannot be credited as successful rejections. A whole valid packing-row
permutation is an allowed positive case. Reassociation of numeric fields
with an unchanged key set remains outside this key-only contract.

Every mutation writes only inside one `TemporaryDirectory` owned by the
invocation; a `finally` restores the changed file, and the outer context
removes the owned fixture even on a failed case. Real frozen inputs are
read-only. A failure during writing is likewise inside the cleanup context.
The receipt is constructed only after that context exits. A launch still
needs the ordinary public `check` route, default shared parent lock and
unchanged resource/watchdog gates; this source review does not exercise
or waive them. Source/input before-and-after identity remains an outer
clearance obligation, especially before a result is credited to a fixed
proof chain.

The reviewed sources contain two positive and eleven negative structural
cases, not thirteen mathematical proof checks. There is no source-level
reason to reject this corrected package at its explicit reference/numeric
boundaries. Root's acceptance, an actual admitted execution and numerical
tail verification are all distinct further events.
