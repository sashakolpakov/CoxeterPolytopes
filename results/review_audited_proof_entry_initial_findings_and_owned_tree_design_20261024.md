# Audited proof entry: initial findings and owned-tree design

24 October2026. Source-only review and implementation note. This record
preserves the first review's findings and describes the new owned-tree
implementation. It is **not** a mathematical proof, a successful execution,
or final publication clearance for the controllers.

No worker, scientific import, enumeration, fixture, test, Lean, TeX or Sphinx
job was run by this author. The first three unfrozen controller versions are
byte-preserved in
`results/source_archive/proof_entry_first_static_review_20261024/manifest.json`.
Every gzip was decoded and compared byte-for-byte. They were not silently
given a passing review after their defects were corrected.

## First source versions and findings

|First reviewed file|SHA256|Bytes|
|---|---|---:|
|`scripts/run_proof_batch_20261024.py`|`7e2ea3420eeb86b41909a95b4a34917cf3162042ba9f48185644d949dbac3dd9`|3731|
|`scripts/audit_existing_lean_module_20261024.py`|`eea8b0c6a327924ebfe06c7a6d837fb63172a004c2d41706b652d316a8289612`|2956|
|`scripts/verify_manuscript_audited_20261024.py`|`cd37938fb41feebc6637e9efa8a53316ab1804f56b5cae27b36b735dc5ed6cc2`|13065|

The initial review requested these corrections before source freeze:

1. The audit helper checked that literal `#print axioms` commands existed,
   then accepted a zero-exit Lean invocation without parsing its actual
   output. Missing results, duplicate results, unapproved axioms or sorryAx
   therefore did not prevent the new helper's PASS status. Validation must
   occur before the inherited batch publishes an audit ledger.
2. The direct public `check` route passed arbitrary wrapper options through.
   Omitting `--parent-lock` for a nonlocking endpoint acquired no shared
   lock. The public route must default to the outer lock, with only a
   small explicitly audited set of self-locking workers exempted.
3. The successor inherited `run_job` cancellation which terminated only its
   direct Python child. Nested Lean starts another session/process group.
   The outer runner's0.3-second group cleanup could kill the inner Python
   monitor before its5-second compiler cleanup finished. The separately
   grouped compiler could then outlive its monitors. Direct public checks
   also needed an independent owner-death cleanup path.
4. Per-job provenance originally bound only the new public entry plus the
   target source. The new audit helper and other controller implementations
   needed an explicit before/after source closure too.

During correction, exact comparison of raw short axiom queries with Lean's
fully qualified result names was found to reject valid existing modules.
The old ResidualBellman and ConcentrationDescent sources supply concrete
counterexamples without executing Lean. One-to-one unambiguous name
resolution is needed. Dropping empty comma-separated axiom tokens is also
too permissive for a malformed-output rejection contract.

These are execution-contract findings. They do not establish that any
historical proof run was optimized or that its mathematical conclusion is
false. Frozen workers, legacy runners and the old driver must remain
byte-exact; the25 current producer-equivalent units pin the old driver.

## Authorship and review boundary

Root authored the new public entry, audit helper and initial successor
driver. This author supplied the fresh static findings above. At Root's
request I then authored
`scripts/monitor_owned_proof_tree_20261024.py` and the new driver's
owned-tree cleanup and controller-closure integrations. My review of those
changes is a self-review. A different author must review that implementation
before final acceptance; my earlier findings do not count as an outside
review of my later code. Root owns final integration and any admitted tests.

## Owned-tree implementation

The new module is both a small standard-library monitor and a command-line
bridge. It imports no scientific worker. The actual child remains an
unchanged bounded research or TeX runner; their admission, RAM/disk floors,
RSS/time checks and shared-lock rules are not replaced or weakened.

The bridge accepts only the current Python executable with `-E -B` and one
of the three explicitly named legacy bounded runner/builder paths. It
preflights process-table and birth-identity access before creating a child.
The child starts in its own session but waits behind a pipe gate. It cannot
execute a legacy runner until the bridge has recorded its birth, created
the ownership tracker and checked that its own owner is still alive. Gate
closure on failure exits before proof work.

On Darwin, birth identity uses `proc_pidinfo(PROC_PIDTBSDINFO)` and both
start-time seconds and microseconds. The layout was checked against the
installed primary SDK header `sys/proc_info.h`: struct size136, PID/PPID/PGID
offsets12/16/100 and start-time offsets120/128. Permission or ABI failure
raises an error instead of falling back to weaker timestamps. Linux uses
the process start ticks from `/proc/PID/stat`. Other platforms are refused.

The process table supplies candidate parent/group relationships. A new
ownership record requires a presently verified owned parent or membership
of a still-trusted owned group. Previously recorded births remain tracked
after reparenting. A reused numerical PID does not inherit the former
birth's ownership. Group leaders are recorded with their birth identities;
a replaced numerical leader invalidates the group. The coordinator's own
group and IDs at most1 are never signal targets.

Before a group signal, at least one observed member's birth and group are
read again, and the group is revalidated. TERM and then KILL are ordered by
the largest observed descendant depth, so inner compiler groups are
targeted before their Python monitors. Cleanup keeps sampling and therefore
retains observed children even if an old monitor exits during cancellation.
If verified live owned processes remain, cleanup fails. A nominally
successful child with observed stragglers is cleaned up and then rejected;
that is not converted into a proof PASS.

The bridge checks owner identity throughout execution and has an additional
59-second elapsed limit. Normal exit, failure, owner disappearance and
cancellation all pass through cleanup. The new manuscript driver's
`run_job` independently tracks the direct public-entry tree and cleans it
on output-cap failure, nonzero exit or interruption. Repeated cancellation
is ignored only while that bounded cleanup is in progress. Direct-child
fallback applies only to the unreaped direct Python child when tracking
could not be initialized; public-entry descendants have their own bridge.

## Explicit limitations

This is a polling monitor, not an atomic kernel process capability. A
hostile short-lived process could fork a new-session child and disappear
between samples before the relationship was observed. Darwin also lacks
the Linux pidfd interface used for atomic identity-bound signaling; the
birth check and numerical group signal have a kernel scheduling interval
between them. The implementation must not be described as proving absence
of every adversarial fork/reparent or PID-reuse race. It conservatively
refuses ambiguous/reused identities and never deliberately signals an
unidentified group.

The concrete existing Lean hierarchy has a Python parent alive while its
compiler runs. Tracking its observed compiler group across the parent's
termination addresses the identified0.3/5-second cleanup race. Positive
runtime interruption/straggler tests and the platform metadata probe remain
separate evidence. None were run by this author. A parser or negative
preflight alone is not a runtime stress test of the owned-tree implementation.

## Controller source closure

The shared `execution_sources()` function binds four new files: public
entry, audit helper, successor driver and monitor. It also binds the six
legacy controllers: bounded research runner, signal-safe Lean runner, base
Lean runner, both allowed TeX builders and the old verification driver.
Each source is hashed with a before/after filesystem identity check.

The bridge captures that ten-file map before its first child and compares
it after cleanup. Its receipt records the map, assertions-enabled flag,
effective optimization flag and cleaned child optimization environment.
The new driver records the same map in plans, per-job records and its
overall receipt, checking it before and after every job and at the end.
This is controller provenance, not a separately hashed Python/OS/Std
dependency closure or proof of any mathematical input.

## Plan scope and remaining evidence

The successor retains the same136 finite57 script selections and9 Lean
roots. Exactly five roots contain actual printed-axiom commands. Splitting
their audits creates151 full jobs and16 core jobs, rather than the older
146/11 inventory. This is a source-derived plan count, not execution.
The four roots without commands must not be called actual axiom audits.

The new actual audit must require an already matching capped target object,
refuse dependencies needing rebuild, obtain the shared lock, run separately
and validate every requested result before ledger publication. Ordinary
theorem parameters remain external even if only standard logical axioms
are printed. No source correction adds a geometric formalization or a new
mathematical bound.

Final source hashes, a complete other-author review of the monitor/cleanup,
negative preflights and admitted positive executions remain separate
requirements. This design note supplies none of those execution receipts.
