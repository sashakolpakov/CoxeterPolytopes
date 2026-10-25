# Corrected public proof entry: source review

24 October2026. **PASS for the corrected Root-authored public entry and
actual-axiom helper at the source boundary below.** The plan-generation and
original orchestration portions of the successor driver were also reviewed.
The new owned-tree monitor and its driver integrations are my implementation:
their treatment here is explicitly a self-review. Root's complete
different-author review and coordinating acceptance of that part remain
separate records.

No scientific worker, Lean, TeX, mathematical import, enumeration or test
was executed by this reviewer. I read all four new controller files, the
complete14-test infrastructure file, the relevant unchanged runner/cache
functions, and the seven self-locking public entry blocks. I also read
Root's actual infrastructure stdout. Source pins are in
`results/review_audited_proof_entry_source_bindings_20261024.json`.

|Corrected artifact|SHA256|Bytes|
|---|---|---:|
|Public entry|`78bfafca1791fcdfc217860ea923cb15c16fbd6eb0fe6943586e41849121e4bb`|6319|
|Actual-axiom helper|`1399fbc6fbc7f0c92496b7f7284fcc8389573c544ee8bb427739666c1af6f7d5`|6370|
|Successor manuscript driver|`54e57bdb9cd0e9a9c3bb160743f22c4da1728c353f1e5a1cf9db1b397e8ef968`|15406|
|Owned-tree monitor|`517c0c78347c0aa81cd16489afcf7f11b7bce064490919b9cf5ff94c88e59cda`|15758|
|Infrastructure tests|`cebc970ae53f0da292337397f808e594538903ac2fdd484b839c48d618519915`|12446|

The earlier three versions and their original findings are preserved in
`results/source_archive/proof_entry_first_static_review_20261024/manifest.json`
and
`results/review_audited_proof_entry_initial_findings_and_owned_tree_design_20261024.md`.
This corrected verdict does not erase those findings or give their original
versions a retrospective PASS.

## 1. Optimization and public routing

Every new executable controller has an ordinary conditional rejecting
`not __debug__` or a nonzero `sys.flags.optimize` before repository imports.
The entry also rejects a nonempty/nonzero `PYTHONOPTIMIZE` setting. These
checks survive Python optimization; they are not assertions. Child
environments have the setting removed and use the current interpreter with
`-E -B`. Every scientific Python descendant of the old bounded runner
therefore inherits a cleaned optimization setting. Existing old assertions
remain enabled along the preferred route.

The entry constructs argv arrays and uses exec, not a shell. Exec preserves
the public entry's parent relationship. Actual check, single-module Lean
compile, separate Lean audit and TeX modes first enter the new owner bridge,
then an unchanged bounded runner. The manuscript route enters the new
coordinator directly; that coordinator launches each actual child through
the same public entry. List/plan modes return an inventory without launching
a proof child.

The proof runner's40GiB/6GiB admission,8GiB/3GiB abort floors,
59-second/1GiB sampled process-tree limits and one-thread environment are
not lowered. The new bridge adds another59-second deadline; it does not
replace those checks. TeX routes admit only the two named existing bounded
builders. The additional monitor's process-metadata preflight is not a
proof-admission substitute. If process monitoring is denied, it fails before
the gated runner is released; it does not disable monitoring.

The scientific worker source and its mathematical-input bindings remain
the old runner/checker's responsibility. The new controllers neither
generate expected certificates nor reinterpret a failure as a success.

## 2. Shared locking, including the hidden worker counterexample

The direct check route requires an explicit `--` boundary and an existing
script directly under the repository's scripts directory. It defaults to
the outer shared lock. The only exemptions are seven exact source/argument
tuples, matching the current catalogue's `owns_lock` entries. I inspected
their public main blocks: the selected tuple reaches their inner lock
before the actual requested scientific calculation. Explicitly requesting
the outer lock as well for an exempted tuple is rejected.

A filename-only exemption was insufficient. The self-hub reconstruction
source has an internal `--worker` branch before its own lock. The corrected
entry exempts that source only with the exact public `--verify` argument.
The internal branch and every other nonmatching tuple receive the outer
lock. This closes the concrete bypass without editing the frozen worker.

Single-module Lean modes require a named `CoxeterProofs.X` target and reject
the whole-project root. They deliberately omit the outer `--parent-lock`
because the compiler or actual-axiom helper obtains that same shared lock
internally. The script remains inside the outer RSS runner. The successor
plan therefore never intentionally acquires the lock twice or removes the
inner lock to avoid an error.

## 3. Actual axiom audit is a separate operation

The actual-axiom helper requires the bounded-runner owner relation before
work, checks the paused marker and obtains the shared lock. It asks the
unchanged cache verifier for the complete local dependency order. The target
must have a matching `KERNEL_COMPILED` receipt/object; any dependency needing
rebuild is rejected. The helper never calls the compiler's main loop and
never falls back to compiling a missing target. Legacy-cached dependency
provenance, where accepted by the old cache verifier, is reported explicitly
and is not relabeled as freshly capped dependency compilation.

The old cache verifier binds target source, actual dependency object hashes,
compiler/configuration/environment identity,2048MiB memory setting,
one-thread setting and target object hash. Consequently a stale target
receipt or a changed dependency cannot silently satisfy the helper.
The audit uses the signal-safe cleanup function and executes the unchanged
batch with `audit=True`; no object-compilation flag is added.

The helper captures the original `base.guarded_run` function before
wrapping it. Its wrapper validates the complete bounded raw Lean output
before the inherited `batch` function can publish its audit ledger. It does
not call `safe.guarded_run` through the replacement and therefore avoids
recursive wrapping.

Actual result validation now checks:

- a supported identifier in every requested `#print axioms` command;
- no duplicate request, actual result or axiom token;
- no missing, additional or ambiguous requested declaration;
- multiline dependency lists and the no-axioms form;
- rejection of malformed canonical result headers, empty comma items,
  malformed identifiers, sorryAx and any axiom outside `propext`,
  `Classical.choice`, `Quot.sound`;
- no truncated output.

Fully qualified queries match exactly. Namespace-local queries must resolve
to one unique actual fully qualified suffix; the resulting full names are
recorded. This is required by existing ResidualBellman and concentration
sources. A raw short/full set comparison would incorrectly reject them.

Only after these checks does the helper print
`PASS_SEPARATE_ACTUAL_AXIOMS_PARSED`. It records that no compilation was
performed and that ordinary theorem parameters remain external. Such
parameters do not appear in `#print axioms`, so the helper does not claim to
validate the geometric or enumeration premises.

## 4. Ownership and controller provenance

The initial owned-tree design note explains the concrete nested-session
cleanup defect, gate, Darwin birth metadata, conservative group ownership,
TERM/KILL ordering, retained reparented children and failure semantics.
Those implementation details were supplied by this reviewer, so their
correctness needs Root's independent critical review.

In self-review, I checked every normal/error/finally branch. No runner is
released before process metadata and initial ownership are available. A
zero-exit leader does not prove all descendants exited; observed stragglers
are cleaned and then cause a failure. Repeated cancellation does not skip
the new driver's cleanup. Controller source changes fail the job. An
ambiguous or recycled identity is not deliberately signalled.

The code and design note disclose the polling and check-to-signal limits.
It is not an atomic pidfd facility and does not prove resistance to every
adversarial instantaneous fork/reparent race. Zero remaining *verified
observed live owned processes* must not be expanded into an unconditional
claim about all processes on the machine.

`execution_sources()` binds all four new controllers and six legacy
controllers. The bridge checks that map before its first child and after
cleanup. The driver records it in metadata plans and execution receipts,
checking before/after every job and at final completion. This is a
controller-code closure, not a Python/OS/toolchain transitive dependency
claim. All frozen originals remain at their old paths and byte hashes,
including the old verification driver pinned by existing B13/frontier
workers.

## 5. Infrastructure evidence and its limits

Root's retained stdout
`results/proof_entry_infrastructure_preflight_20261024.stdout` records
**14 passing infrastructure tests in1.107seconds**. I read the complete
test source and stdout. The test domains are:

- optimized-startup refusal before downstream work;
- default/exact-tuple locking, double-lock refusal and hidden-worker case;
- Lean outer-monitor routing and whole-project rejection;
- plan-only full151/core16 phase counts;
- real preserved axiom stdout plus positive/missing/duplicate/ambiguous/
  malformed/unapproved/sorryAx cases;
- synthetic reparenting and recycled-PID ownership cases;
- two actual tiny owned Python stand-ins in different sessions, both
  ignoring TERM, cleaned with observed inner-before-outer signalling.

The live fixture is now gated before it can create its leaf. It constructs
the ownership tracker before release, bounds readiness through selectors,
and kills only the still-unreaped gated direct child if initialization fails.
The earlier unsafe fixture fallback is not the executed current version.

These are infrastructure tests. They launch no scientific worker, Lean or
TeX and discharge no geometric, numerical-certificate or mathematical
algorithm premise. Testing the audit parser on preserved stdout is not a
new actual Lean audit. No successful unchanged proof module was rerun by
these tests.

## 6. Planning and claim boundaries

The current finite57 selection still contains136 scripts and9 Lean roots.
Only five roots have explicit actual-axiom commands. Splitting those five
audits makes151 full jobs and16 core jobs. The old146/11 plans and their
receipts remain historical evidence; they are not silently rewritten.

The successor retains catalogue coverage strings and explicit formal gaps.
A successful future run would prove only the declared checks under those
listed premises. Individual cached modules, partial command coverage,
plans, parser tests and controller hashes do not establish a combined
manuscript runtime PASS or complete transitive regeneration.

No remaining source blocker was identified in Root's corrected public entry,
axiom helper or planning logic at this boundary. Actual scientific launches
still require Root's separate owned-tree acceptance and fresh resource
admission. This review supplies neither a new mathematical bound nor an
end-to-end Lean formalization.
