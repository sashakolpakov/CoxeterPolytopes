# Multi-anchor paid-foot kernel: actual execution and binding review

21 October 2026. **PASS for this module's saved bounded compilation,
six actual axiom outputs and current source/object/dependency bindings.**
This is an independent read-only evidence review. No compiler, Lean
import, axiom audit, scientific worker or prior kernel was rerun.
Only small standard-library metadata parsing, streamed file hashing,
archive inspection and a nonblocking lock probe were used.

## Exact current evidence

| File under results/ unless otherwise stated | SHA256 |
|---|---|
| formal/CoxeterProofs/MultiAnchorPaidFootKernel.lean |801ccfc1dd19afee9ab65451ddc63640d95efc7a4c0d802d6433dfe4ff897d4f|
| multi_anchor_paid_foot_compile_20261021.json |91e0777467c3738062f1939b8c510859e2a09ba6240daa73fb70db4b70d5e75e|
| multi_anchor_paid_foot_axiom_audit_20261021.json |fcb594f7342d9dce897e1bd9b442753676ced5fdddbaff35d3d6ac686ae64972|
| multi_anchor_paid_foot_actual_axioms_20261021.txt |210c96f6defca314f2d1ee7a8e3ad11626424158dc9f6550049734cf2524290a|
| multi_anchor_paid_foot_axiom_summary_20261021.json |39f1f40d0c7cfac321290122197afe2e1e815a3ffea9f8524637efe2fd5aa7d1|
| review_multi_anchor_paid_foot_source_replay_reviewer_20261021.md |119e195b91cdc54ee2b37e6a9f66a3c5665e4b33bc27ac03fd11986c84add7ed|
| docs/multi_anchor_paid_foot_logic_20261021.md |f130e546fec57e788755e16a8b6ab999fff11224e11172cf59466dc236de6ac9|

The 7989-byte Lean source and the pre-execution review are unchanged.
The compile and audit snapshots are byte-identical to their respective
operational ledgers under formal/.proof-batches. They bind the same
source, compiler/configuration/environment and five local dependency
objects. Their status fields are KERNEL_COMPILED and KERNEL_AXIOM_AUDIT.
The reviewed runner publishes these statuses only after successful child
exit and unchanged input snapshots. The target object exists and its
streamed digest is exactly

    a1fea19989fe3a235442cca91b5557f37955a7adb699623aa5e5b733a1bf59aa.

## Actual printed declarations, not a reconstructed expectation

I parsed the entire 772-byte raw output, including wrapped axiom lists,
and checked there is no extra unparsed nonwhitespace content. Its six
unique names occur in exactly the order of the six source commands:

* payment_coefficients
* weighted_payment_decomposition
* joint_anchor_height
* multi_anchor_paid_pressure
* multi_anchor_paid_energy
* multi_anchor_paid_contradiction

Each is in namespace CoxeterProofs.MultiAnchorPaidFootKernel. Every
actual list is exactly {propext, Classical.choice, Quot.sound}; none
contains sorryAx or any additional axiom. The parsed summary reproduces
every full name and axiom list, source/raw digest and output byte count.
The raw digest and772-byte length match BOTH compile and audit ledgers.
Thus the six printed outputs are backed by retained actual bytes, not
just a list prepared before execution.

## Local dependencies and compiler identity

The local import closure contains exactly the following five modules;
I streamed each current object and checked its own current capped compile
ledger, source digest and complete local dependency map as well.

| Local module after CoxeterProofs. | Current object SHA256 |
|---|---|
| RightAngledSharpDimensions |f2bd3858fee23605e95212d21784afcc8f6a66f148e11f524ea3cbeb077d7f91|
| RightAngledCubicScalarKernel |8d06c91e92439233c0469782bd382f5cb87fc7056488800e40298bf4e00c781c|
| RightAngledProjectionHeightKernel |b2e265c8445767735f7ea149cfcf8344df074229bf62bafe158b0b482d7ef2a3|
| SummedActualFootKernel |4010db9c32b04e5abad9938d2e61dc72ba3eb93cf4092c4b3b5b74f3ecb8a468|
| AnchorAndPaidFootKernel |b36389f191b1df69ec8f4f2a8c5c441d53d68007408fe938b197603bc6486bdd|

All five have KERNEL_COMPILED provenance with2048MiB/one thread and
matching current inputs. No dependency is supported merely by a legacy
uncapped cache label. This review does not claim they were newly executed.

The compiler snapshot matches current lean-toolchain, lakefile.toml,
lake-manifest.json and all four listed binary files: the pinned
Lean4.32.2 executable and three shared libraries. The recorded ELAN_TOOLCHAIN, LEAN_PATH
and LEAN_SRC_PATH values also match. This is the recorded toolchain and
local-object closure; it is not a separately hashed closure of every
transitive Std object.

## Resource and owned-cleanup boundary

The unchanged signal-safe wrapper has SHA256
`c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d`;
its base runner has SHA256
`97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c`.
I read the launch, lock, input-snapshot, child-exit and cleanup paths.
The command uses Lean -M2048 and -j1. The runner holds the shared
exclusive lock, checks parent ownership and continuing resource floors,
enforces its timeout and reaps its owned child group. No watchdog or
resource gate is disabled by this module.

The compile ledger records2.19s and preflight158.24GiB free disk /
6.75GiB available memory; the separate audit records.38s and
158.23GiB /6.55GiB. Both exceed the40GiB/6GiB admission gates.
These ledgers record Lean's memory setting and preflight availability;
they do not contain a measured peak process RSS, which is not asserted
here. Successful record publication and the preserved actual output are
consistent with the reported terminal zero exits.

The owned formal/.proof-batches/work path is absent, the shared lock
can be acquired and immediately released, and the target's only build
artifact is its retained .olean. I removed nothing and did not inspect,
terminate or clean any foreign workspace process or work directory.

## Map, counts and mathematical scope

The original preparation map is preserved at
results/source_archive/multi_anchor_paid_foot_map_before_actual_compilation_20261021.md.gz,
archive SHA256
`b337bcdc0d611143959343a4d655274e9f833e096e4cd17d82785a1f8cf3afe7`.
Its decompressed bytes hash exactly to the original b909940f map bound
by the unchanged source review. The current map adds the actual records
and accurately identifies the one-end A20 coefficient application.
I flagged a stray plus sign and an indexing typo before final freeze;
both are now corrected. In the zero-based array Omega1..Omega20,
Omega6 is entry5 and Omega16 entry15. No Lean source changed.

Exactly six newly checked conditional declarations are added. The stated
aggregate12 kernels/99 declarations inherits the previously accepted
11/93 and adds these six; this audit reruns none of that prior work and
does not replace the preserved historical seven-kernel master receipt.
Actual geometry, available old floors, feasible foot pressure, reference
flatness, the full residual square and concrete strict scalar estimates
remain application inputs as explained in the source review. The
execution formally checks the generic internal payment/energy logic,
not whole A20 or an end-to-end geometric dimension theorem. No new
scientific replay, bound or example follows from this verification.
