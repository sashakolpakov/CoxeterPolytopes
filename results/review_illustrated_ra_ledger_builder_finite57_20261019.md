# Independent source review of the ledgered manuscript builder

19 October 2026. **PASS for the intended serialized local build contract.**
Read-only source review; no TeX, subprocess, resource probe, AST/checker body
or build was executed by this review. No source edit is requested.

Reviewed all154 lines of
`scripts/build_illustrated_ra_manuscript_with_ledger_20261019.py`, SHA256
`ae521225307be682f120860c1e0a90279b206b24191a084be19ccedec0c12c8e`,
and the relevant resource/cleanup functions of its imported signal-safe
runner and resource helper.

The source and ledger paths are restricted to local manuscript TeX and
results JSON. A preexisting ledger is rejected, and successful publication
uses exclusive-create mode. An existing PDF is refused by default. The
explicit `--replace-unfrozen` override is allowed only after checking the
PDF path against all existing `*freeze.json` manifests; frozen evidence
therefore cannot be replaced through that option. For a historical PDF
that must remain unchanged, use a different manuscript stem or
`--check-only`, not the replacement option.

The builder itself owns the shared nonblocking exclusive flock. It applies
the existing40GiB/6GiB start gates before spawning and the8GiB/3GiB runtime
gates while the child is active. The process-tree RSS probe must work
before launch. Runtime checks enforce the one-GiB sampled process-tree
limit,59-second deadline, bounded512KiB child output and parent identity.
Numerical library thread variables are fixed to one. The three TeX passes
are sequential and stop after at least two when the rerun request clears.

Every spawned TeX group is stopped and reaped in `finally`, using the
signal-safe helper's final group kill even if its leader has exited.
Temporary compilation is isolated in a `/private/tmp` temporary directory,
whose context removes intermediates on success or failure. The input
source hash is checked after compilation. The resulting ledger binds the
source, final PDF and builder hashes, page count, final warnings, pass
count, resource/time measurements and the explicit absence of any
mathematical-verification claim.

The status means a successful TeX compilation, not zero warnings or proof
verification. Warnings are retained rather than implicitly waived.
`--check-only` prints a result with `check_only=true`; it does not publish
the temporary PDF or create the requested ledger. Failed/prelaunch-deferred
attempts likewise do not create a success ledger. The coordinator must
keep their separate provenance if needed. These boundaries agree with the
source and do not obstruct the intended new-PDF build.
