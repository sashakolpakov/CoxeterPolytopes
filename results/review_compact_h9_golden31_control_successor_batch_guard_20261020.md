# Golden31 bounded-certificate control successor: outside code review

20 October 2026. Control-delta review only. PASS; no correction requested.
No import, scientific execution, generation or replay was performed.

Bound successor worker:
`scripts/check_compact_h9_golden31_unit_bisector_20261020.py`,
SHA256 `ffffc330b97510ea700aed057390cbe74bf3525ad6e364d16f3cd41169303576`.
Bound successor contract:
`results/compact_h9_golden31_fixed_checker_preparation_20261020.md`,
SHA256 `550a64e17a4008434060c46766c3cdeade08be575a9dbef86591e5a545711294`.

The previous complete mathematics/source review remains unchanged:
`results/review_compact_h9_golden31_unit_bisector_batch_guard_20261020.md`,
SHA256 `35a0ef17e17b65102620577714df405279030d8fcbe6d3cd791511d5a18b9619`.
The written mathematical proof remains unchanged at SHA256
`7259c05e994fd48d8bb5f70c172e98f9c0120bbfa13fc58f3a6b08174e5f9ea1`.

I independently read all changed control blocks. The new certificate
reader requests at most MAX_CERT+1 bytes and rejects an oversized result.
Verification snapshots the bounded saved bytes BEFORE recomputation,
then checks both exact new payload equality and equality with a fresh
bounded read AFTER recomputation. The old unbounded read_bytes path and
its time-of-check/size-check gap are removed. Neither read can allocate
an unexpectedly large certificate if the file changes during the check.

The actual conditional RuntimeError guard rejects disabled assertions
before parsing or checking; the additional assert is consistent, and
cannot be the sole protection under -O. The preparation SHA pin matches
the successor bytes. Owner enforcement, exclusive output creation,
postcheck source hashes, alarm and outer-runner resource requirements
are unchanged. The mathematical check body and literal-host scope are
unchanged. The contract describes the new snapshot behavior accurately.

This addendum supersedes only the prior worker/contract control-byte
binding. It does not assert an execution result or change the earlier
mathematical verdict. Root remains the sole guarded executor.
