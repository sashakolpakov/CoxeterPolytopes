# User-authorized completion of the Hexagon queue

The user explicitly authorized the side task to proceed with the remaining
submissions, reporting a five-per-day account allowance. The official API
schema was fetched again: version 1.13.0-beta, 36,413 bytes, SHA256
`46e761cf8bfafc54589139426026ab7875bae5f1a4f0edb48243cddda8c092f3`.
It exposes no account-quota endpoint. The local cap therefore follows the
user's account report; server quota checks remain active.

The existing calendar worker had already submitted the ideal-volume paper
as `2611.00958v1` at 00:06 UTC on 24 November. After reconciling that actual
receipt and verifying no unknown requests or unresolved markers, the client
was updated to permit five accepted commits per UTC day and remove the
obsolete different-day eligibility requirement. The approved-payload
snapshot, release manifest, metadata, source files, reference PDFs, unknown
mutation holds and preview requirements were retained. Original client and
configuration bytes are preserved in ignored owner-only local backups.

The LP companion was actually uploaded, previewed, compared and committed
once as `2611.00962v1` at 02:00 UTC, with HTTP 202. The remote 36-page
preview has exactly the retained PDF's normalized full text and bibliography,
BibTeX ran, and the final pass has no undefined references. Only the two
approved TeX/BibTeX files were uploaded. The original rejected `.bbl` attempt
and previous accepted receipts remain intact.

The publication watchdog recorded 67.032 seconds for the LP pipeline,
81,952,768 bytes peak process-group RSS, exit code zero, and no stop reason.
Initial free-disk and available-memory fractions were 35.289075% and
36.950429%, with the unchanged percentage policy source hash
`51ff53ece114b805c50a79f1e69a72ba4f275296b4239679232fda2c3430114f`.
The shared scientific/typesetting lock covered the remote preview pipeline;
the parent coordinator kept its lane free. No scientific, Lean or local
TeX process was launched by this task.

The completed calendar label is disabled. Its actual macOS output says
`disabled`, rather than the older `true` representation initially expected
by the new public-report exporter. That export failure is recorded
separately and the parser was corrected. No submission was retried.
The accepted source/metadata bytes were unaffected. The client's preview
summary also omits private draft and preview identifiers from future stdout.

The [actual completion summary](../release/hexagon_submission_schedule_20261123/completion_five_per_day_20261123.json)
records a fresh remote status check: the sharp paper `2611.00913v1` is
published; the ideal-volume and LP papers are under review. Acceptance with
HTTP 202 is not a publication or correctness claim. The held geometric
dimension paper was not submitted because its sources still differ from
the recorded release/build package. No bulk payload, credential or third-party
reference PDF was uploaded.

This is an execution/source-delta account, not a new independent review of
the complete client. The earlier independent source review remains bound
to its historical client/configuration bytes.
