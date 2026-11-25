# Hexagon submissions and local schedule

**Later LP correction.** The existing LP work was subsequently updated,
under a separate explicit user instruction, to **2611.00962v2** with the
concise Vinberg E7 remark. HTTP 202 was received at 03:32:58 UTC on
24 November; its immediate state was processing. The
[separate correction record](../hexagon_lp_vinberg_update_20261123/README.md)
preserves v1's source and receipts and documents the new exact payload.
The initial queue records below remain unchanged as submission history.

**Current status, 24 November 2026 UTC (23 November in Toronto).** All
three approved ready papers have been accepted for processing, and the
completed calendar queue is disabled. The
[actual completion record](completion_five_per_day_20261123.json) refreshes
their remote states:

| Paper | Version | Current state |
| --- | --- | --- |
| Sharp right-angled dimensions 8 / 4 | [2611.00913v1](https://hexagonmath.org/abs/2611.00913) | Published |
| Second ideal right-angled four-volume | 2611.00958v1 | Under review |
| LP certificates for 57 / 17 / 15 | 2611.00962v1 | Under review |

The ideal-volume paper was submitted by the existing calendar worker at
00:06 UTC on 24 November. The user's later explicit instruction authorized
up to five deposits per day, and the unchanged approved LP companion was
submitted at 02:00 UTC that day with HTTP 202. Its 36-page remote preview
matched the retained PDF's complete normalized text and bibliography.
No new mathematical or local TeX job was run. The geometric dimension
paper remains held for the source/build/package mismatch described below.

The current client removes the superseded one-per-day local delay while
retaining payload hashes, metadata declarations, private mutation locks,
unknown-result holds, preview checks and server enforcement. The official
API schema is unchanged and has no account-quota endpoint; the account's
five-per-day allowance is recorded as the user's current account report.
The old client and pinned configuration bytes remain in owner-only ignored
backups. Their original source-review and installation records remain
historical evidence. A report-export parser mismatch with macOS's
`disabled` enum was [preserved](completion_report_initial_failure_20261123.json)
and corrected without repeating a submission.

The remaining sections preserve the **original installation record**;
their one-per-day rule, future eligibility dates and early statuses are
superseded by this update and the linked completion record.

Prepared on 23 November 2026 as a separate publication task. The sharp
right-angled paper was **accepted for processing with HTTP 202** at
04:24:13 UTC as `2611.00913v1`. Its subsequent status is `under_review`;
publication has not yet been confirmed. Actual private receipts remain in
ignored owner-only `.local` files. The local calendar job for the other two
ready papers is **installed and loaded**. Its first invocation exited
successfully and deferred the next paper until 24 November UTC. See the
[installation receipt](installation.json) and [queue](queue.json).

The first paper's private preview and submission compile both completed.
The final status check also found all reported screening jobs complete; the
submission itself still had status `under_review`.
The preview has six pages, the same bibliography, and the same extracted
text as the approved PDF except for its automatic title date: the unchanged
`\date{\today}` source produced 23 November instead of 4 November.
The final BibTeX/LaTeX pass has no unresolved references. This comparison
does not rerun the mathematical checks or replace a rendered-layout review.

## API result

The [official API entry](https://hexagonmath.org/api/v1) and
[OpenAPI specification](https://hexagonmath.org/api/v1/openapi.json)
provide authenticated submission operations. The inspected specification is
version `1.13.0-beta`; its bytes and SHA256 are recorded in `queue.json`.

Automation uses an account-issued bearer credential scoped to
`submissions:write`. Hexagon explicitly says its ORCID browser-session cookie
is not a credential for external automation. The account holder has supplied
the local path `hexagon_key.local`. The repository-wide `*.local` ignore rule
excludes it and nested local settings files. No `*.local` files are tracked,
and no history entry for this credential path was found. Its permissions are
`0600`. Neither its contents nor any credential-derived value is included here
or in the upload manifest. Browser cookies and password stores were not opened.

Credential validity was established by the successful private draft creation.
The published API has no
credential-inspection or authenticated account endpoint; reading a public
endpoint would not validate the key. The client reads it inside its process,
sends it only as an authorization header to the documented HTTPS origin, and
never includes it in command arguments, receipts or exception output.

The API supports private, idempotent drafts, source uploads, private compilation
previews and a commit that queues screening. No future-date scheduling field
or endpoint appears in the complete published specification. A new account
has a one-deposit-per-UTC-day limit; private drafts last 24 hours. The
local queue therefore creates each draft on its execution day. Actual
publication follows moderation rather than a date chosen by the queue.
These are current service conditions, documented
in the [API specification](https://hexagonmath.org/api/v1/openapi.json) and
[submission policy](https://hexagonmath.org/about).

## Selected payloads

| Earliest submission date | Paper | Status |
|---|---|---|
| 23 November UTC | Sharp right-angled dimensions 8 / 4 | Submitted as `2611.00913v1`; under review. |
| From 24 November UTC | Second ideal right-angled four-volume | Scheduled locally; all four release hashes match. |
| From 25 November UTC | LP certificates for 57 / 17 / 15 | Scheduled locally; all four release hashes match. |
| Held | Geometric dimension bounds 57 / 17 / 15 | Current TeX and bibliography differ from the recorded build/package inputs. |

Each ready Hexagon upload contains exactly `main.tex` and `references.bib`,
sent individually. Hexagon rejected registration of the generated `main.bbl`
with HTTP 400 because that supporting file type is not accepted. Its original
bytes remain in the canonical release checks; none of the TeX, bibliography,
compiled bibliography or reference PDFs has been changed. The remote preview
must generate and display the bibliography correctly before commit.

The PDFs are local review references, not upload payloads. The queue retains
the original three preparation descriptors and records the two-file Hexagon
selection separately, including the actual rejection. This avoids
repacking archives. No local mathematical or typesetting job was run for this
task; the first paper's remote preview and submission compilation are recorded
separately.
The bulk research archives, Lean products, logs, third-party PDFs and private
material are excluded. Repository links remain in the public comments.

The geometric paper's current source is `c8e02529…`, whereas the release
manifest records `8f9e5181…`; its bibliography is `932bf447…` rather than
`1783c386…`. Its retained PDF still matches the recorded `9b42a06a…` PDF.
The held entry records the complete hashes. It must be reconciled by the
coordinating publication workflow before activation; this side task did not
compile it, replace it with old source, or modify a frozen manuscript.

Before submission, the public author-name query for `Kolpakov` returned zero
works and no next page. This checked that public filter at that time only,
not private drafts or name variants.

## Authorized metadata

The `.metadata.template.json` files contain exact manuscript abstracts,
mathematical categories and the user's completed declarations. The first
authenticated response accepted the contributor object naming Alexander
Kolpakov; no ORCID was invented. The same contributor metadata is used in
the remaining two templates.

The account holder selected **CC BY-SA 4.0** and **mixed human–AI text** for
all three ready papers. Their templates now contain the exact API values
`"license": "CC BY-SA 4.0"` and `"authorshipMode": "mixed"`.
The credential path is also supplied. The user chose **some** for human
understanding and confirmed correctness, mathematical interest, citation
effort, distribution rights and the hosting licence for all three papers.
The provenance statement records referee-level review and distinguishes
ideas supplied by the author, whose provenance is known, from ideas
originating with the AI, for which the author cannot always trace how the
model arrived at them. This distinction concerns provenance; it does not
identify proofs as unchecked or poorly understood. The statement does not
claim complete formal verification or guarantee correctness.

No arXiv identifier is invented: these packages have not been submitted to
arXiv according to the repository's current submission record.

## Request sequence

1. Recheck the current API, account allowance, completed metadata and all
   selected local hashes. Create one stable idempotency key from that finalized
   metadata and exact file manifest; reuse it only for that same request
   within the API's documented 24-hour idempotency window.
2. `POST /api/v1/submissions` with the completed metadata and
   `Idempotency-Key`, using the scoped bearer credential. Save the actual draft
   ID and response in the ignored local receipt. The first paper's draft has
   been committed; no other paper has a draft yet.
3. For each of the two accepted source files, `POST /api/v1/submissions/{draftId}/uploads`
   with `kind`, `filename`, `relativePath`, `contentType` and `sizeBytes` from
   the queue. Send the bytes using the returned upload's numbered-part
   endpoint, then call its `/complete` endpoint. Preserve the actual upload
   receipts and recheck the staged file list.
4. Use the documented private `/preview/check` and `/preview` workflow before
   commit, coordinated so it does not violate the repository's serial
   typesetting policy. Check the selected entry point, compiler and resulting
   PDF. The first paper's preview passed these checks; the remaining previews
   will run only after their dates and the shared execution lane permit it.
5. `POST /api/v1/submissions/{draftId}/commit` once ready. A successful `202`
   queues asynchronous screening; it is not publication. Poll
   `GET /api/v1/submissions/{draftId}` and retain actual status/identifier
   receipts. Never replace a timed-out response with a blind duplicate commit.
6. Start the next paper on the next eligible UTC day. Unpublished deposits can
   be canceled through the documented lifecycle endpoint; published versions
   are retained as scholarly records under Hexagon's policy.

## Installed calendar job

The installed service is
`gui/501/org.coxeterpolytopes.hexagon-submissions-20261123`, with its plist in
`~/Library/LaunchAgents/`. It checks at minutes 5 and 35 of each local hour
and on load. Eligibility is evaluated in UTC, in the order above, with at
most one accepted submission per UTC day. The first actual calendar
invocation exited 0, made no network request, and deferred the second paper
until its eligibility date. The job disables its own label after all three
accepted-submission receipts have been validated. These are earliest dates,
not guaranteed publication dates or execution times: the computer must be
available, and an occupied scientific lane or unavailable resources defer
the attempt.

The [client](client.py) runs through the publication-I/O watchdog. The
[independent source review](../../audit/hexagon_submission_scheduled_queue_independent_source_review_20261123.json)
binds its exact code, watchdog, calendar configuration and
[approved payload snapshot](approved_payloads.json). The scheduled worker
also verifies the release manifest, source/PDF hashes, metadata and the
first actual accepted receipt before any queued remote action. A changed
paper or configuration stops the queue; it does not authorize a revised
submission. The reviewed configuration still records its preparation state
because its bytes are pinned. The separate installation receipt records the
subsequent successful installation.

Each scheduled submission holds the same private mutation lock as manual
operations and the repository's shared scientific/typesetting lock for its
whole pipeline. Busy locks defer before a draft is created. Exact current
upload IDs, sent-byte hashes, preview digest, successful terminal response,
bibliography and PDF content are checked before the single commit attempt.
An unknown mutation outcome stops for reconciliation rather than a blind
retry. Credentials stay out of arguments and logs; redirects are refused.
Receipts and worker logs remain in ignored owner-only `.local` files.

If `unresolved-preview.local` exists in this directory, the coordinator
must also hold new scientific work until the actual remote job is
reconciled, even if process termination released the execution lock. The
marker is cleared only after verified terminal evidence has been saved.
Neither that marker nor `unresolved-commit.local` was present at installation
verification. The default `prepare` action makes no network request. Remote
previews run on Hexagon; this task does not run local TeX.
