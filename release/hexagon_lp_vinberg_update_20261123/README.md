# Correction of the existing LP submission

The user's explicit instruction to update the Hexagon submission was
executed against the existing work **2611.00962**. The corrected package
was accepted with HTTP 202 at 03:32:58 UTC on 24 November 2026 as
**2611.00962v2**. Its immediate status was **processing**; publication is
not yet confirmed. No duplicate paper was created.

The [authorization](authorization.json) binds the exact canonical source,
the matching 36-page PDF, the new build/package receipts and the unchanged
metadata. Removing the single remark labeled
`fv:vinberg-weight-obstruction` restores the original v1 source's exact
SHA256 and byte count. That v1 source, bibliography, compiled bibliography
and accepted receipt are retained in the ignored owner-only
`payload_history.local/` directory. The old queue's `approved_payloads.json`
and original accepted receipt were not rewritten.

The official API schema was re-read and remained version 1.13.0-beta,
SHA256 `46e761cf8bfafc54589139426026ab7875bae5f1a4f0edb48243cddda8c092f3`.
The documented workflow was used: pause the latest never-published package
with `/corrections`, open its editable correction draft with `/revise`,
stage the updated sources, preview them, then commit. The server assigned
v2 under the same identifier. Other papers were not updated.

The revision endpoint cloned the old staged files into a distinct editable
draft. The [first adapter stop](initial_stage_reconciliation_20261123.json)
was retained before replacing them. The cloned inventory was checked,
then only those editable staged copies were removed; both deletion calls
returned HTTP 200. The old committed package and its history remain intact.
The same two allowed files, `main.tex` and `references.bib`, were then
uploaded with the new exact byte bindings. The metadata, including title,
author, abstract and personal declarations, was preserved.

The remote preview succeeded. All 36 pages' normalized text and the
bibliography exactly match the locally rebuilt edition. BibTeX ran and the
final pass has no undefined references. The
[actual completion record](completion.json) includes those checks, the
new source/PDF hashes, the accepted version and its immediate remote state.
Raw API receipts and preview artifacts remain in ignored owner-only files.

The successful correction pipeline used the unchanged publication-I/O
watchdog and the shared scientific/typesetting lock: 73.936 seconds,
91,226,112 bytes peak process-group RSS, exit code zero. Initial free-disk
and available-memory fractions were 35.124743% and 41.67544%, with the
unchanged percentage resource policy. The parent coordinator kept its
scientific lane free. No local TeX, Lean or mathematical checker was run.
Both global unresolved-operation markers were clear at completion and the
lane was released.

The retained bounds are 57/17/15. The bibliography and compiled bibliography
are byte-identical to v1. No bulk research files, credentials or third-party
reference PDFs were uploaded. The original submission and source/configuration
records remain historical evidence; this correction has separate authorization,
source bindings and receipts.
