# Staged publication audit: source review

21 October 2026. **SOURCE PASS at the publication-hygiene boundary.** I read the complete audit and publication watchdog, including all corrections below. I did not author these utilities. I proposed the corrected raw-marker and manifest checks; shared review contributions are disclosed. No utility import, Git command, test, proof worker, scientific job, Lean invocation or document build was executed by this reviewer.

| Source | Bytes | SHA256 |
|---|---:|---|
| `tools/audit_public_tree.py` |7210|`0a25fb8e6981a5e1ed8d77e57a4035a810e719535e4a7055a95c6012c9e71bfd`|
| `tools/run_publication_io.py` |4394|`d683e0577b428e768f860586cca923337d0d0c9981c1c1cd5d199623ceff858d`|

The audit obtains the staged index entries and reads their exact object IDs through one `git cat-file --batch` subprocess. It sends one request, drains that blob and its terminator, then sends the next request; it does not fill both pipes with an unbounded batch. Reads are at most one MiB, except the separately bounded four-MiB staged bulk manifest. Empty/unmerged indexes, invalid blob responses and truncated streams fail. Unsupported modes, known cache/compiler/native-library/temporary paths, blobs reaching100MiB, and possible secret markers prevent a successful hygiene result. Diagnostics give file paths and reasons, not secret contents.

The secret scan now covers raw binary as well as text bytes, carrying256bytes across block boundaries. The private-key/GitHub/AWS marker patterns are heuristics, not an exhaustive credential detector. The PGP header is split across adjacent source literals: the compiled expression still detects the full header, while the audit source no longer flags itself. I verified that this split is the sole delta from the previously fully read corrected audit.

The staged bulk manifest must have the supported version, nonempty parts/files, distinct safe archive names, unique canonical relative `data`/`results` paths, exact nonnegative integer lengths,64-hex SHAs and a matching original-byte total. Listed release-only bulk paths must be absent from the index. Thus the test uses staged metadata rather than a possibly different working manifest. Malformed inputs fail rather than produce a successful audit.

All six selected PDF paths must be present. **Presence is not PDF identity, source/PDF consistency, or a mathematical proof check.** Those remain the separate frozen-artifact receipts and written reviews. Likewise this scan cannot recognize every arbitrarily renamed binary or prove the absence of every possible secret.

The watchdog successor differs from the reviewed publication monitor only by adding `tools/audit_public_tree.py` to its explicit allowlist. It retains the owned-process-group, sampled256MiB RSS/publication resource policy; it does not relax the separate40/6GiB proof admission gates or authorize a proof runner. Root reports a separate actual staged-audit success; its execution/output binding belongs in the separate runtime receipt, not this source review. No further source blocker was found. The user's requirement to finish the pending Lean compilations and actual axiom audits before commit/push remains in force.
