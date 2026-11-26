# Search cleanup, 25 November 2026

The user requested removal of obsolete and failed searches while the current
global LP/SDP was being repaired. The coordinating agent approved the exact
legacy families and the subsequent list of superseded numerical payloads.

[`result.json`](result.json) records the actual removals. The cleanup deleted
5,443 obsolete untracked search artifacts, notes, scripts and audit records;
14 byte-restorable raw duplicates; and 124 superseded global numerical
payloads, including 16 tracked files. The total is 11,435,280,853 logical
bytes. Logical file lengths are not a claim about physical disk savings.
No scientific job, solver, Lean build or TeX build was launched by cleanup.

The exact candidate/dependency inventory is in [`inventory.json`](inventory.json).
The actual path, byte length and SHA-256 of each deleted legacy file is in
[`removed_files.jsonl`](removed_files.jsonl); the later numerical payloads
are in [`removed_global_payloads.jsonl`](removed_global_payloads.jsonl).
[`dependency_exclusions.json`](dependency_exclusions.json) distinguishes
protected proof inputs from conservative historical dependency exclusions.
Those exclusions do not assert that every old helper default remains needed
by a current runtime checker.

## Historical execution records and restoration

The compact historical result/status/source-binding records remain unchanged.
They describe what actually ran at their pinned source versions. A path in
such a record does not promise that a superseded numerical cache remains
present in the current checkout, or that the failed attempt proved anything.

The 16 deleted **tracked** payloads matched the pushed recovery tag
`pre-search-cleanup-20261125` at commit `70b23827` before deletion. Their
exact paths are in
[`tracked_payload_live_text_references.json`](tracked_payload_live_text_references.json).
No live Markdown/RST document linked directly to those exact payload files.
To restore a historical tracked payload into the worktree, use its recorded
path, for example:

```sh
git restore --source=pre-search-cleanup-20261125 --worktree -- results/global_shared_cusp_sdpa/input.dat-s
```

Run historical procedures with their recorded source version, rather than
interpreting old commands as the current solver workflow. Restoration of
bytes is not a mathematical verification.

The 14 removed recent-run raw duplicates remain restorable from their
retained `.gz` archives. Before deletion, both the current raw and archive
hashes were checked against already-passed complete byte-round-trip receipts.
The exact archive paths and receipts are recorded in `result.json`.
Generated `native_initial.dat` duplicates were also checked against exact
balanced-array extraction from their archived source states.

The recovery tag covers committed history. Obsolete **untracked** generated
search payloads were deliberately deleted under the cleanup instruction;
their recorded hashes do not constitute backup copies. Their deletion
removed no current accepted proof input or bulk-manifest asset.

The latest `global_edge_closure_sdpa07/native.result`, the archived original
candidate used for new scaling, every current frozen model, accepted exact
witness/certificate, final PDF, six accepted TeX editions and all 936 bulk
manifest paths remain preserved. Newly generated cone-reduction/pipeline
work was outside the cleanup scope.
