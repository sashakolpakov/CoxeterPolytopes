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

A [follow-up source cleanup](obsolete_entry_points.json) removed six superseded
search/preparation entry points after their consumers were reviewed. The current
mathematical helpers remain; obsolete production bodies in retained helper
modules were removed and their execution functions reject old solve routes.

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

The archived `global_edge_closure_sdpa07/native.result.gz`, the archived original
candidate used for new scaling, current original primal/dual, accepted exact
witness/certificate, final PDF and six accepted TeX editions remain preserved.
Newly generated cone-reduction/pipeline work was outside the initial cleanup scope.
No removed path is among the 936 bulk-manifest entries. This cleanup does not
restore or verify those archived historical inputs.

## Superseded representations, 26 November

The [follow-up receipt](obsolete_representations_followup_20261126.json) records
six removed untracked caches from `global_edge_closure_prepare02` and
`global_edge_closure_conditioned01`, totaling 866,316,600 logical bytes. Their
exact lengths and SHA-256 hashes were recorded before deletion, after checking
disjointness from the current `ready03` bound artifacts and bulk-manifest entries.
The new preparation, proof replay and post-scaling audit had completed; current
bound-artifact paths and the explicitly retained old receipts remain present.

These obsolete representations were deliberately discarded without new archives.
The retained native07 archive does not include their old reduced model or inverse
coordinate maps, so those caches can no longer supply that reconstruction. The
failed candidate03 remains in original coordinates and remains checkable against
the frozen original input. This follow-up ran no mathematical check or solver.

After the corrected pipeline and native run finished, a separate
[owned-cache cleanup](../global_constraint_pipeline_local_cleanup_20261126.json)
removed 18 further files, totaling 2,368,162,349 logical bytes. These were
intermediate compiled models, completed partial-replay caches, byte-verified raw
duplicates of archived checkpoints and reproducible restart files. All current
launch-receipt dependencies, reduction proofs and the terminal native state
remain present. Intermediate models can be reconstructed by replaying the
retained proofs; the terminal checkpoints have complete byte-verified archives.
