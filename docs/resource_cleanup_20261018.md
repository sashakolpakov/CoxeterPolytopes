# Disk cleanup, 18 October 2026

Cleanup preserved every certificate, report, source file, and session transcript.
Static large JSON files were transparently compressed by APFS; their paths and
decompressed bytes remain unchanged. No active Codex database was touched.

Three bounded batches handled 34 files: 33 workspace JSON files and one inactive
Codex session older than seven days. They reclaimed **4.37 GiB of allocated file
blocks**. Streaming SHA-256 comparisons passed for every source/copy pair and
again for every final replacement. File size, modification time, permissions,
and ownership were checked. Free space after cleanup was approximately **75 GiB**.
An earlier, unrelated increase of approximately 50 GiB was observed and is not
attributed to these batches. Detailed measurements are in
[resource_cleanup_20261018.json](resource_cleanup_20261018.json).

The reusable [cleanup script](../scripts/reclaim_generated_space_20261018.py)
processes one file at a time with 1 MiB hash buffers, a maximum 512 MiB source
file, batch limits, and a minimum 4 GiB free-space reserve. It skips symlinks,
hardlinks, already compressed files, immutable files, recently changed files,
and files reported open by `lsof`. Temporary copies are removed on ordinary
failure. It replaces a source only after byte and metadata checks and only
when physical allocation decreases. The APFS compression flag makes normal
reads return the original data; copying to a filesystem without compression
may require the original logical space.

The workspace threshold is ten minutes of inactivity; session history requires
seven days. Session compression needs write access outside the workspace. No
additional history batch was run after sufficient free space was established.
The pilot and batch logs record actual free-space deltas separately from
per-file allocation savings because unrelated system writes can affect free
space concurrently.

No proof or search task ran as part of cleanup. Generated TeX intermediates and
Python cache directories were separately removed by the root agent while
retaining PDFs and research inputs. Subsequent proof work must use the bounded
batch runner and remove its temporary files after each batch.
