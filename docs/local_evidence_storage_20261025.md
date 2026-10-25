# Local evidence storage

A stale168MiB temporary prefix was removed after every byte was compared
with the preserved historical result. The full archive was also verified
in compressed and decompressed form before removal. The
[preservation record](../results/temporary_prefix_preservation_verified_20261025.json)
gives the exact restoration procedure; the
[removal receipt](../results/temporary_prefix_removal_20261025.json)
records1.03seconds and28.4MiB sampled group RSS. No scientific input,
mathematical check or proof admission rule changed.

Eight further parked H9 enumeration outputs were compressed with complete
byte round trips, reclaiming **310.6 MiB**. The
[receipt](../results/parked_enumeration_json_streamed_archive_20261025.json)
records each original path, content hash and archive. Eight finite57
input registers were checked for nonintersection before removal. This took
2.386seconds with31.5MiB sampled group RSS; no mathematical worker ran.
Restore the recorded bytes and verify their hash before using an original
path. The ignored storage is
`formal/.proof-batches/parked_enumeration_json_archive_20261025/`.

A second operation archived115 further outputs from parked H9/H20/H3 and
√5 searches, reclaiming **2.29 GiB**. Its
[execution receipt](../results/parked_search_json_streamed_archive_20261025.json)
binds every original and compressed object. All complete gzip round trips
passed before any original was removed. It took11.706seconds, with sampled
peak group RSS33,636,352bytes. These files are ignored local storage under
`formal/.proof-batches/parked_search_json_archive_20261025/`; the same exact
restoration requirement below applies. The current110-file terminal
manifest was unchanged. This operation was not a mathematical check.

The 25 October storage operation compressed 43 old √5 search outputs and
reclaimed **5.39 GiB**. It used streamed gzip, checked each original SHA256,
then verified every complete decompression before removing any raw copy.
The [execution receipt](../results/parked_sqrt5_json_streamed_archive_20261025.json)
lists every original path, size, content hash, compressed path and hash.
The operation took23.396seconds and sampled a peak process-group RSS of
33,718,272bytes. Its disk/RAM floors and exclusive lock remained enabled.

The current 110-file finite57 terminal input manifest was checked before
and immediately before removal. None of its files was archived. This does
not mean every historical script can run without restoration: a script
requiring one of the listed raw paths must first receive that exact file.

The gzip files and byte manifest are local, ignored storage under
`formal/.proof-batches/parked_sqrt5_json_archive_20261025/`.
To restore a file, stream its gzip archive to the receipt's original path,
check the uncompressed byte count and SHA256, and only then use it as an
input. The archival helper is a fixed, single-use maintenance operation,
not a transparent data loader or a mathematical verifier.

Earlier byte-preserving maintenance has separate receipts for
[one historical JSON output](../results/historical_json_streamed_archive_20261025.json)
and [unused Lean cache objects](../results/lean_generated_cache_byte_preserving_cleanup_20261025.json).
Sources, final PDFs, certificates, execution records and history remain
preserved. These operations do not rerun any successful proof or validate
the mathematical contents of an archived file.

No archive was uploaded. The separate public bulk manifest still binds
936 scientific files in33 archives, and its upload remains subject to the
existing exact-payload authorization requirement. Local byte restoration
and mathematical verification are distinct operations.
