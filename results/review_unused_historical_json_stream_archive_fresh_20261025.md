# Independent review: one unused historical JSON streamed archive

Reviewer audit_finite57_algorithms,25October2026. Root authored the utility.
I read the entire first and final sources, requested archive-identity and
final-manifest revalidation, then checked those repairs. No utility, scientific
worker, import, enumeration or test was executed by this reviewer.

Final source: `tools/archive_unused_historical_json_20261025.py`,
SHA-256`06eba30d757f8c90245de40c04f9b2a619ac9347aa8fd67e9d38291156832e35`.
The predecessor f3e9cb01 is preserved by root. The earlier transparent-ditto
attempt is a separate failed storage operation, not a successful compression
or scientific result. Root reports it stopped safely at the256MiB RSS guard;
this source review does not independently execute or certify that attempt.

The new utility has exactly one allowlisted historical result, fixed raw size
1283771506bytes and SHA-256
`4328a9600faec034bd59069e609a9d9323aceea4b9cc00383472c85f883c9437`.
It requires the raw path to be a regular nonsymlink file, currently untracked
by Git, and absent from the exact pinned110-input terminal manifest8a0973d7.
It creates a previously absent ignored owned archive directory and an exclusive
gzip target; it never expands its fixed FILES list dynamically.

Compression and decompression are streamed in1MiB blocks. The original stream
must have the fixed count/hash, and a complete gzip round trip must reproduce
that count/hash. The compressed stream itself receives a SHA-256. Archive
identity is recorded before the round trip and checked after the round trip,
after compressed hashing and immediately before raw removal. Original identity
is checked after the reads and immediately before removal. The exact110-input
manifest is re-read with its byte cap after the forced pre-removal resource
probe and must still equal the initial pinned bytes.

The archive manifest, containing original and archive paths, raw identity,
compressed identity and completed round-trip assertion, is written before any
raw unlink. Thus no raw content is removed until its authenticated restorable
archive exists. Failure before that point retains the original; an owned
partial archive may remain for explicit cleanup. Interruption or final probe/
receipt failure after removal may leave no external receipt, but the completed
archive and its earlier manifest remain. This is a nonhostile-filesystem byte
preservation argument, not a power-loss durability claim.

The archived object is deliberately **not** claimed to remain readable at the
original raw path. The receipt documents restoration by streamed gzip decoding
back to that path, followed by fixed SHA verification. No scientific JSON is
decoded and no numerical verification occurs. No reference licence, public
upload, publication packing or regeneration claim is introduced.

The existing exclusive proof lock is mandatory. Storage probes run before any
mutation, every0.2s while streams advance, and at forced pre-removal/final gates.
Every1MiB block checks owner and the59-second deadline. Disk20GiB, RAM3GiB,
sampled process-group RSS256MiB and bounded3s OS probes remain enforced.
The only external processes are the bounded ps/vm_stat probes and mapped-cache-
limited read-only Git inspection. Monitoring refusal stops the task; there is
no lowered mathematical admission gate or uncontrolled scientific subprocess.

At this stated scope I found no remaining material defect in the final source.
This is independent source clearance for root's fixed single-object storage
operation, not proof that it ran, reclaimed disk space or changed any
mathematical verification status. Root owns execution and its receipt.
