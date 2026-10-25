# Independent review of the two-object Lean cache archive utility

25 October 2026. **Source clearance accepted at the exact fixed-object
and declared unused-cache boundary.** Root authored the utility and owns
execution. This reviewer read its entire source and did not execute it,
query resources/processes, delete files or import any worker.

Source: `tools/archive_unused_lean_cache_20261025.py`, SHA-256
`cf505749268b706dc8cc59eb26ceab8ae074d07a95ee15eaec69bcc1c68f8c53`.

The operation is limited to TWO fixed, untracked generated cache objects:

|Object|Bytes|SHA-256|
|---|---:|---|
|Generated/Ideal59.olean|107735496|`4f941c500c527933c6de005af23c373718f2f830dcf991e5c24fed8e72ff6588`|
|Generated/Finite59.olean|63746968|`ed42354fb5332a5fb470fa5a21529971311edbc7fed1bfb06b66f9627153b500`|

These live only below the fixed `formal/.lake/build/lib/lean/CoxeterProofs/`
cache prefix. The script requires regular non-symlink files of the exact
sizes and Git's untracked result, using small Git pack/delta caches.
The operational premise that no pending selected job needs these objects
is root's separately declared dependency assessment; this script does not
discover or prove it. Their scientific sources and execution history are
not targets. Restoration does not require recompiling an old module.

The utility takes the existing exclusive proof lock nonblockingly and
requires new archive and receipt paths. It streams in1MiB blocks to
fixed gzip filenames in the ignored proof-cache archive. It computes
the ORIGINAL fixed hash while copying, closes the stream, then separately
decompresses and hashes the entire restored byte stream. Count and digest
must match exactly. The original stat identity is checked after each
round trip and for BOTH objects again before EITHER unlink. Thus a copy,
identity or round-trip failure cannot silently authorize deleting an
unverified original. If interruption occurs between the two unlinks,
both verified archives already exist. No atomic power-loss durability
or hostile concurrent-filesystem guarantee is inferred.

The source performs no upload, publication packing, mathematical import,
Lean/TeX job or signal to any process. Its Git checks and ps/vm_stat calls
are read-only. It preserves enough fixed names/hashes to restore each
cache by streaming decompression and verification. The receipt clearly
distinguishes byte preservation from proof verification.

The dedicated STORAGE guard requires20GiB disk and3GiB available RAM and
samples its own process-group RSS below256MiB, before touching inputs
and during every block/verification phase. Monitoring failure refuses
the operation. This is not a change to the40GiB/6GiB proof admission
policy, a new scientific runner, or an enforced heap allocation cap.

Limit precision:59 seconds is checked at guard entry. Each intervening
vm_stat/ps subprocess has a3-second timeout, so the utility alone does
not guarantee an independent kill at exactly59 seconds. An owner
watchdog supplies a hard wall if required; otherwise this is a sampled
maintenance deadline. It cannot start a proof below its original gates.

Actual archive hashes, reclaimed bytes, resource measurements, cleanup
receipt and post-operation proof admission remain unexecuted here.
No storage result or scientific PASS is credited by this source review.
