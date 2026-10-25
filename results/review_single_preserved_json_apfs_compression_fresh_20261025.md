# Independent source review: one preserved JSON APFS compression

Reviewer: audit_finite57_algorithms,25October2026. I did not author the
storage utility or its successive repairs. I authored the earlier reviewed
OwnedTree monitor; root's separate full review and the preserved infrastructure
regressions supply its independent review. This note does not convert my
monitor self-review into another-author review.

Final utility: `tools/compress_one_preserved_json_apfs_20261025.py`,8519bytes,
SHA-256`90d1a575cf1275cdf0e4f8f1012e7e55ba83cde9077350fae92835bd86ec7c79`.
Pinned monitor:
`scripts/monitor_owned_proof_tree_20261024.py`,
SHA-256`b07f84fbcc5b35e484696f0515b8e233ea975f13c2da6ff75c0d3b54b99ed5db`.
No utility, import, copy, scientific worker or test was run by this reviewer.
Only source reading and source hash/size metadata were performed.

I read the entire first utility, entire repaired versions and the final bounded
read delta. The first findings were missing copied-file identity bracketing,
source revalidation immediately before replacement, PID-only cleanup versus a
claimed owned-group scope, and a normal module import after hashing. The
intermediate final monitor reread was also unbounded. Root preserved those
predecessors and repaired each issue; they earn no execution credit.

The final utility operates on exactly one fixed historical path,
`results/fresh_sqrt5_rank13_all_selected_nonade_separation_20261017.json`,
of size1283771506bytes. It decodes no scientific JSON. It hashes the original
stream in1MiB blocks, records inode/device/size/time identity, and copies only
to a newly owned same-directory path using the fixed official ditto invocation.
The original remains present throughout copying and verification.

The copied stream must have exactly the same byte count and SHA-256. Its file
identity is recorded before that hash, checked after it and again immediately
before replacement. The original identity is also rechecked after the final
resource probe immediately before replacement. Mode and numeric owner/group
must agree, and allocated blocks must strictly decrease. Only then does one
same-filesystem atomic replacement retain the original path and verified bytes.
A failed pre-replacement copy is removed only at the newly owned temporary path.
The utility makes no power-loss durability or hostile concurrent-filesystem
claim. An interruption or receipt-write failure after replacement may leave no
receipt, with the already verified byte-identical source in place; the emitted
receipt explicitly records this boundary.

The owned monitor is loaded from its exact bounded hash-checked bytes with a
private ModuleType registered in sys.modules, so neither a cache nor a second
import reread determines executed monitor code. Its main is not called. A
ProcessReader probe precedes copying. The ditto subprocess starts a new session;
its observed birth identities and descendant groups are tracked through the
reviewed OwnedTree API. A normal successful copy requires no cleanup signals.
Failure cleanup uses that ownership monitor; if initialization itself fails,
the tracked Popen child is stopped and the original is not replaced. This
retains the monitor's disclosed polling/observed-ownership scope and does not
claim atomic discovery of unseen descendants.

The shared proof lock excludes concurrent proof batches. Storage checks retain
20GiB disk,3GiB RAM, sampled256MiB process-group RSS and a59-second deadline,
with bounded OS-probe timeouts. The source includes its own group plus observed
owned copy groups in RSS. Those storage-only checks neither lower nor satisfy
the scientific40GiB/6GiB admission gate. Failure of a monitoring probe refuses
the operation; the code does not bypass it. Thread-count environments are
restricted for the copy. No publication upload or mathematical process occurs.

I found no remaining source defect material to this explicitly bounded,
nonhostile, storage-only byte-preservation scope. This is source clearance for
the final pinned utility, not a claim that compression ran, reclaimed space,
regenerated any certificate or proved any mathematical statement. Root owns
execution and the resulting receipt.
