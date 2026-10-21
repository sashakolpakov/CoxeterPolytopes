# Committed 1,097-core snapshot and witness replay

**Final/frozen: generation and saved replay pass.** The [checker](../scripts/fresh_sqrt5_nonsingular_attachment_search_20261015_batch1097.py)
and [certificate](fresh_sqrt5_nonsingular_attachment_search_20261015_batch1097.json)
bind the immutable
[full snapshot](../data/enumerations/fresh_sqrt5_nonsingular_attachment_search_20261015_batch1097.json).

The unchanged general producer was stopped by SIGINT while evaluating
core 1097. Its last committed cursor is 1097, so cores 0 through1096 are
complete. The snapshot retains all 1,281 states and 184 pending cores, all 311
literal seed maps and every accepted-child map. No database transaction was
edited and the unfinished attempted core was not promoted to a result.
The saved batch-history field still describes the completed earlier 610
boundary; the atomic next-core cursor records the additional 487 completed
transactions. This distinction is intentional and documented.

The archive's `database` and `database_sha256` are **historical metadata
only**. The checker binds the archived JSON bytes and never reopens that
database. A recursive archival audit must exempt precisely these two
historical fields in this archived JSON from any live-database equality
test, as for the separately frozen batch610 archive. It must still verify
the snapshot bytes, all source hashes, literal maps and replay certificate.

| Checked item | Count |
|---|---:|
| Committed marked cores |1097|
| Complete terminal one-white profile occurrences |13025|
| Literal seed and child maps |2438|
| Exact affine exclusion witnesses |890|
| Literal inherited closed-subcore inclusions |98|
| Accepted child occurrences |2127|
| Pending cores preserved |184|

The checker reconstructs each complete finite-corner/no-short profile
domain. It verifies every stored H4, affine, ground-field, physical-relation,
physical-inertia and conjugate-inertia exclusion, every accepted child and
map, all old-rule stopping payments and full reachability. For a paid
subcore it checks the displayed full-label inclusion and, when sealed,
orthogonality of every other present white root. It does not repeat the
factorial search for that already supplied witness or reproduce the count
of unsuccessful subsets. Those counts play no role in the verified proof.

The committed outcomes are 642 expanded nonsingular cores, 76 sealed
nonsingular rank exclusions, 98 inherited subcore payments, 221 new sealed
old-rule payments, 50 sealed unpaid cases and 10 expanded singular cases.
There are 379 processed singular states in total. An older payment is not
automatically valid for a changed forest-sensitive rule.

This is a proof-witness and complete-domain replay using the inherited
reviewed Fraction field, not an independent arithmetic implementation or
the unchanged producer's full saved replay. It makes no all-size closure,
actual quadrilateral realization or dimension-bound claim. The next
[dimension-specific fork](fresh_sqrt5_nonsingular_rank29_exhaustion_20261015.md)
starts afresh from the same 311 literal roots and records every singular
zero-touch weight obligation, including all cases paid by old rules.
