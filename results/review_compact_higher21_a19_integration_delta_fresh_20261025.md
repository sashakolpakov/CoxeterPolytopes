# Higher tail through exact21 and whole A19: integration-delta review

25 October 2026. SOURCE ONLY. No proof/test/import/build job; no edits
to any manuscript or frozen source. This is a focused integration
review following the two complete mathematical reviews, not a second
full reading of unchanged bodies.

Reviewed receipts:

* `results/compact17_higher_tail_exact21_working_integration_20261025.json`,
  SHA256 `53d620b9e1398b74ce64c0e2f21bf6b9010cfd867a022617f88fd7afde00d6d6`.
* `results/compact17_whole_a19_working_integration_20261025.json`,
  SHA256 `723dfed2acb7cd1f4ed5d74361bc1790c17100ec01803b468e3c007545cb5a28`.

The complete non-append hunks and complete append-to-frozen-source
differences were read. The source comparisons were ordinary textual
diffs, not mathematical execution.

Higher-tail patch
`results/compact17_higher_tail_exact21_working_integration_20261025.patch`,
SHA256 `782b37491dc78d90eacdd76cb940ff6f5c61d94274c0f64a9d1508ae00d3ced4`,
maps main 01309ce5 to 4bcef9f1. The added mathematical fragment is
identical to the reviewed 7515d8db source except its source-status
comment and the single beta-domain clarification. This is confirmed
by a full textual diff of the extracted append against the immutable
gzip's restored source.

That clarification correctly restricts beta5>.1, beta6>.012 and
beta7>.0095 to the retained G_i-floor reference when five is not a
weak helper. The separate weak-helper-five interpolation charges
only dotted six, whose actual coefficient exceeds .10 and therefore
.012, and does not charge dotted seven. It explicitly denies an
arbitrary-concave-vector interpretation. This is exactly the domain
distinction reconstructed in the complete outside review
`results/review_compact_higher_tail_through_exact21_fresh_20261025.md`.
No scalar, branch, price, pressure or other proof body changes.

The overview and remaining-obligation hunks now say that dimensions
29 through 21 are supplied, exact20/19 remain to be integrated and
jointly reviewed, and the published compact bound is used only at
dimensions at least 30. They do not claim a complete printed proof
of the higher-tail theorem. The general-only/arithmetic separation
and the final joint dependency/reproducibility/build obligation remain.

A19 patch
`results/compact17_whole_a19_working_integration_20261025.patch`,
SHA256 `68000d62fc114efd924a51d9f030c3d40cc15f629bcdb2fb90387ab04870ddd7`,
maps main 4bcef9f1 to 483453fb. It has ONLY the 904-line append;
the extracted append is byte-for-byte the reviewed 0eb7982 source
apart from its single source-status comment. The actual whole-A19
proof and shared N20 window reference are unchanged. The reviewed
higher-tail source supplies that shared lemma at N20, without using
whole A19 or a larger-rank whole exclusion as its premise.

The receipts bind successive predecessor/successor hashes correctly:
01309ce5 -> 4bcef9f1 -> 483453fb, the latter 688987 bytes and SHA256
`483453fbd0ad6298e48fc217d5c9c209b28709d673f384e75518d1d419dee3c0`.
Preserved predecessor gzip archives and byte-round-trip claims are
recorded by root; this review does not credit their I/O as proof runs.

Verdict: these exact integration deltas are supported at the scopes
of their complete preceding reviews. No new defect found. This
does not review later main edits, exact20/19, whole B19/D19, the
whole compact17 composition, a compiled PDF or a Lean kernel.
