# Full outside-author source review: four-profile D19 census

21 October2026. **PASS for the current source; no remaining correction.**
This review binds:

* Worker `scripts/check_d19_nonisolated_quadratic_fork_time_compatibility_20261021.py`,
  SHA-256 `f1b3f5097d60cc8dc3dc2e498409d513cd8dc1adb0306c527be2f1b279b34e89`
  (9697B/203lines).
* Mathematical specification
  `results/arithmetic_d19_remaining_dotted_fork_physical_bounds_and_census_20261021.md`,
  SHA-256 `80d67e3b48b33534f7e716e2011323123d4b05a4f1aa6fe9fd770919bb9b930f`.

The review reads the entire new worker and specification and the imported
base source. No checker was imported, executed, fixture-tested or sampled.
Only source reads and stdlib metadata/hash/archive comparison ran. No source
was edited by the reviewer. The reviewer authored parts of the inherited
59e/d8d inventory and supplied some physical-core ideas; this shared
mathematical input is disclosed. The new worker implementation is by a
different agent and is independently reviewed here.

## Resolved pre-run blocker and exact bindings

The first frozen worker befca0ee hardcoded the predecessor specification
hash86d68adc and required the CLI hash to equal it. It could not accept
the corrected current80d specification. This was reported BEFORE execution.
The author preserved those exact9697 bytes at
`results/source_archive/check_d19_nonisolated_quadratic_fork_before_proof_scope_pin_fix_20261021.py.gz`,
gzip SHA-256 `4c64d3c5954da222b00de87cb3566365900913697d826a376e9d275d0d2b0154`.
I independently decompressed and hashed it: the round-trip is exactly
`befca0eebc1b7d2841d5a3c11bc49b82b66494ec2b2e4e8498b130c0b6d73129`.
The ONLY successor difference is the64-character PINS[PROOF] value,
86d68adc to80d67e3b. All formulas, loops and output logic are unchanged.

All four direct input pins match: corrected specification80d67e3b,
base arithmetic/prefix worker cac3bc1d, actual quadratic-field theorem
790b70eb and complete weak-prefix inventory d8d80322. The worker's own
current source hash also matches. These are direct source bindings,
not a new replay or recursive hash validation of every inherited theorem.

## Complete four-family mathematical domain

The specification supplies ARITHMETIC .61 prefix floors and a genuine
actual feasible-core proof of gamma<5. Quadraticity is needed for the
subsequent trace/norm census, not for that physical bound. I checked the
far-core cross term and norm expression, the decreasing-in-w endpoint
comparison, its three rational seals, and the flat near-helper3 core.
The corrections to the former nonarithmetic scope claim are respected.

The code implements every integral rectangle, with inclusive endpoints:

| Family | Trace | Norm | Rows | Exact physical and raw-time tests |
|---|---|---|---:|---|
| FAR3 | 2..10 | -5..3 | 81 | 2<t<10; H=(19t^2+4t-4)/4 |
| NEAR3 | 1..9 | -13..-1 | 117 | 2<t<10; H=(19t^2+34t+11)/4 |
| FAR_DPOS3 | 2..9 | -4..-1 | 32 | 2<t<10; H=t(19t+8)/4 |
| FAR4 | 3..14 | -5..-1 | 60 | t>0,8<t^2<200; H=t(19t+8)/4 |

The total is290. Every row passes through the unchanged exact
positive-root routine; nonpositive and square discriminants are rejected.
This is correct for irreducible ACTUAL quadratic parameters. Every actual
physical parameter is the larger root, because its other conjugate lies
in the strict small interval dictated by H'<0. FAR4's positivity is
tested BEFORE using its squared interval; no negative-root ambiguity is
introduced. All families then require H>0 and H'<0 by exact field signs.
The loose specification strips are used only to prove the integer bounds;
the worker uses the stronger exact H test and therefore loses no actual
anchor. FAR3 norm0 is harmless overlisting and is rejected as reducible.

The raw parameters and times match the specification, including NORM4
and t=2sqrt2*gamma for FAR4. No square-root extension, unit-time substitution
or helper3 root-norm assumption is silently imposed on that family.

## Every prefix pair and the complete field-square test

The only mathematical functions imported from pinned cac3bc1d are its
exact quadratic root/arithmetic routines, COMPLETE prefix-candidate list
and full K-square test. The new script never calls the old isolated
anchor generator, old census, regressions or main, and never opens its
old JSON output. Normal import executes definitions only; the base's
`if __name__ == '__main__'` guard prevents its old run.

Every one of the inherited257 necessary prefix rows is retained. Their
canonical squarefree fields and positive radicals identify the physical
embeddings. The code tests EVERY same-field anchor/prefix pair, with no
early break on a successful witness. Physically and conjugately positive
U/H is checked, then the full exact square routine is called. That routine
includes both norm signs and rational/pure-radical cases; it is not a
norm-only rejection heuristic. Every returned witness is multiplied back
and compared exactly with U/H before recording a match.

Every compatible typed prefix row, index, parameter, time, ratio and square
witness is saved. Per-anchor possible-weak-index sets are correctly deduped
from ALL matches; the global union is their union, without losing the
full underlying match records. Every zero-match anchor is also saved.
Family-domain and candidate counters, total square-test count and reason
counts have the stated meanings. Duplicate time models in distinct role
families remain distinct typed anchors, correctly.

## Termination, input integrity, output and execution boundary

All loop domains are fixed finite sets: at most290 anchors and257 prefixes,
hence at most74530 pair visits before same-field filtering. The imported
root factorization/sign/square tests use exact integers and Fractions and
terminate; the discriminants are bounded by these small domains. There
are no symbolic packages, floating tolerances, root scans, large JSON
loads, subprocesses, parallel workers or unbounded external searches.
This proves finite termination and bounded structure, not measured RSS
or wall time.

All four direct inputs are stream-hashed BEFORE the normal base import.
The imported module's resolved path and source hash are checked immediately,
and the worker plus all four direct inputs are checked again after the
census. The current proof path/hash is enforced by the required CLI.
The JSON destination must resolve directly under results, have suffix.json,
not already exist and be created in exclusive mode. Existing evidence is
therefore preserved. The complete certificate is written there; stdout
omits the potentially large survivor/zero-match lists and prints bounded
counters. The status and scope correctly describe necessary local-data
compatibility, without assuming an empty result or a geometric realization.

The script itself does not replace the external start gates, process-tree
watchdog or shared lock: any actual run still requires the unchanged
root-owned bounded runner and fresh resource admission. This PASS is the
requested BEFORE-RUN source verdict only. There is no actual runtime,
output count, empty-list, branch-exclusion or dimension conclusion in it.
