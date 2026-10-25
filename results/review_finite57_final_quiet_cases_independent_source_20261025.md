# Independent source review: one-case final quiet receivers

Reviewed source only on 2026-10-25. No worker was imported or executed; no
syntax check, test, enumeration, build or Git operation was performed. The
new worker and unchanged parent were read in full relevant sections. The
parent helper-loading chain was inspected as text. No blocking source defect
was found at the parent's explicit mathematical premises.

## Frozen review target

- `scripts/check_finite57_final_quiet_cases_20261025.py`: 152 lines, 8,532 bytes;
  SHA256 `9e4e56ce9b290e0b660f8b3c7c87442a986c9d01a1da93ed40bd6b8f36d335b4`.
- Unchanged parent
  `scripts/check_finite57_current_universal_receivers_20261025.py`: SHA256
  `29ddbcf65cae2a7c7ff7703a1b4a81bad3ca1fdb0ade4bdbffda26fb9c8a1631`.
  Comparison is specifically with `final_quiet`, lines 498--548, and its
  `quiet_physical`, `current`, `role`, `domain`, `assumed`, `application` and
  `context` helpers. This review does not extend their accepted premises.

## Early partition and copied body

`CASE_INDICES` is exactly `(3+6*i+j)` for `0 <= i < 8`, `0 <= j < 3`:
24 distinct original indices in the 48-row source. Batch names use these
original indices, rather than ordinal positions in the 24-row subset.
`quiet_case` checks the full original length and the exact, duplicate-free
repair-index domain, then selects `p = by[index]` and `q = original[index]`.
Selection precedes `quiet_physical`, `E.matrix`, `current`, `role` and
classification. Before selection, the two `c.field` calls load pinned JSON
projections and `domain` checks list/key identities; they do not evaluate the
other 23 physical cases. The worker never calls the parent's `final_quiet`.

The numerical and application body equals the selected iteration of the
frozen parent, allowing whitespace, qualified helper names and the local
`out` placement. In particular it retains:

- Original ordered marks, left/right side and exact
  `isolated4`/`mate4`/`center3` profile associations; complete reconstructed
  Gram equality and the original white set; the same two distinct actual
  quiet and paying corners.
- All eight ordered triples in `product(('A','D'), repeat=3)`, with three
  separately labelled whole rank-15 copies. The quiet receiver is `1/2`
  with complete rank-58 forest and no ignored unsealed component.
- `isolated4` paying value `1`, and `center3` paying value `1/2`, with the
  same complete rank-58 forest checks and the same three witness copies at
  both corners. The whole-star E0/full-Z/both-entire-columns seal remains
  an explicit assumption with unchanged text and bridge source.
- `mate4` paying intrinsic value `0`, ideal whole `B13`, exact images
  `[[10,11],[10,12]]`, the named independent `quiet-b13` provider and exact
  `[index,*kinds]` consumer key. Its total remains `None`; no zero is
  credited as a unit.
- The same theorem, repair-selector position, original-source selector,
  original-source context hash, reconstructed physical context hash,
  ordered copy identifiers, receiver records and structured assumptions.

## Coverage, termination and replay

The case checks exactly eight distinct application keys, equal to
`{(index,*kinds): kinds in product(('A','D'), repeat=3)}`. It claims
`case_complete_at_declared_premises=True` and
`quiet24_complete_at_declared_premises=False`. The 24 disjoint first-key
coordinates therefore specify the original 192-key product. Every case
advertises the same ordered original-24 domain. An eventual union must
check all 24 actual case receipts and their full keys; this worker emits no
aggregate receipt or assertion that the other cases ran.

The eight-iteration loop and shared engine's unchanged fuel cap of
1,000,000 bound numerical traversal. The per-case worker makes 16 `current`
calls, rather than the parent's 384 across all cases. Each case reads the
same two source JSON projections and the same preparation inputs; there
are no branch-specific source reads. Consequently unchanged-source case
receipts have equal `consumed_inputs`, including the added parent pin.
Replay retains original selectors and application keys, not a renumbered
or synthetic parent receipt. Runtime fit within the existing limits has
not been established by this source review.

## Loading and refusal guards

The new worker requires an exact allowed batch and `--parent-lock`, rejects
optimization and the optimization environment, pins the parent before
loading it in a private module, and calls only its existing `prepare`.
The private module name prevents the parent's main entry point from
running. Its preparation chain still calls the collector's shared-lock
check before source input loading or numerical work. It does not create a
replacement lock or start a runner.

The parent's loader `76022d95...`, route `43a0c4e...`, collector `930ffff5...`,
engine `39dd55d2...`, literal library `f6eba7af...`, source ledger
`68359555...`, extra manifest `a79e1b3d...` and base input manifest
`8a0973d7...` match their existing pins. The new worker adds the exact
parent source to consumed pins, rejects a conflicting pin, checks every
consumed digest after computation, and rechecks its own bytes, parent,
ledger, loader, extra manifest, route, collector and base manifest before
emitting output. The unchanged helper caps source reads and shared fuel.
The output cap is 500 KiB including its final newline, with refusal rather
than truncation.

All independent-S8, protected74, terminal-bridge and global-finite57 flags
remain false, as does independent quiet-B13 runtime verification. The
existing external launch contract remains one monitored process, 59
seconds, 1 GiB and one thread under the shared lock, with unchanged
40/6 GiB admission and 8/3 GiB abort floors. Those runner requirements
are not replaced by the worker's lock or fuel checks.

The final consumer, controls, aggregate receipt validation and resource
admission are outside this review. No source modification is requested.
