# Current D16/D13 output audit: source preparation and coverage contract

19 October 2026. **SOURCE PREPARED, NOT EXECUTED.** The only executed
operations in this lane were source/JSON reading, streamed input hashes,
filesystem metadata, source writing and an AST syntax parse. No audit body,
matrix calculation program, scientific producer, census, solver, Lean or
build ran. The coordinator owns source review and any guarded execution.

Source: `scripts/check_d16_d13_physical_output_retention_20261019.py`,
394 lines,21,039 bytes, SHA256
`36c3b674e296faf6cea5da48d2c1caaf586ec503763377d27089247bf4f15d70`.
It imports only the reviewed standard-library saved-interface reader,
not a mathematical producer. That reader is pinned at SHA256
`328c473c3fa69e547050793485f51cf399dcd667770d642431de25da5f50c908`.

## Read-only input plan before execution

The complete exact path/byte/hash list is
`results/d16_d13_physical_output_audit_inputs_20261019.json`, SHA256
`6b6af762943a9279f747fcceca67de75522d6864201e067aba22630b619cab78`.
It was created from the reader's literal source inventory and the
22-packet registry by AST inspection, without importing an audit or
producer. Each input SHA was streamed in64KiB blocks.

There are54 distinct JSON files and106,205,768 total bytes. Each file's
exact name, bytes, hash and read mode is listed in that small plan.
Only four exceed the4MiB whole-load cap:

|Input|Bytes|Read method|
|---|---:|---|
|`fresh_distance8_high_overlap_20261017.json`|20,858,417|stream to the existing unique1142/common-core markers only|
|`fresh_d17_a9_a2_context_zero_ambient58_frontier_20261018.json`|18,099,126|stream declared `local_proof` section|
|`fresh_order16_d_quarter_epsilon_ideal57_scout_20261018.json`|10,610,727|stream declared `local` section|
|`fresh_order17_d13_shortleaf_seven_thirtyseconds_scout_20261018.json`|10,124,067|stream declared `local` section|

Every other JSON is loaded alone below4MiB and released after its
collection is traversed. No giant ambient57 JSON is read. The54-file
inventory includes the741KB selected D16 packet and168KB order17
partial packet used for local models. Input bytes count distinct JSON
sizes, not repeated streaming passes. Source/report pins are additional
small files. The worker checks equality of its dynamic reader inventory
with the pinned plan, then rechecks every size and hash.

## Recognizers and fail-closed acceptance

The D16 recognizer reconstructs the ENTIRE marked connected16-root
component, its unique fork, its two short leaves and its long endpoint.
The marks must be precisely long-end/short-leaf. The D13 recognizer
requires TWO different entire marked components, an A4 path with an
endpoint mark and a D13 fork with a short-leaf mark. Each target corner
must be positive definite.

Recognition is invariant under positive root rescaling: every connected
edge is tested by g_ij<0 and4g_ij^2=g_ii*g_jj, rather than a literal-1
test. The accepted complete core models, however, use the archived
norm-two normalization and exact full Gram, white set, marks and corner
root set. A differently scaled copy is still recognized as a target but
does not silently pass the exact model check; it becomes unknown and
stops the audit. No subcore is treated as a whole actual component just
because its principal Gram is positive definite.

For D16 the worker accepts only:

* the complete selected1644/1664/1931/2324 radical geometries, with the
  exact displayed A1 and their inherited actual whole-A1 proof; or
* the complete18-root order16 geometries at the four literal corners
 663(1,13),663(3,11),676(3,12),829(2,13), now covered by the independent
 short2/rank-obstruction proof.

It freshly checks all four selected radical relations, their complete
saved allowed-contact entries and shortcut bounds. For the order16 core
it computes an exact rational inverse and verifies its product with the
Gram, its contact block[[4,-4],[-4,4]], the positive affine kernel,
the14-white principals and all nine original actual black sets. Unknown
target geometries, even with an unproved displayed singleton, stop.

For D13 the worker independently constructs all four complete19-root
Gram models in a rational coordinate space:13 D coordinates, five A4
coordinates and a time coordinate with metric coefficient-1/5. Its two
omitted walls have time coordinates5/2 and1, representing physical
coefficients sqrt5/2 and1/sqrt5. The full19-by19 Gram is compared with
each1974/2367/2504/2511 saved core. It checks the positive affine kernel,
radical, white/edge principals, exact marked coupon role, and the rational
capacity products-5/4,3/2,4/5,8/5. The actual-facet projection and capacity
two conclusion remain the separately written geometry.

The worker also reads the saved distal `local` section and checks its
complete140/71/488 domain, four changed intrinsic pairs, all four selected
ideal marks/payments, all four finite coupon marks and the historical
7/32/25/32 split. This links the modeled finite geometry to its actual
paired ideal output rather than matching only the finite corner.

Every D13 target is classified separately as one of:

* a complete matched distal geometry: its whole-face proof can be
  replaced by the parametric identity(1/2-alpha)+(1/2+alpha)=1;
* an explicitly unchanged half-only request; or
* an unknown, which stops the audit.

Per-corner scalar fields `payment`, `weight`, `lower`, `finite_half`,
`new`, `new_weight`, `new_payment`, `current_weight`, `new_lower` are
recorded. Aggregate `total` and unrelated connected `base` fields are
not misread as demands on the marked disconnected angle. A matched
complete distal geometry permits a new parametric proof up to a unit,
rather than freezing the historical25/32 number. Other requests must
be at most1/2. Where no numeric request is stored, acceptance as
half-only also requires the explicit historical/disjoint-role source
registry recorded by the worker. That missing-amount source semantics
is a written premise for root review, not a claim proved from hashes.
An unknown post-distal/marginal source is not automatically half-only.

## Complete saved physical collections

Both modes cover the same explicit physical interface:

* all22 original later packets, with exactly564 reconstructed corners;
* all216 current literal configurations from the seven original
  collections, including the four order16 target corners;
* all19 newer selected physical sections, including rank15 sparse
  `selected_proofs`, the guarded A16/D16/A17 and rank13--16 seals,
  order16/D10/distal17 and every order11 section;
* all five scalar witness extensions, also traversing any copied
  physical corners they contain;
* the280+52 positive chosen-D12 outputs. Their original constructive
  source only outputs connected small/D14 payments or the specified
  A-small/A11 and A-small/D11 roles. The worker verifies they contain
  neither a D16 output nor a marked A4/D13 combination.

For those chosen outputs it first requires a nonempty list of one or
two actual role records. Each must have an explicit canonical type,
unique integer node list of the matching rank, and a mark in that list.
The exact producer grammar is A_r(r>=1), B_r(r>=2), D_r(r>=4) orE6;
missing/None, `unclassified` and unknown labels fail before the target
exclusion test. A missing type can never masquerade as a non-target.

Each declared newer physical collection must yield corner rows; a
present hash alone is insufficient. Targets found in copied records
must still match a complete actual geometry or its half-only route.
At the end every D16 selected model and every distal D13 model must
have been encountered in the physical scan. The two completed later
A15 and D17 physical refinements introduce no new positive output:
their marked changes and retained alternatives are already among the
audited existing interfaces. This compositional fact remains explicit.

This list is the claimed interface coverage; the bare54-file count is
not substituted for it. Nonfinite displayed targets are counted and
excluded from the finite-vertex refund. Unknown finite targets abort
before a certificate is written, with at most twelve concise diagnostic
records in the error output.

## Symbolic alpha and unclaimed global optimization

The frozen distal theorem proves the entire0<=alpha<=1/2 domain with
the same selected/contact geometry. The saved maps, radical seals,
basic/short/triangle unit sets and endpoint role do not depend on alpha.
Only the selected ideal weight1/2-alpha, ideal D total113-4alpha,
rank16 envelope(113-4alpha)/17 and finite marginal alpha change.
Other ideal envelopes are unchanged; competing B111 meets D only at
alpha=1/2. Terminal margins therefore gain4n16*alpha/17 relative
to the no-distal baseline. Their exact counts/minimum still require
the appropriate profile arithmetic; this worker computes neither.

The same modeled ideal output is not later debited below1/2-alpha.
The exact current chronology/role argument is:

|Stage|Relevant action|
|---|---|
|`fresh_order17_d5_alternate_half_zero_transfer_20261018`|Final prior order17 table, D total113; all four distal entries still1/2. Earlier order17 changes are already included here.|
|`fresh_order17_d13_shortleaf_seven_thirtyseconds_scout_20261018`|Changes those four distance12 entries to7/32, interpreted parametrically as1/2-alpha, and only their endpoint/short finite counterpart to1/2+alpha. The new worker checks this paired local record.|
|`fresh_a16_rank17_spectator_extension_scout_20261018`|Subsequent extension changes finite A16 and named finite restorations only. Its D13-short A4 restoration uses INNER marks, not the endpoint role. Ideal weights are inherited unchanged.|
|`order11_fractional_cusp_closure_research_20261018` and its local packet|Subsequent ideal changes concern order11 B10/D10/cycle11 only; new finite outputs use B8/B10/D10. No order17 weight changes.|
|Completed A15 absence and D17-node2 capacity|Finite marked outputs/physical activation only; no ideal changes. The proposed D16 absence step is also finite-only.|

For each matched full four-black core the other three corner principals
are positive definite; the seal adds the same positive forest E. Those
are therefore actual finite corners, so the later order11 ideal changes
cannot affect an untracked second ideal output of this quadrilateral.
At the one ideal corner the17-node affine factor is whole: its positive
kernel forces every other incident wall orthogonal to all its roots.
It cannot secretly enlarge to an order11 or larger affine component.
Thus the parametric whole-face identity uses the retained current ideal
weight at that SAME actual marked corner. This chronology and geometry
remain explicit written inputs; the saved local-stage check alone is
not claimed to prove absence of arbitrary future ideal changes.

The half-only source registry likewise has a specific scope: all22
historical packets, seven literal source collections and pre-distal
newer sections predate this marginal alpha. The five scalar extensions
only widen already-sealed witnesses or change disjoint typed coupons;
the sole post-distal A16 rank17 extension has the inner/end distinction
stated above. Post-distal order11 sections are deliberately omitted from
the missing-amount half-only fallback. Any target there needs an
explicit amount at most1/2 or a matched full distal geometry. This is
why the registry is a reviewable written source premise rather than an
unconditional inference from a filename or absent numeric field.

The physical capacity replaces uniform4alpha*a4*d13 by at most
2alpha*[a4*d13>0], valid throughout the same parameter domain.
At alpha=1/2 that extra is at most1. The old specific D13 forest
would then cost at most11429/32+1=11461/32=358.15625, but this is
NOT a complete all-forest maximum. The D16 candidate saves
d16*[a1=0]. Both proposed steps are pointwise nonincreasing relative
to their respective old uniform payments and add no positive family.
The independent A1+2D16+D24 formal370 obstruction remains.

Any PASS here establishes only this saved-output/model interface,
conditional on complete actual-face routing, whole-factor meanings,
full arbitrary-contact seals, the rank lemma and physical capacity.
It does not prove these geometric premises by enumeration or optimize
the complete finite forest domain. No new dimension bound is claimed.

## Coordinator-only generation and separate replay

The worker owns NO inner flock. Use the reviewed parent-lock wrapper:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- \
  scripts/check_d16_d13_physical_output_retention_20261019.py --write
```

Default `--mode both` creates
`results/d16_d13_physical_output_retention_20261019.json` exclusively.
Separate replay uses `--verify`. Optional `--mode d16` or `--mode d13`
retains the same full input scan but separates local-model/target work;
its certificate has the corresponding mode suffix. Modes are not
silently mixed at replay. Canonical JSON avoids tuple/list or timestamp
differences, and the certificate is capped below1MiB.

The wrapper retains one shared lock,6GiB/40GiB start gates,3GiB/8GiB
abort gates,59-second outer limit,1GiB live process-tree RSS cap and
parent-death handling. The child adds a45-second alarm. Timing/RSS
have not yet been measured for this new source; the previous analogous
54-file saved-interface audits took under a second, and this version
adds four fixed rational19-root models and one18x18 inverse. A few
seconds is an estimate, not execution evidence. No audit is launched
by preparing this source.

## Pre-execution review revision

The root read the preceding362-line source at hash6f77b66c...c2356 and
requested two improvements before any execution. The current source
adds strict positive chosen-output type/node/mark schema validation and
checks the exact paired ideal-stage marks/payments for every D13 model.
The ideal chronology and half-only source semantics above address the
second mathematical review point explicitly. No certificate or ledger
ever bound the earlier prepared worker, and neither version has run.
