# Corrected M6 batch: generation and separate replay

20 October 2026. **PASS for the reviewed finite interval and exact
matrix calculations.** These are research results inside the open
dimension9 candidate, not a new compact example or dimension bound.

The corrected preparation is
`golden_dyadic_m6_actual_short_edge_preparation_20261020.md`, SHA256
`2afd545c3c3166cd07381949b972fb012280b18ed673a090025b6aac30537dd2`.
The fixed worker is
`scripts/check_golden_dyadic_m6_actual_edge_20261020.py`, SHA256
`9e48ab80a337f43925a1be7ee6761b02da5eaa0c60df6dfb649bc640ec5e6949`.
Paths without a directory prefix are under `results/`.

Root ran the worker through the existing parent-lock wrapper, first
with --write, then --verify in a fresh process. Both PASS, with
unchanged source hashes. All nine profiles and 924 coefficient tests
are covered. Exactly h and g survive as integral acute normals;
`g=2h+beta1` and their common collision parameter are verified.
The reviewed redundancy proof selects h as the unique simple wall.
The endpoint, saturated B6 Gram, full integral future involution,
trace and fixed doubling identities all pass.

|Object|SHA256|Size / execution|
|---|---|---|
|`golden_dyadic_m6_actual_edge_20261020.json`|`869c2887ff4a8a3cb21571c2a1ccef109991dbe53845623c4cd6fc5759af52b0`|15,349 bytes|
|`golden_dyadic_m6_actual_edge_generation_20261020.json`|`25140a65477b1116f3764efadec3c2788e3c9148e8af5128f3cf406ee62f34c6`|0.113s; peak55,738,368 bytes|
|`golden_dyadic_m6_actual_edge_replay_20261020.json`|`10228212340df2e5dcc63d94d03a2555a23a1459de1701bbf033d3a13eb4a7ae`|0.112s; peak56,410,112 bytes|

Both jobs held the shared exclusive lock, used one thread and the
59-second/1-GiB process-tree limits. Initial available RAM was 8.00
and 8.02GiB; free disk was177.07GiB. No resource override was made.
The earlier two fail-closed attempts and the nonsimple-root correction
are recorded in the preparation and source reviews; they produced no
certificate or success ledger. No failing case was filtered away.

The complete written proof has root review, outside-author review
`review_golden_dyadic_m6_edge_resource_20261020.md` (SHA256
`6cca91d7fa7072912e3fc4c6d15e2d6b9c7739686e5d12237c5bca24ca7072ec`),
and correction review
`review_golden_dyadic_m6_simultaneous_root_correction_batch_guard_20261020.md`
(SHA256 `85dcf93805d6d2aa14994914cb69d62e2fde23a59cf66a4285791bf117544316`).
The initial pre-execution review is historical and explicitly superseded
where it expected one accepted normal rather than one simple wall.

## Geometric closeout and limits

The reviewed finite reduction plus the executed census prove h actual.
The checked integral involution maps the complete actual B6+H3 germ
to the actual endpoint germ and extends over the exact perpendicular
splitting, so it preserves the full ambient chamber. It makes w actual.
Thus all nine native M6 walls are actual. The explicit known compact
P6-double lineage then gives a compact nine-wall cone. Acute projection
places all its vertices in the full chamber; convexity proves equality.
This closes the FULL M6 wall identification without assuming that any
compact integral-root cone is automatically a full reflection chamber.

The finite-type inventory, actual-germ argument, sealed H3 geometry
and known-example compactness are complete written mathematical inputs,
not consequences of the finite arithmetic certificate. No end-to-end
Lean or fresh primary-literature classification check is claimed.
The involution has order2 and supplies no infinite chamber period.
The ambient R_5 remains OPEN for compactness and finite volume.

No Lean or TeX was run. The worker created only its retained certificate
and compact execution ledgers; -B suppressed bytecode and no temporary
build products were created. Sources and earlier evidence are preserved.
