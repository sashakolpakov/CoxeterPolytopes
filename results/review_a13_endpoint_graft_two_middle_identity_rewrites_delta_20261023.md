# A13 scalar kernel: complete two-rewrite successor review

23 October 2026. **SOURCE/DELTA PASS; no compilation or audit PASS
is claimed.** I read the entire current214-line source and159-line
map, the preserved predecessors, exact complete diffs, archive manifest,
actual failed-compile record and its full6359-byte Lean diagnostic.
No compiler, worker/import, test, plan, build or Git operation ran.

## Exact current binding and preservation

| Object | Bytes | SHA256 |
|---|---:|---|
|formal/CoxeterProofs/A13IEndpointGraftScalarKernel.lean|9974|211a8c1ab98f3e820a957ed5e2966a1c818400d1082701b6abdd2278e0410792|
|docs/a13_i_endpoint_graft_scalar_kernel_logic_20261023.md|9277|d75dedd025dc0db321a9215eb546d7a121a140bcb3c839fa20d6c9230f1023b2|

The preservation manifest is
`results/source_archive/a13_endpoint_graft_before_middle_identity_rewrites_20261023/manifest.json`.
It is1226 bytes, SHA256
`7d1e16a8e91fe63cb9e1b6eee8916de258555e686cc2e52fa19c1d6b9662f6ef`.
Its source gzip is2410 bytes, SHA256
`d33df05e87ff1c5e57603c2655218be36679a4a2fb62d3e5da94104718066c92`;
decoding independently recovers EXACTLY9892 bytes/212 lines with
SHA256 `69bd0c29d6e9661ec6428ed2ec0995eecb72aa8a62d611b88bc76577376b637e`.
Its map gzip is4186 bytes, SHA256
`e01caf333e12eb89512ec9b01eac1d63c918ef5418c0404374ace39908d50c33`;
decoding independently recovers EXACTLY8728 bytes/151 lines with
SHA256 `fa482d15f8743d96394af4e89ed845db207503770e486121ff5cc38d4db882b2`.
Both decoded lengths and hashes match the preservation manifest.

I authored the original source/map; Sharp supplied their full different-
author static review. ROOT authored the present two-line repair. This
report independently reviews THAT repair and its precise successor/map
boundary; it is not a replacement independent review of my original code.
The existing full Sharp review remains its inherited pre-repair gate.

## Entire source delta and algebra

The ENTIRE source diff consists of adding the SAME line in exactly TWO
middle branches, immediately after obtaining the scalar inequality:

    rw [← middle_identity] at scalar

In `nonhelper_all_ranges` this rewrites the squared deficit to
2500000*data.middleGain, so the local fact becomes

    217000 < 2500000*s^2 + 2500000*data.middleGain.

Positive scaling by1000 gives exactly the next requested statement
`217 < 2500*(s^2+data.middleGain)`. In `helper_all_ranges` the analogous
local fact becomes

    385000 < 5000000*s^2 + 5000000*data.middleGain,

whose positive scaling by5000 gives exactly
`77 < 1000*(s^2+data.middleGain)`. Thus the rewrites expose identical
linear atoms to the existing ordered-ring tactic. They neither change
the scalar identities nor add a premise about the desired inequality.
The local equality is already an explicit theorem argument in each case.

Every theorem signature, `PaidGraftData` field, branch domain, square,
paid-pressure/norm interface, final contradiction and all ELEVEN print
declarations are byte-unchanged. The Std-only generic ordered-ring
coordinate model is unchanged; there is no Nat proxy, geometric instance,
fallback theorem, extra assumption, declared axiom, proof hole,
`native_decide`, Mathlib or scientific-worker import. Actual geometry,
field/inventory facts, original-wall feasibility, legal individual prices
and full joint weak-foot absorption remain external and explicit.
This source-level check cannot guarantee a future elaboration result.

## Failure chronology, complete map delta and verdict

The actual admitted failed compile is recorded at
`results/lean_a13_i_endpoint_graft_verification_20261023/attempt_2_failed_compile.json`.
That1203-byte record has SHA256
`7e19e8fb0ffbffb1465cfa5d4b6415a8d7352ec1e7fa6ca674441ceb3400fdc0`.
Its full raw diagnostic is6359 bytes, SHA256
`7db0dd78c335aaf677a146f9d85325cc9ac6ca97475e96b5317794963fd1a725`.
It reports the two failed `grind` gain scalings at predecessor lines123
and147. Their scalar facts and intended gains have exactly the mismatch
the new explicit rewrites address. The diagnostic's error-propagated
`sorryAx` in the two failed all-range declarations and final descendant
earns ZERO verification credit. The unchanged unused-variable warning is
retained, not suppressed. No separate actual axiom audit was started.
The earlier4.7GiB refusal was BEFORE a Lean child and remains preserved.

The ENTIRE map diff updates the source binding/size/linecount, records
the refusal then failed admitted compile, explains the two rewrites,
and distinguishes requested prints from failed diagnostic dependencies.
All mathematical input bindings, theorem meanings, eleven-name print
list and external premise inventory remain unchanged. It correctly
requires renewed root preclear, an admitted guarded compile and a
SEPARATE successful actual audit. It asserts neither a compilation
success nor verified standard-axiom-only dependencies at this stage.

The exact two-line repair and complete map/status delta receive PASS.
No theorem statement or mathematical proof scope changed. Root's eventual
guarded run remains the execution gate; this report supplies none and
does not alter a frontier, ceiling, paper, successful module or old receipt.
