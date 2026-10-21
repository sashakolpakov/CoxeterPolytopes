# LargeInertParentArithmetic: execution and interface review

21 October 2026. **PASS: actual bounded compilation and all four actual
axiom outputs are correctly bound.** No compile, audit, mathematical job,
worker import or source edit was performed by this reviewer.

The unchanged source is
`formal/CoxeterProofs/LargeInertParentArithmetic.lean`, SHA-256
`a62fe70788b10ee6b1f47d046ab732b79b3ad2ccc5c77b069efc53d42f919c7d`.
The historical precompile review remains unchanged at
`review_large_inert_parent_arithmetic_before_compile_20261021.md`, SHA-256
`4b330b0eb211001542ff7415efedc8dd9003db20115e6a5d7fd1dc3b3a33a017`.

I read both complete execution records and streamed the source, object,
all recorded compiler/configuration files and retained output. Bindings:

| Artifact | SHA-256 |
|---|---|
| `results/large_inert_parent_compile_20261021.json` | `9c445b529831e2f988e5c30a7cae35e2f3b69966ff215ebe8fdc07bfd250c514` |
| `results/large_inert_parent_axiom_audit_20261021.json` | `2d7367d15f23a2dae764c8e82c0719d34c45bcc9397f89643eb4550bf2546ba2` |
| Current 80,520-byte module object | `7827d930595fa2fd365462f2de39f077c69b33f671623df8937f2e6cea75ae9a` |
| `results/large_inert_parent_actual_axioms_20261021.txt` | `ceb369a30de9fcaf8eac16eb79fcda90c903c0b74ee7784c1f7304666b296bc1` |
| `results/large_inert_parent_axiom_summary_20261021.json` | `539a8581a7f0e7c856677ba72ddb8456ad0cbd265c18175d2b2f9c232baa7097` |

Both immutable execution snapshots are byte-identical to their current
operational ledger counterparts. They have identical source/compiler
inputs. The pinned toolchain, lake configuration/manifest, Lean executable
and all three recorded shared libraries match current bytes. The compiled
object matches its compile record and is newer than the source. The empty
dependency-object map correctly represents **no local imports**; the only
import is Std. It is not a separately hashed transitive Std-object closure.

The retained raw output is exactly375 bytes and matches **both** child
stdout hashes and lengths. Every line has the exact form “does not depend
on any axioms.” Parsing gives precisely the four requested source commands,
in order:

* `coefficient_table`;
* `canonical_surplus2_table`;
* `remaining2_table`;
* `remaining2_lower_bound`.

All are in `CoxeterProofs.LargeInertParentArithmetic`; there are no missing,
duplicate or extra declarations, warnings, error diagnostics, `sorryAx`
or other axioms in this output. The parsed summary agrees with the raw
record and source hash. This verifies the actual outputs, rather than
inferring an axiom result merely from the source commands.

The compile took0.55s with153.33GiB free disk/6.03GiB available memory at
admission; the separate audit took0.41s with153.33GiB/6.58GiB. Both records
specify2048MiB and one thread. These are admission resources and a memory
cap, **not a measured peak RSS**. The signal-safe and base runners retain
their frozen `c0ddc98c…` and `97822ecf…` hashes and the existing40/6 start
and8/3 abort gates. The source is rechecked before a ledger is published.
The observed module products contain only its required `.olean`; temporary
work and pending ledger files are absent, and the shared lock was idle
when checked. No existing local dependency module was required by this run.

The prepared map archive has compressed SHA-256
`5fb40ed0580404feef0d171d40df1f56853f9036f98357e352ec3bc3d4991b96`
and decodes exactly to the reviewed4,141-byte map
`144b67fd3a446b4ac758bbb1cc80571c288ed97a7e92bc8b6bc6f665672d0bf4`.
The current map is4,909 bytes, SHA-256
`8508a3cb63ab7146a2c11faba71297c22415250c66165049c5eba380f30e2580`.
Its exact diff updates execution status, evidence links, theorem-table
heading and reproduction wording; no formula, domain or geometric premise
changed. It honestly distinguishes successful execution from the earlier
prepared status and from any fresh visual or geometric verification.

The compiled result is four exact, axiom-free arithmetic facts about the
fifteen signed integer formulas. Coefficient interpretation, link geometry,
convex order, geometric charge envelopes, actual incidences, finite forest
coverage and global allocation remain outside Lean, as the map states.
Neither the selected-family financing theorem nor a dimension56/global
polytope conclusion has been formalized. The independent written-proof
review remains a separate result. No correction is required in the
execution packet at this scope.
