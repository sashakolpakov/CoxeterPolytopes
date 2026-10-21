# Current ideal5/6 logical verification

19 October2026. The latest one-module Lean compilation and its separate
printed-axiom audit **PASS**. This is conditional logical verification,
not a formalization of hyperbolic polytopes or proof of their existence.

Source:
[RightAngledElevenVertexCounts.lean](../formal/CoxeterProofs/RightAngledElevenVertexCounts.lean),
SHA256`79846f2cf509e7262ca109e0a7ea2b79fa8cc421df772e9b77dc77116823e04d`.
The [source preparation and assumption inventory](right_angled_eleven_vertex_counts_lean_preparation_20261019.md)
remains the historical handoff; the execution below supersedes its
uncompiled status.

|Kernel-checked implication|Explicit premises still requiring geometric or published interpretation|
|---|---|
|Carrier has at least37 vertices and25 facets|Inherited primitive +3 cusp-defect counts and large-three-facet bridge; new internal equality exclusion when that three-facet has11 vertices.|
|Ideal5 has at least63 facets and211 vertices;3V>=631|Actual facet injection and the primitive facet/cusp propagation identity and inequalities, not the old final62/206 bound.|
|Ideal6 has at least275 facets,272V>=6039F andV>=6106|Complete list of actual five-facets with the preceding interfaces, actual injection, exact coface incidence and the universal reciprocal-transport inequality.|
|Mixed11/12 fractions47I5<4N5 and248I6<15N6|Complete ideal/nonideal face lists, the internal ideal-face structures, published other-face minima16/27 and strict NK averages20/42. Positive denominators are derived.|

Compilation took0.45s; the separate axiom audit took0.40s. The
[compile ledger](../formal/.proof-batches/CoxeterProofs_RightAngledElevenVertexCounts.json)
and [axiom-audit ledger](../formal/.proof-batches/CoxeterProofs_RightAngledElevenVertexCounts_axiom_audit.json)
bind the final source, compiler and three capped cached dependencies.
Start resources were187.13GiB disk and7.43/7.44GiB available RAM.
The batch used one thread and2048MiB, removed temporary outputs, and
printed only`propext`, `Classical.choice` and`Quot.sound`; exact rounding
uses no axioms. There is no warning or`sorryAx` in the successful run.

The first attempt stopped at the resource gate. A subsequent compilation
found a tactic atom mismatch between an abbreviation and its expanded
field, plus an unused simp argument. An explicit type on the equality
exclusion and removal of that unused argument fixed both. These were
elaboration changes; no premise or mathematical conclusion changed.
The failed run is not verification evidence. Only the final source hash
and successful compile/audit ledgers above establish the current scope.

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.RightAngledElevenVertexCounts --audit-axioms
```

The actual cusp, spinor, small-facet bridge, eleven-vertex equality and
transport geometry remains reviewed written mathematics. Internal
lemmas are explicit hypotheses here, separately from published inputs;
they are not relabeled external axioms. No final63/211 or275/6106 size
bound is simply assumed. The module proves those consequences from the
displayed interfaces.
