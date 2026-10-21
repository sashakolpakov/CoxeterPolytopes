# Manuscript verification

The [result–manuscript table](result_manuscript_correspondence.md) assigns
each result one primary paper. The machine-readable source is
[manuscripts/catalogue.json](../manuscripts/catalogue.json). Older drafts
and large certificate atlases are supporting archives, not additional
publication claims.

There are six current papers. Their division by technique is explained in
the [paper architecture](paper_architecture.md). Every proof-relevant
computation should follow the
[mathematical algorithm contract](computational_proof_contract.md): define
the domain and algorithm, prove termination and soundness/completeness,
then state the exact computed proposition. A successful run alone does
not discharge the preceding obligations.

From the repository root, list the available IDs and inspect a plan:

```sh
python3 -B scripts/verify_manuscript.py --list
python3 -B scripts/verify_manuscript.py ideal_ra4_volume_and_structure --plan
```

Use the actual ID printed by `--list`. A plan starts no proof process and
does not claim verification. The preserved predecessor of the combined
ideal-four manuscript has two full written consolidation reviews and a
16-page/five-figure PDF with zero TeX warnings. Its four core script checks
and conditional Lean interface pass in the
[final execution record](../results/ideal_ra4_unified_core_verification_final_20261018.json).
The current publication source adds explicit finite-model algorithms;
its [fresh full source audit](../results/review_ideal_ra4_publication_full_finite57_20261020.md)
passes at the declared inherited-input boundary. Its source-bound PDF
build now passes: 17 pages with zero warnings, as recorded in the
[six-paper PDF receipt](../results/six_paper_pdf_completion_20261020.json).
A fresh combined script replay for this publication revision is not
claimed. The earlier execution record applies to the preserved source,
not automatically to this revision.
The [ideal4 proof map](ideal_ra4_volume_and_structure_proof_map.md) gives
the result-to-section correspondence, published inputs and the dependency
that makes the chain theorem applicable to the second-volume equality case.

Run the declared checks serially:

```sh
python3 -B scripts/verify_manuscript.py ideal_ra4_volume_and_structure --mode scripts
python3 -B scripts/verify_manuscript.py ideal_ra4_volume_and_structure --mode logic
python3 -B scripts/verify_manuscript.py ideal_ra4_volume_and_structure --mode core --pdf
```

`scripts` replays the core solver-free checks. `logic` compiles the listed
Lean modules and runs their explicit `#print axioms` commands where
present; the plan records which targets contain such an audit. `core`
does both; `full` also
includes historical replays explicitly marked full in the catalogue.
The runner can reuse an older object with the explicit status
`LEGACY_CACHED_UNCAPPED`. The historical eight-target legacy cache set retains
this provenance; a successful cached exit is not a fresh capped compilation
or a declaration-level axiom audit. The seven recently verified kernels
have actual compilation and parsed axiom-output receipts. The
[selected-command source audit](../results/review_six_paper_selected_verifier_modes_20261020.md)
lists these targets and checks that the current scientific commands
cannot select a plan or generation branch instead of their declared
verification. Compact23's source-integrity check and compactRA4's
typesetting-only selection retain their explicitly narrower scopes.
For finite57 this covers31 inherited higher-tail obligations: the
ambient59--163 layers and the independent all-rank support, failure,
packing and compact-facet checks. The current full selection has113 script
jobs, nine conditional Lean roots and optional PDF compilation (123 jobs
with PDF), recorded in the actual
[current shape-counting verification plan](../results/publication_shape_counting_verification_plan_20261020.json).
Catalogue `3bc62bbd…` replaces the primary finite57 standalone
`ShapeCurvatureArithmetic` root (ambient59 arithmetic) with the verified
dimension58 shape kernel; the 11 core/123 full job counts are unchanged.
The old root's source, evidence and historical `finite57_atlas` selection
remain preserved.
The [preceding 123-job plan](../results/publication_audited_forest_tail_verification_plan_20261020.json)
added the two ancestor-forest phases and tail scalar kernel. The earlier
[120-job plan](../results/publication_complete_partitioned_verification_plan_20261020.json)
is preserved.
That preceding plan's [reviewed catalogue delta](../results/finite57_distance9_catalogue_integration_20261020.json)
replaces only the old108 selection with a scanner contract and56 actual
whole-interval commands; it preserves the old registry entry and every
other paper selection. The preceding64-job plan is preserved as history.
The earlier [44-job publication plan](../results/six_paper_publication_plans_20261020.json)
predates the D13 parser/19-part replacement.
The saved [43-job plan](../results/finite57_full_replay_tail_plan_20261020.json)
predates the addition of the verified concentration kernel. Of the 31 tail
additions,30 have bounded runtime PASS records in the
[higher-tail ledger](../results/finite57_higher_tail_replay_batches_20261020.json):
65.182 seconds total and256,278,528 B maximum RSS, including the
[remaining-check ledger](../results/finite57_higher_tail_remaining_replay_20261020.json).
The remaining tail108 obligation (covering109–116) now has a complete
56-batch replacement replay; its old59-second timeout remains preserved.
Each retains the
59-second/1GiB script cap; neither timeouts nor low memory authorize
relaxing it. The current core has eleven jobs with PDF. This is selected
certificate replay, not regeneration of
the entire discovery history; the
[interface audit](../results/review_six_manuscript_verification_interfaces_resource_20261020.md)
identifies the boundaries.
The [two ancestor-forest phases](../results/finite57_ancestor_forest_replay_execution_summary_20261020.json)
now actually PASS: 8114 tables and 10,195,241 transition inequalities in
5.106s total, with 318,799,872 B maximum process-tree RSS. Their
[contract](../results/finite57_ancestor_forest_replay_contract_20261020.md)
replays the complete coarse and first increasing-family calculations;
local geometry, actual refinement and some ambient59–65 forest/local-census
computations remain inherited finite-data inputs. No combined 123-job runtime
PASS or whole transitive regeneration is claimed. The
[fresh composition review](../results/review_finite57_theorem_composition_adversarial_20261020.md)
checks the mathematical implications at that declared boundary.
The [original 44-job full run](../results/finite57_publication_full_replay_20261020.json)
passed its first endpoint job, then failed at the eager D13 replay's 1 GiB
RSS limit. There is no full runtime PASS. Its
[memory-bounded 19-part replacement](finite_volume57_lemma_map.md#memory-bounded-d13-replay-contract)
preserves the original certificate and complete field/pin coverage. Its
[current runtime ledger](../results/finite57_d13_batched_replay_runtime_20261020.json)
records PASS for all19 parts:24.408 seconds total,184,041,472 B maximum
RSS. The current D13 parser fixture `2f6b63da…` passed; the
[earlier fixture execution](../results/d13_json_projection_contract_execution_20261020.json)
is preserved as history. The historical distance9 fixture `abc3b918…`
also passed. Worker `126d5fe3…` later actually passed its scanner fixture
(1.721s,89,882,624 B maximum RSS) and metadata (13.832s). Parser tests do
not prove mathematical replay.
The [coverage validation receipt](../results/partitioned_replay_coverage_validation_20261020.json)
records prior verifier `41f3c82b…` self-test PASS (14 adversarial
mutations) and complete D13 partition coverage19/19, without rerunning
mathematics. The [distance9 count repair](../results/distance9_symbolic_boolean_count_repair_20261020.json)
fixes a SymPy Boolean sum TypeError with the identical predicate. The
`126d5fe3…` domains replay then failed after2.709s on three expected
tuple/list key comparisons. The [source-repair review](../results/review_distance9_expected_key_json_representation_fix_finite57_20261020.md)
records those three expected tuple-to-list changes and the complete20-call
inventory, without changing any mathematical domain, predicate or ordering.
Current worker `f5f0b470…` now has an actual [scanner-fixture PASS](../results/distance9_json_key_lists_scanner_fixture_execution_20261020.json)
under fixture `f6fd1267…` (1.732s,89,915,392 B), without mathematical
interval credit. Current metadata and domains phases actually PASS
(13.792s and5.369s), as do full `overlap-domain` generation, disjoint joins,
inertia and scalar profiles. All384 classified overlap owners and949
isolation cases also PASS:2196/2196 units in56 capped batches,
388.008s total and159,432,704 B peak RSS. Complete108 replacement replay
passes. Old-source
fixture and metadata passes earn no current coverage. Coverage checker
`456abc2b…` and coordinator `cb776f59…` have metadata-only self-tests.
No full-manuscript runtime PASS or raised cap is claimed.
The current combined D13/distance9 metadata check PASSes:
D13 covers19/19 and current distance9 covers2196/2196 units, including
all overlap classifications and isolation cases. The old-source pilot is correctly rejected
as current coverage. This is receipt validation,
not a mathematical rerun.

For reproducible108 replay planning and one-batch resumption, use the
[coordinator](../scripts/resume_distance9_replay_20261020.py):

```sh
python3 -B scripts/resume_distance9_replay_20261020.py \
  --ledger results/finite57_distance9_partitioned_replay_runtime_20261020.json --plan
python3 -B scripts/resume_distance9_replay_20261020.py \
  --ledger results/finite57_distance9_partitioned_replay_runtime_20261020.json --max-jobs 1
```

Optional `--width PHASE=N` changes interval width only and requires timing
pilots; the serial gates, locks and59-second/1 GiB caps are unchanged.
A job-limit exit is not complete proof coverage, even with exit code zero.
Only same-source PASS receipts earn coverage; stale and failed records
remain preserved. The coordinator itself is not a scientific check. The
[batch contract](finite_volume57_lemma_map.md#resuming-the108-distance-nine-replay),
[scheduling source note](../results/distance9_replay_chunk_schedule_source_note_20261020.md)
and [assumption-boundary audit](../results/review_six_paper_formal_assumption_boundaries_resource_20261020.md)
explain ownership, pilot sizing and the remaining geometric inputs.

`--pdf` compiles a fresh temporary PDF and discards it, preserving the
existing final PDF and its ledger. To publish a first PDF, use the bounded
[builder](../scripts/build_illustrated_ra_manuscript_20261018.py) without
`--check-only`. No mode runs an optimizer or launches a parallel proof farm.

The current right-angled dimension package is `right_angled_sharp`.
One five-page paper proves the sharp ceilings8/general and4/entirely
ideal. Its preserved [four-job core run](../results/right_angled_sharp_core_verification_20261020.json)
passes two exact workers, one capped Std Lean module with12 printed
theorem axiom audits, and a temporary PDF build. The cubic scalar kernel
adds ten actual axiom prints and the projection/height kernel adds eight.
The current core selection has three Lean roots and six jobs with PDF. Use
`right_angled_sharp --mode core --pdf`; the
[proof map](right_angled_sharp_dimensions_publication_proof_map_20261021.md)
records all geometry and external statements that remain explicit premises.
The [preserved predecessor six-job execution](../results/right_angled_sharp_projection_height_core_verification_20261020.json)
PASSed in5.439s total outer job time. Its
[actual output validation](../results/right_angled_sharp_projection_height_core_axiom_validation_20261020.json)
checks all30 printed declarations against the three source printer lists
and accepts only the standard axioms. These repeated audits add no new
declarations to the separate seven-kernel/62-print count. The temporary
five-page PDF has zero warnings; the published PDF is preserved.
The [fresh adversarial written review](../results/review_right_angled_sharp_fresh_adversarial_finite57_20261020.md)
adds a full mathematical read at that external-input boundary, without
new executions. A later
[foundations review](../results/review_right_angled_sharp_foundations_adversarial_20261020.md)
independently audits the same geometric and published-input boundary.
The [published-input audit](../results/review_sharp_ra_published_inputs_dufour_20261020.md)
checks the mixed-vertex 27-facet consequence in Alexandrov's primary text
as the paper declares; Dufour's original remains paywalled. It also checks
Vinberg's original all-wall implication. The
[ideal-branch source audit](../results/review_sharp_ra_published_inputs_ideal_branch_20261020.md)
checks Musin's cardinality bound and Kolpakov's facet minimum with metric
rigidity in their original papers.
The ideal4 successor also has a focused
[equality/structural-chain review](../results/review_ideal_ra4_structure_batch_guard_20261020.md).

The preserved earlier dimension10 package is `mixed10`; its12-page
PDF,four exact checks and one conditional count module have historical
records in the [proof map](right_angled_dimension10_proof_map_20261020.md).
Its actual six-job run passes. Use `mixed10 --mode core --pdf`;
the older `mixed11` package is preserved as a predecessor.
The current `compact20` package uses the universal-only successor
`manuscripts/compact_arithmetic_dimension_bound.tex`: three exact checks
and six unchanged universal-theorem Lean roots. Its current 17-page PDF
compiles with zero warnings and matches its source hash in the
[PDF receipt](../results/six_paper_pdf_completion_20261020.json).
This build does not claim a fresh combined core replay.
The [field-closure source audit](../results/review_compact_arithmetic20_field_closure_batch_guard_20261020.md)
(`994b4965…`) passes at its written-proof and declared finite-data
boundary, without parked19 inputs or a new finite/Lean execution. The original
combined source, both parked field19 appendices,
20-page PDF and source-bound custom builder remain available as
`compact20_with_parked_refinements`. That historical package keeps its
original eleven-job scope and76-input manifest; it is excluded from the
current universal paper.

The compact23 package, `compact_farthest_bounds`, now selects its
conditional accounting kernel and optional PDF build: two core jobs with
PDF, or one without. Its full tier additionally checks source bytes,
not mathematics. The `compact_ra4` core still performs typesetting only;
without `--pdf` it has no jobs and correctly refuses to claim a pass.

The seven kernels in the
[publication logic map](publication_logic_formalization_20261020.md)
have now compiled under the bounded runner. The
[logic verification receipt](../results/publication_logic_verification_20261020.json)
records 62 actual axiom prints: ten cubic, ten concentration, eight
compact23, eight residual Bellman, eight projection/height and nine
finite-volume tail, plus nine dimension58 shape prints. The former six/53
receipt is preserved. Only `propext`, `Quot.sound` and `Classical.choice` occur,
with no `sorryAx`. They are included in verified conditional catalogue
coverage. Their geometric and script premises remain explicit.

The [residual Bellman certificate module](../formal/CoxeterProofs/ResidualBellmanCertificate.lean)
(`cd2a05a6…`) compiled in0.54s and passed its axiom audit in0.41s under the
unchanged bounded runner. Its [eight actual printed declarations](../results/ResidualBellmanCertificate_axiom_summary_20261020.json) use only
`propext` and `Quot.sound`. Its
[assumptions map](residual_bellman_certificate_logic_20261020.md) describes
the local integer-table conditions and derived exactness/unreachability.
Concrete alphabet/context completeness, actual table validation and
geometric interpretation remain inputs; this compilation is not a new
concrete table replay.

The [projection/height module](../formal/CoxeterProofs/RightAngledProjectionHeightKernel.lean)
(`2fac8028…`) compiled in3.47s and passed its separate3.15s audit. Its
[eight actual prints](../results/RightAngledProjectionHeightKernel_axiom_summary_20261020.json)
match the recorded1,057-byte stdout and current compiled object. It derives
coordinate Lorentz projection identities, strict wall signs, squared-height
recurrence, arbitrary finite descent and the scalar radial gap. Positive
normalizers, actual incidence/facet and future-sheet interpretations,
ridge-neighbor completeness, face inheritance and radial/tensor realization
remain inputs; no end-to-end geometric Lean proof is claimed.

The [finite-volume tail kernel](../formal/CoxeterProofs/FiniteVolumeTailKernel.lean)
(`b8f03fd4…`) has actual compilation and nine parsed axiom prints in its
[execution record](../results/finite_volume_tail_kernel_execution_20261020.json).
It checks all 105 bounded-tail cap comparisons by kernel `decide`, derives
scaled payment subtraction, and closes the near/far Lorentz scalar
contradiction. Its [assumptions map](finite_volume_tail_kernel_assumptions_20261020.md)
retains the geometric/table interpretation, loss estimates and
compact-facet alternative; it is not a full formalization of the n≤58 tail.

The [dimension58 shape kernel](../formal/CoxeterProofs/FiniteVolume58ShapeKernel.lean)
(`40d12b57…`) compiled in 0.59s and passed its nine-print audit in 0.55s;
the [summary](../results/FiniteVolume58ShapeKernel_axiom_summary_20261020.json)
and [execution receipt](../results/finite_volume58_shape_kernel_initial_execution_20261020.json)
bind the unchanged source and eight harmless lint warnings. Its
[assumptions map](finite_volume58_shape_kernel_assumptions_20261020.md)
separates derived oriented-link/extrema counting from geometric incidence
and degree/denominator interpretation. The older `shapeLower` record is
not automatically instantiated; the [original Khovanskii audit](../results/review_finite57_published_shape_inputs_20261020.md)
checks the published inputs separately. All seven kernels have raw output
and fourteen immutable compile/audit snapshots in the master receipt.

The current core selections, including PDF, are:

|Manuscript ID|Exact script jobs|Lean roots|PDF jobs|Total|
|---|---:|---:|---:|---:|
|`finite57`|1|9|1|11|
|`compact20`|3|6|1|10|
|`compact_farthest_bounds`|0|1|1|2|
|`right_angled_sharp`|2|3|1|6|
|`ideal_ra4_volume_and_structure`|4|2|1|7|
|`compact_ra4`|0|0|1|1|

These are current plan counts, not a claim that every combined selection
was rerun after the catalogue changes.

Add `--ledger results/my_new_verification.json` to retain one compact
execution report. The path must not already exist. A failure or resource
deferral is recorded as such; it cannot produce a successful manuscript
verification. Source, catalogue and checker hashes identify the material
checked. Sources and old frozen evidence are never overwritten by a
verification run.

## What a successful check means

There are three separate layers:

|Layer|What is checked|What remains assumed|
|---|---|---|
|Exact scripts|The declared finite tables, rational inequalities, model incidences or provenance checks|The written proof that the enumeration covers the intended geometric problem; correctness of any trusted enumerator|
|Conditional Lean|The stated implications and arithmetic from typed hypotheses|Published facts, explicitly named script outputs, and geometric lemmas not yet formalized|
|PDF build|The manuscript typesets; the build report records warnings for review|Every mathematical assertion|

The [assumption register](manuscript_assumptions.md) names published
inputs, script outputs and remaining internal geometric interfaces.
Its source-handoff compilation note predates the successful ideal-four
and compact20 batches; current per-module status is in the catalogue and execution
records. The premise text remains the version hashed by that run.
Internal geometry is not relabeled as a published fact. Assuming a correct
enumeration is allowed, but its output and its bridge to geometry must be
separate, precise hypotheses. Assuming the desired main conclusion itself
would make the verification empty and is not an acceptable interface.

Lean theorem parameters are assumptions even when `#print axioms` lists
only the usual kernel axioms. Read the theorem signature and the register
together. A successful conditional proof does not establish an end-to-end
Lean formalization of the geometric manuscript.

## Resource control and current execution status

Each script batch has a59-second limit, a1GiB live process-tree RSS limit
and one computation thread. Older non-locking checkers inherit a lock
held by the common wrapper. Lean uses the existing serial module runner
with a2048MiB cap and dependency/source ledgers. PDF batches use the
bounded builder. Temporary build outputs are removed after each batch.

Every batch requires40GiB free disk and6GiB available RAM to start and
stops below8GiB disk or3GiB available RAM. Below20GiB disk, the coordinator
starts a cleanup agent. These limits are not relaxed to finish a build.

All six current PDFs match their source-bound build ledgers in the
[audit-edition receipt](../results/six_papers_audit_edition_20261021.json):
finite57 has48 pages, compact20 has17, compact23 has36, sharpRA has5,
ideal4 has17 and compactRA4 has14. Finite57 source `c300db69…ca82a0` /
PDF `1f3c020c…` has main exposition pp.1–12 and essential appendices/table
documentation pp.13–48, with21 underfull spacing warnings and no overflow
or unresolved references. The other five have zero warnings. The44-page
`89bd7658…` source / `ebf996a3…` PDF and41-page phase remain archived.
Compact23 source `530d964e…` / PDF `1e6b7155…` retains the accepted
correlation correction and removes stale workflow/audit-pending prose;
its [editorial review](../results/review_compact23_editorial_scope_cleanup_20261021.md)
and [new build](../results/compact_farthest_audit_edition_editorial_build_20261021.json)
record the unchanged mathematical scope. SharpRA source `0c29a786…` /
PDF `609a9fc0…` adds face-inheritance and height-notation explanations,
with a [five-page zero-warning build](../results/sharp_ra_audit_edition_readability_build_20261021.json).
These source successors have no new combined script/Lean runtime claim.
[Selected finite57 pages5,16,24,25](../results/finite57_publication_delta_visual_review_20261020.json)
passed visual inspection; no all-page inspection is claimed. Typesetting
does not establish mathematical claims or add script/Lean passes. The
separate logic receipt above records the seven kernels' actual compilation
and axiom audits.
The [finite57 repair receipt](../results/finite57_full_budget_bridge_repair_20261020.json)
binds the independently reviewed full-budget classical unit-charge bridge
before the first-failing-factor split. The preceding 41-page source/PDF
are archived. This is an exposition repair within the existing proof,
not a new dimension bound or a new finite-enumeration execution.
The later [publication bridge receipt](../results/finite57_publication_bridge_completion_20261020.json)
records two further independently accepted insertions: current triangle
exclusions and complete coarse forest-screen domination. Current source
`c300db69…ca82a0` now has an actual bounded rebuild PASS after the earlier
no-child resource deferrals. Its stored-positive-entry tight-witness
condition and rank-zero exception are explicit. No full replay or
additional dimension theorem follows. D13 is now19/19 by the separate
bounded replay ledger; all30 other original tail commands and the complete
108 replacement have actual scoped PASS records. The TeX builds themselves
add no mathematical script or Lean pass.

### Earlier verification history (18--20 October2026)

The orchestration and legacy-lock extension passed Python syntax and
catalogue plan checks. A real below-threshold run recorded
`DEFERRED_RESOURCE_GATE`, without launching a proof. When memory rose,
the ideal-four core run passed all four exact replays, compiled
`ManuscriptIdeal4`, and printed its theorem axioms: only `propext`,
`Quot.sound` and `Classical.choice`. The older imported objects retain
their explicitly reported legacy kernel-cache provenance. A toolchain
recursor-name error in the first attempt was corrected before the final
successful run; no failed object was published.

The16-page ideal4 PDF build took1.914seconds and84.3MB peak process-tree
RSS, with zero warnings. The compact20 PDF first deferred below the
memory gate, then compiled when memory recovered:8 pages, zero warnings,
82.2MB peak. Its [core run](../results/compact20_core_verification_20261019.json)
passes all ten jobs: three exact replays, six declared conditional Lean
module checks, and a temporary PDF rebuild. The new compact20 theorem
audits list only propext and Quot.sound, with no sorryAx. Older imported
caches retain their explicitly reported legacy provenance. All temporary
build outputs were removed. These records concern the preserved revisions;
the current publication PDFs and new logic kernels have the separate
receipts above.

The preserved combined compact arithmetic predecessor expanded to20 pages:
universal20, AppendixA for eight actual fields at19, and AppendixB
adding actualK7 to that scoped19 theorem. The
[new source-bound build](../results/compact_arithmetic_nine_field19_included_build_20261020.json)
passes with zero warnings and all76 combined input pins checked before
and after execution. The new K7 conditional kernel and its separate
six-theorem axiom audit also pass. Its
[proof map](compact_arithmetic_nine_field19_proof_map_20261020.md)
lists all internal geometric imports and remaining Lean premises.

The historical combined-version plan is:

```sh
python3 -B scripts/verify_manuscript.py compact20_with_parked_refinements --mode core --pdf
```

It declares eleven serial jobs: three exact checkers, seven Lean roots
and a source-bound temporary PDF rebuild. K7 imports the eight-field
kernel as a checked dependency; its previous audit remains preserved
instead of being redundantly overwritten. The first new combined run
passed the three exact checks and three cache checks, then deferred
before the next audit when available memory fell below6GiB. This
[deferred record](../results/compact20_nine_field_core_verification_20261020.json)
is retained as such. After available memory recovered above8GiB, the
[final combined run](../results/compact20_nine_field_core_verification_final_20261020.json)
passed all eleven jobs. The source-bound temporary PDF had20 pages
and zero warnings; it was discarded and the published PDF was preserved.
All actual geometry remains at the explicit conditional Lean boundary.
