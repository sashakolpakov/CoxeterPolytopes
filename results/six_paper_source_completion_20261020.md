# Six-paper source completion and first audits

20 October 2026. The six current TeX drafts contain their intended arguments at the finite-data and companion-theorem boundaries stated in each text. All six have internal mathematical source audits and matching current PDFs. The48-page finite57 edition includes the two written bridges and the saved-table wording clarification; its44-page predecessor is byte-preserved. The compact23 edition includes the accepted correlation-bound precision and remains36 pages. The [PDF receipt](six_paper_pdf_completion_20261020.json) distinguishes these bindings; the [bridge receipt](finite57_publication_bridge_completion_20261020.json) records the current finite57 revision, and the [logic receipt](publication_logic_verification_20261020.json) records selected conditional Lean verification. No dimension bound or known example changed in this writing phase.

| Paper and current source | Completion and audit | PDF status |
| --- | --- | --- |
| [General finite-volume57](../manuscripts/finite_volume57.tex) | Full corner rules, geometric routing, finite-forest and cusp algorithms, analytic concentration and both higher-dimensional tails integrated. The currentc300db69 revision also supplies the independently reviewed current-triangle exclusions and complete numerical coarse screen; coarse representative counts are distinguished from actual multiplicities. The residual certificate also requires witnesses at every stored finite positive rank and treats zero separately. Recorded finite enumerations remain explicit inputs. | [48-page current PDF](../manuscripts/build/finite_volume57.pdf): main1–12, essential appendices13–48. Spacing warnings only;44-page predecessor preserved. |
| [Compact arithmetic20](../manuscripts/compact_arithmetic_dimension_bound.tex) | Universal-only text, complete golden-field local rules and seed/queue/forest interfaces. [Fresh exposition audit PASS](review_compact_arithmetic20_exposition_finite57_20261020.md). Earlier field reductions remain precisely stated internal companion theorems with attached proof homes; they were not re-proved by this audit. | [17-page universal-only PDF](../manuscripts/build/compact_arithmetic_dimension_bound.pdf), zero warnings; combined20-page predecessor preserved. |
| [General compact23](../manuscripts/compact_farthest_vertex_bounds.tex) | Ordinary components, debits, actual pressures, continuation/caps and final forest argument integrated. [Full source audit PASS](review_compact_farthest23_full_batch_guard_20261020.md), followed by [fresh classical/forest audit](review_compact_farthest23_classical_forest_finite57_20261020.md) and the separate coupled-D23 audit. The correlation endpoint is correctly non-strict; strict demands preserve the contradiction. | [36-page current PDF](../manuscripts/build/compact_farthest_vertex_bounds.pdf), zero warnings. |
| [Sharp right-angled8/general and4/ideal](../manuscripts/right_angled_sharp_dimension_bounds.tex) | Complete short proof with its two branches. [Fresh full source audit PASS](review_right_angled_sharp_dimension_bounds_full_resource_20261020.md); named published inputs checked at the stated scopes. | Existing source-bound5-page PDF, one figure. |
| [IdealRA4 volume and structure](../manuscripts/ideal_ra4_volume_and_structure_publication.tex) | Gap, equality reduction and restricted chain theorem, with exact fixed-model algorithms and completeness proofs. [Full successor source audit PASS](review_ideal_ra4_publication_full_finite57_20261020.md). Published classifications and recorded exact outputs remain explicit inputs. | [17-page current PDF](../manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf), zero warnings; verified16-page predecessor preserved. |
| [CompactRA4 volume constraints](../manuscripts/compact_ra4_volume_constraints.tex) | Capacity, neighboring excess, zero-excess recognition, positive-excess gap and degree-two constraint. [Fresh full root audit PASS](review_compact_ra4_volume_constraints_full_root_20261020.md). Compact minimum and second-volume questions remain open. | Existing source-bound14-page PDF, four figures. |

The common methods and exact points where cases separate are in the [paper architecture](../docs/paper_architecture.md). There are still six paper families and32 owned TeX files; preserved revisions, input fragments and generated atlases are not extra proposed papers.

## Repairs and proof interfaces

The57 transcription was corrected to retain the D5/D11 cross-payment already present in the inherited rule. Its all-rank tail now explicitly permits the repeated omission needed for an affine rank-one factor and cites Nikulin's own unrestricted narrow-face statement. In the compact23 A-cap argument, the audit found that the possible weak middle index at N=22 needed to be included in the omitted-floor count. The written repair includes it in the same six-hole loss bound. These repairs do not alter either ceiling. Each fresh review records the mathematical source version it read. Three subsequent typography-only revisions preserve their statements/data; the [PDF receipt](six_paper_pdf_completion_20261020.json) binds those editorial deltas to byte-verified archives and current sources.

Proof-relevant algorithms now specify their finite objects and rejection tests, why the search halts, why the geometric cases are covered, and the exact recorded conclusion. Raw rows, hashes and historical atlases remain in the reproducibility package. Checking a certificate does not by itself prove that its domain covers every polytope; the source audits checked those bridges at the explicitly stated boundaries. They did not rerun the recorded enumerations. The arithmetic20 field-reduction companion theorems remain separate technical inputs, rather than being presented as published results or as consequences of a hash check.

The parked nine-field19 refinements are excluded from the universal arithmetic20 source. The old combined source/PDF, both appendices and their manifests remain byte-preserved. The [initial metadata receipt](six_paper_source_completion_metadata_20261020.json) records the initial source-completion hashes, static reference/environment diagnostics, preparation archives and preservation checks. The later PDF receipt binds current sources and PDFs. Static diagnostics are not a TeX compiler or a mathematical proof.

## Builds, reproducibility and Lean

During the initial source-completion phase, no scientific, TeX or Lean worker was started. Disk space remains well above the20GiB cleanup trigger, but available RAM has stayed below the6GiB start gate. Temporary proof/build files were therefore not created. Revised page counts will be recorded only after actual serial compilation and figure inspection.

The initial [six verification plans](six_paper_source_completion_plans_20261020.json) were generated without execution. The [verification guide](../docs/manuscript_verification.md) provides the serial entry point. Its core plans contain8,10,1,4,7 and1 jobs respectively; the two geometry-only papers currently have PDF jobs rather than scientific checkers or Lean modules. These plans describe available checks, not completed runs.

Existing conditional Lean records retain their original source and assumption scopes, including the sharp right-angled12-theorem audit. This phase adds no Lean coverage. Those steps were then completed in the publication follow-up described below, under the shared lock and resource limits. External mathematical facts and exact-script conclusions must remain explicit premises wherever the formal interface assumes them. None of these internal reviews is outside refereeing or end-to-end geometric formalization.

The publication builder enforces the same live resource gates,59-second cap,1GiB process-tree limit and exclusive lock, and removes its temporary directory. The four commands used for initial publication were:

```sh
python3 -B scripts/build_illustrated_ra_manuscript_20261018.py manuscripts/finite_volume57.tex --replace-unfrozen
python3 -B scripts/build_illustrated_ra_manuscript_20261018.py manuscripts/compact_arithmetic_dimension_bound.tex
python3 -B scripts/build_illustrated_ra_manuscript_20261018.py manuscripts/compact_farthest_vertex_bounds.tex
python3 -B scripts/build_illustrated_ra_manuscript_20261018.py manuscripts/ideal_ra4_volume_and_structure_publication.tex
```

The first replacement is allowed only for the unpinned preparation PDF whose bytes have already been preserved in the archive receipt. The other three are first publications at their new paths. Build failures or warnings must be resolved and the resulting PDFs inspected before the catalogue's PDF status is advanced.

The subsequent [logic/replay follow-up](publication_logic_and_replay_preparation_20261020.json)
prepared three source-reviewed, uncompiled kernels and added31 existing
higher-tail leaves to finite57's full plan (43 serial jobs including PDF;
the eight-job core is unchanged). The bounded57 builder was attempted
and refused at its start gate with5.4GiB available RAM, before any TeX
child or temporary build directory was created. No PDF or verified Lean
coverage changed. The [logic map](../docs/publication_logic_formalization_20261020.md)
states what those drafts derive and what remains a premise.

## Completed PDF and selected-logic follow-up

The completed earlier PDF phase matched all six source hashes:44/17/36/5/17/14
pages in the table order. The subsequent finite57 revision below is not
yet typeset. Five PDFs have zero warnings. The44-page finite57 has spacing
warnings only; its mathematical proof appendices are essential parts
of the argument, not optional raw data. Root inspected selected current
PDF pages, and a separate agent checked every source/PDF/build-ledger
binding, all seven parked-preservation pins and five typography archives.
This is not a claim that every rendered page received visual inspection.

Ten successful bounded TeX batches, including the necessary typography
repairs, produced the four current revisions. Temporary TeX products
were removed after every batch. The57 preparation PDF and the old
combined arithmetic/ideal editions remain preserved. The new builds
do not alter any dimension bound or rerun a mathematical enumeration.

Three selected Lean kernels now pass bounded compilation and28 actual
axiom prints: right-angled cubic10, concentration10, compact23 accounting8.
Initial elaboration errors were repaired without changing their theorem
statements. Parsed actual output uses only propext/Quot.sound/Classical.choice,
with no sorryAx. Source/object/compiler/import bindings and the precise
remaining assumptions are in the
[verification receipt](publication_logic_verification_20261020.json) and
[logic map](../docs/publication_logic_formalization_20261020.md).

The current core plans include these modules:9/10/2/5/7/1 jobs including
PDF. Finite57's original full plan has44 jobs, including31 selected
higher-tail leaves. Twenty of those leaves now have fresh bounded passes;
the108 replay timed out at59 seconds, and ten other leaves await successful
execution. The original full replay passed the final57 endpoint check,
then failed at the D13 loader's1GiB memory cap. It has no full PASS.
Plans and metadata are not proof executions; these partial kernel results
do not formalize the entire geometric manuscripts.

## Subsequent adversarial audits and bounded replay

The finite57 audit found that its unbounded classical tail used a
full-budget unit-charge estimate without supplying that estimate in the
text. The repair proves it before invoking the smaller residual budgets:
it gives the polarization identity, the complete B/C/D diagonal table,
the small-A finite coverage, and positive parity polynomials for the
remaining A ranks. The [independent repair review](review_finite57_full_budget_classical_tail_repair_finite57_20261020.md)
and [integration/build receipt](finite57_full_budget_bridge_repair_20261020.json)
bind the written repair to the current44-page source and PDF. This was a
missing written bridge; the audit and rebuild do not constitute a fresh
execution of the580 finite A cases.

A [fresh sharp right-angled adversarial review](review_right_angled_sharp_fresh_adversarial_finite57_20261020.md)
and a [fresh idealRA4 structural review](review_ideal_ra4_structure_batch_guard_20261020.md)
found no mathematical gap at their stated external boundaries. The latter
checks actual disjoint cuts and the passage from local octahedra to a
chain, including finite and ideal contacts; it does not classify arbitrary
ideal right-angled polytopes.

The [partial tail ledger](finite57_higher_tail_replay_batches_20261020.json)
records the20 actual passes and the108 timeout. A19-part D13 replacement
and a partitioned108 replacement preserve the original scripts and
certificates. Their complete coverage and runtime remain to be established;
both current parser/scanner fixture suites now pass, including optimized
Python rejection and the108 root-array trailing-whitespace regression.
Five D13 parts passed (the four input-pin partitions and metadata);
fourteen mathematical replay parts remain. Both successors have
complete [D13](review_d13_transfer_batches_batch_guard_20261020.md) and
[108](review_selective_distance9_batched_replay_root_20261020.md)
source-equivalence reviews, distinct from mathematical runtime passes.
RAM-gate deferrals create no scientific child and count as no mathematical
pass. Scientific, TeX and Lean jobs remain serial with unchanged limits.

The reader-facing review of the preceding edition found no
publication-readiness blocker. Its six source/PDF/catalogue/build-ledger
bindings and page counts were checked after the first bridge repair.
The reviewed insertion is byte-preserved in the repository archive and
its owned temporary original was removed; the shared proof lock is idle
and the proof work directory is empty.

## Current triangle and coarse-screen revision

A further mathematical audit identified two places where the printed
definitions did not yet justify inherited conclusions. The current
source1d8332ec supplies both: an explicit no-triangle proof for the three
reduced small affine corner types, and the full numerical upper problem
used by the3298-to39 coarse forest screen. The latter includes every
coefficient, the item generator, reservations for both remaining contexts,
the separate-family scope, and the distinction between representative
rank-two counts and actual A2 counts. Coarse k12 also differs explicitly
from actual A12 multiplicity. Both fragments and their integrated delta
have independent review acceptance; see the [completion receipt](finite57_publication_bridge_completion_20261020.json).

The source and PDF of the44-page edition were byte-preserved before this
edit. The bounded rebuild was deferred before any TeX child at3.7GiB RAM.
Five other current PDFs retain their previous matching bindings. The
current57 source therefore has no new compiled page count.

The [coverage validator receipt](partitioned_replay_coverage_validation_20261020.json)
records passed synthetic metadata tests and a correctly incomplete5/19
D13 result. It does not rerun mathematics. A source-reviewed residual
Bellman Lean kernel now derives exact maxima and unreachability from all
transitions and tight predecessors, but is UNCOMPILED; the verified
three-kernel/28-print count is unchanged. A separate [arithmetic20 audit](review_compact_arithmetic20_field_closure_batch_guard_20261020.md)
found no missing field-reduction proof interface or dependence on parked19.

## Final certificate wording and fresh compact audits

The current finite57 source is `c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`
(156120 bytes,3322 lines). A final saved-table wording correction requires
a tight predecessor at every positive **stored finite** value and limits
the predecessor criterion to positive ranks, with zero as the reachable
base. The exact two-hunk delta against the archived1d8332ec source has
[separate source-review acceptance](review_residual_bellman_certificate_source_batch_guard_20261020.md).
The recurrence, data and both earlier inserted proof fragments are unchanged.
The corresponding Bellman kernel remains uncompiled. The44-page89bd PDF
remains a predecessor; the current-source PDF has not yet been built.

The latest45-second admission watch sampled4.00–4.93GiB available RAM
and launched no TeX child. Its observations belong to the preceding1d8
source and are retained in the bridge receipt. The6GiB gate is unchanged.

Fresh direct audits also find no must-fix gap in the
[coupledD23/dimension24 argument](review_compact_d23_coupled_dimension24_adversarial_batch_guard_20261020.md)
and the [compactRA4 volume constraints](review_compact_ra4_volume_constraints_fresh_adversarial_resource_20261020.md).
They reconstruct the actual proofs and disclose prior contributor roles;
they are internal written audits, not external refereeing or new executions.
Two older verification guides now distinguish the actual compact23
PDF/conditional accounting-kernel passes from their earlier preparation
status. No dimension bound or example changed.

A bounded build of the currentc300 source was subsequently attempted and
deferred at5.2GiB available RAM, before any TeX child. This newer
source-bound attempt is in the bridge receipt; older3.7GiB attempts
remain historical. No current-source PDF was produced.

## Current PDFs completed after the resource gate reopened

Both revised papers now compile under the unchanged serial59-second/1GiB
limits. Current pages are **48/17/36/5/17/14**, in the table order.
Finite57 sourcec300 has PDF`1f3c020c…bee926`; compact23 source38c8290b
has PDF`2930f6e0…d1cb1`. All six current source/PDF/build-ledger bindings
match in the [updated PDF receipt](six_paper_pdf_completion_20261020.json).
The successful-batch count is now12, including earlier typography fixes.
Five PDFs have zero warnings; finite57 has21 underfull spacing warnings
and no overflow or unresolved-reference warning. Its actual48-page count
comes from the PDF page tree because the build-log page regex returned
null. This does not affect the successful TeX build.

Root visually inspected revised finite57 pages5,16,24,25: the saved-table
statement/proof, affine triangle argument, and coarse coefficients/formulas
fit and are legible. The [visual receipt](finite57_publication_delta_visual_review_20261020.json)
records this limited inspection and removal of all four temporary PNGs
(1,184,571 bytes). This is not inspection of every page or mathematical
verification. TeX intermediates were removed by the builder.

The compact23 source now says correlation **at most**1/3; the endpoint
can equal the coarse upper bound, while the strict normalized demands
still contradict Cauchy by exact margin43/1875. The one-phrase delta and
byte-preserved predecessor are bound in the
[correlation receipt](compact23_correlation_precision_20261020.json).
The [fresh classical/forest audit](review_compact_farthest23_classical_forest_finite57_20261020.md)
accepts the final source. Bounds, examples and verified Lean coverage are
unchanged. Earlier pending-build statements above record the chronological
pre-build stages and are superseded by this section.

Five further D13 parts now have actual bounded PASS records: support,
rank15, one-sided, quiet and literal. Their total wrapper time is6.854s.
Current coverage is10/19; later0 was deferred at the6GiB admission gate
before a mathematical child. The metadata coverage verifier independently
reports10covered and nine missing. This is not a complete D13 replay.

The final compact23 date-only cleanup removes the obsolete words
“uncompiled source draft” from its title date. Current source is
`2718192d0d34e60418f0d0c7e00c3ed211ed0cdc0613cf5978771963a4a19705`,
PDF`39e665d8…481225`, still36pages and zero warnings. Exact comparison
with the archived38c829 source confirms that all mathematical text is
unchanged. Its actual serial rebuild passed, bringing successful TeX
batches in this follow-up to13. All six current bindings still match.

The complete19part D13 replay now passes under unchanged caps:24.408s
total wrapper time and184,041,472B maximum process-tree RSS. The current
coverage validator independently confirms19/19. The primary full plan
now selects the19parts plus parser fixture in place of the eager loader;
actual plan output is63jobs includingPDF, with the9job core unchanged.
The old44job plan and failure evidence are preserved. Thirty of31
selected higher-tail leaves now pass (65.182s total,256,278,528B maximum);
only the108-range alternative replay remains incomplete.

Its first pilot exposed a SymPy-Boolean counting TypeError. The minimal
integer-count correction preserves the exact predicate and complete
domain and has independent source review. The current reader fixture
and math intervals remain unexecuted after the RAM admission gate
closed. The new metadata validator correctly rejects the old-source
metadata PASS, giving zero current108units, while retaining19/19D13.
Neither this partial package nor the new63job plan is a full runtime PASS.

The residual Bellman module's proof body is byte-identical to the reviewed
c387 version; its currentcd2a05a6 header refers to durable status evidence.
It remains uncompiled. The verified three-kernel/28-print count is unchanged.

The final documentation audit reconciled all ten current guides against
the receipts and actual63job/9job plans. It checked823local file links
and7Markdown fragments, all resolving. Root's final integrity receipt
binds all six source/PDF/build pairs,21latest byte-verified archives,
all44unchanged sharpRA frozen inputs, and the ten guide versions.
The proof lock is idle, unused proof work files are zero, and all root
TeX/scientific jobs are terminal. About171GiB disk remains free.

## 20 October 2026 follow-up: Bellman compilation and partial108 replay

The residual Bellman source `cd2a05a6…` has now actually compiled under the
unchanged one-thread/2048MiB runner:0.54s for compilation and0.41s for the
separate axiom audit. All eight printed declarations use only `propext`
and `Quot.sound`, with no `sorryAx`; the compiled object is `d5c0b98f…`.
The [current logic receipt](publication_logic_verification_20261020.json)
therefore records **four selected kernels and36 actual axiom prints**.
This supersedes the historical uncompiled/three-kernel statements above;
the Lean proof body did not change. The
[assumptions map](../docs/residual_bellman_certificate_logic_20261020.md)
still requires a concrete complete item alphabet, integer weights,
interaction-free context separation, actual memo-certificate validation
and geometric domination. The abstract theorem proves the certificate
implication, not those application inputs. Finite57's
[current metadata plan](finite57_bellman_verification_plan_20261020.json)
now has10 core jobs and64 full jobs includingPDF; the earlier9/63 plan
is preserved. A plan is not an actual combined runtime PASS.

The current distance9 worker is `f5f0b470…`. Its three expected-key
tuple/list corrections are representation-only and have a
[complete source/type inventory](review_distance9_expected_key_json_representation_fix_finite57_20261020.md);
the earlier failed attempts and minimal Boolean-count repair remain
preserved. The [current scanner fixture](distance9_json_key_lists_scanner_fixture_execution_20261020.json)
(`f6fd1267…`) has actually passed, checking13 positive and32 negative
fixtures plus the1262-row frozen projection. The current `metadata [0,1)`
and `domains [0,1)` phases have also passed. The full `profiles [0,8)` pilot passed,
recomputing all756 scalar profiles across link dimensions108–115.
These actual outcomes are retained in the
[partial runtime ledger](finite57_distance9_partitioned_replay_runtime_20261020.json).

Some overlap-domain intervals have also passed, but the complete
distance9/108-range conjunction remains incomplete: remaining overlap-domain,
join/inertia, overlap and isolation units are not credited by a fixture,
a metadata self-test or a scalar-profile pass. Resource-gate deferrals
launch no mathematical child and add no coverage. The current
[coverage validator](partitioned_replay_coverage_validation_20261020.json)
(`456abc2b…`) and reviewed resume coordinator preserve old failed/stale
rows and accept only valid current-source receipts. No full-manuscript
runtime PASS or cap increase is claimed. All six current PDFs and their
matching source/build bindings remain unchanged by this follow-up.


## Projection and height logic audit continuation

The six manuscript sources and PDFs remain unchanged. The new coordinate
projection/height kernel compiled on its first attempt (3.47 seconds) and
passed all eight printed-axiom checks (3.15 seconds), bringing the selected
publication logic to five kernels and 44 actual prints. Only the standard
axioms occur. Actual face incidence and real geometric interpretation remain
explicit inputs; this is not an end-to-end geometric Lean proof. The
[current receipt](publication_logic_verification_20261020.json) and
[assumptions map](../docs/right_angled_projection_height_logic_20261020.md)
record the precise boundary.

The inherited108 replay now has complete generation, disjoint-join, inertia
and scalar-profile phases, plus97 of384 overlap classifications. The
heaviest owner passed in35.723 seconds at139,264,000 bytes peak RSS. All949
isolation cases remain pending; its first pilot was stopped at admission,
without starting a child. The combined receipt check remains incomplete
(19/19 D13 parts and960/2196 distance9 units). The
[current combined plan](publication_projection_height_verification_plan_20261020.json)
is metadata only, not a full runtime pass. Resource gates are unchanged.

## Completed inherited replay and current sharpRA core

All eight current distance9 phases now have actual bounded PASS records:
2196/2196 units, including all384 overlap-owner classifications and949
isolation cases. The56 successful partitions took388.008 seconds in total,
with159,432,704 B maximum process-tree RSS. The
[terminal runtime ledger](finite57_distance9_partitioned_replay_runtime_20261020.json)
preserves every earlier failed, stale and deferred attempt; its original
checkpoint is retained as history. The combined
[receipt validation](partitioned_replay_coverage_validation_20261020.json)
now PASSes19/19 D13 parts and2196/2196 current distance9 units. This
closes the outstanding equivalent replay of the inherited108 tail; it
does not turn the earlier failed monolithic command into a success or
establish a combined full-manuscript runtime pass.

The current sharpRA selection has an actual
[six-job combined PASS](right_angled_sharp_projection_height_core_verification_20261020.json):
two exact checks, three conditional Lean roots with fresh axiom output,
and a warning-free temporary five-page PDF build. Total outer job time
was5.439 seconds. Its
[actual-output validation](right_angled_sharp_projection_height_core_axiom_validation_20261020.json)
checks all30 emitted declarations against the source printer lists and
accepts only the standard axioms. These repeat existing declarations;
the separate five-kernel/44-print total is unchanged. The published PDF
and all six manuscript sources remain unchanged.

The [selected-verifier mode audit](review_six_paper_selected_verifier_modes_20261020.md)
found no actual plan/generation bypass in the current scientific commands.
It explicitly qualifies eight legacy Lean targets as cached objects,
rather than fresh capped compilations. A separate source-only exposition
audit found the algorithm, termination/completeness argument and exact
output stated at the finite57 paper's declared finite-data boundaries.
SharpRA has no essential script-computed input: the selected scripts
check the algebra and model already proved in its text. Neither audit is
external refereeing or an end-to-end geometric formalization.

Root and independent review then accepted the
[catalogue proposal](finite57_distance9_partitioned_catalogue_proposal_20261020.json).
The [applied delta](finite57_distance9_catalogue_integration_20261020.json)
replaces only the primary finite57 selection of the old108 command with
the current scanner contract and56 whole interval commands that actually
passed. Their exact disjoint cover is2196/2196; no command was invented
by trimming or merging successful intervals. The old registry entry and
previous catalogue are byte-preserved, and every other paper selection
is unchanged. The [actual current CLI plans](publication_complete_partitioned_verification_plan_20261020.json)
have120 full/10 core finite57 jobs includingPDF, and6/6 sharpRA jobs.
Plans are metadata only; no full120-job runtime PASS is claimed.

The five-kernel receipt now binds ten immutable compilation/audit record
copies, so later runner cache refreshes cannot replace the bound evidence.
The two recently repeated audit outputs match the previous child bytes
and declaration sets. Earlier receipt/summary versions are byte-verified
compressed archives. Sources, objects and the five-kernel/44-declaration
count are unchanged.


## Independent priority-paper audits and scalar-tail verification

Fresh independent source reviews of the current
[finite57 composition](review_finite57_theorem_composition_adversarial_20261020.md)
and [sharp right-angled foundations](review_right_angled_sharp_foundations_adversarial_20261020.md)
found no internal gap at their declared geometric and finite-data input
boundaries. These are internal mathematical reviews, not external refereeing.
The six TeX sources and matching PDFs are unchanged.

Two additional inherited forest phases now have actual bounded PASS:
8114 Bellman tables and10,195,241 transition inequalities,5.106 seconds
combined,318,799,872 bytes peak process-tree RSS. The
[execution summary](finite57_ancestor_forest_replay_execution_summary_20261020.json),
[source review](review_finite57_ancestor_forest_batches_source_composition_20261020.md)
and [scope contract](finite57_ancestor_forest_replay_contract_20261020.md)
state exactly which native computations were replayed. This closes the
coarse and first increasing-family replay scopes; some ambient59–65
forest/local-census outputs remain declared finite-data inputs. A full
selected plan is not regeneration of every transitive producer.

The new [scalar-tail kernel](../formal/CoxeterProofs/FiniteVolumeTailKernel.lean)
compiled in1.18 seconds. Its nine actual printed declarations pass the
[axiom audit](FiniteVolumeTailKernel_axiom_summary_20261020.json), with only
the standard axioms and no `sorryAx`. It proves105 literal cap comparisons,
scaled payment subtraction, and the near/far scalar contradiction. The
[assumption map](../docs/finite_volume_tail_kernel_assumptions_20261020.md)
keeps the actual geometric interpretation, concrete cost tables and
compact-facet alternative outside this kernel. The
[master receipt](publication_logic_verification_20261020.json) now binds
six kernels and53 actual axiom prints through immutable execution records.
The earlier five-kernel/44-print receipt is byte-preserved.

The independently reviewed [catalogue integration](finite57_audit_completion_catalogue_integration_20261020.json)
selects the two passed phases and new kernel, and reconciles stale status
fields against existing receipts. The
[current actual CLI plans](publication_audited_forest_tail_verification_plan_20261020.json)
have11 core/123 full finite57 jobs including PDF, and6/6 sharpRA jobs.
No combined123-job runtime PASS or end-to-end geometric Lean proof is
claimed. Parked arithmetic19 refinements remain preserved and excluded;
no dimension bound, example or manuscript family was added.

A final evidence audit found that two original kernels retained parsed
axiom summaries and hashes without raw stdout files. Two short cached
[audit recaptures](publication_logic_two_raw_output_recaptures_20261020.json)
restored those files; their2050 and1085 bytes match the original immutable
child-output hashes exactly. No object was recompiled and no declaration
was added to the six-kernel/53-print count.


## Primitive shape counting and primary-source audits

The [new shape-counting kernel](../formal/CoxeterProofs/FiniteVolume58ShapeKernel.lean)
(source40d12b57) passed its first capped compilation in0.59 seconds and
its nine-declaration axiom audit in0.55 seconds. Its
[raw-output summary](FiniteVolume58ShapeKernel_axiom_summary_20261020.json)
binds all3188 child bytes, which match both compile/audit hashes. Only
standard axioms occur and there is no `sorryAx`. Eight harmless lint
warnings remain; the reviewed source was unchanged after compilation.

The new content derives the minimum-index assignment and the two841
extremal gains from primitive below/above index counts, reserves two
**distinct** positions once, and obtains the754/8410 doubled shape lower
bound from the two incidence identities. It permits signed local
curvatures and counts either finite or ideal extrema. No desired shape
inequality is an input. The [assumption map](../docs/finite_volume58_shape_kernel_assumptions_20261020.md)
keeps actual link/index interpretation, degree balance and optional Nat
normalization explicit; the legacy `WrittenGeometry.shapeLower` field
is not automatically instantiated. The independent
[source review](review_finite_volume58_shape_kernel_source_20261020.md)
and [execution/replacement review](review_finite57_shape_kernel_catalogue_replacement_and_execution_20261020.md)
pass.

The [reviewed catalogue replacement](finite57_shape_kernel_catalogue_integration_20261020.json)
selects this source in place of the primary paper's standalone older
ambient59 arithmetic module. The older source/evidence and historical
atlas selection remain intact. Current finite57 plans stay11core/123full
jobs including PDF, with sharpRA6/6 in the
[new actual metadata plans](publication_shape_counting_verification_plan_20261020.json).
No combined123-job runtime PASS is claimed. The
[master logic receipt](publication_logic_verification_20261020.json) now
binds seven kernels,62 actual axiom prints and14 read-only execution
snapshots; the previous six-kernel receipt is byte-preserved.

The [sharpRA input summary](sharp_ra_external_input_audit_summary_20261020.json)
and [Khovanskii shape audit](review_finite57_published_shape_inputs_20261020.md)
check the precise published premise scopes. Original Musin, Kolpakov,
Vinberg and Khovanskii texts were read. Dufour's full original article
was not obtained; the exact mixed-vertex27 consequence was checked in
Alexandrov's original text, as already stated by the manuscript. No
citation or vertex-type gap was found. These are internal and
bibliographic audits, not external refereeing of the new theorems.

All six source/PDF pairs and numerical dimension ceilings are unchanged.
The source-completion and initial priority-audit milestone is fulfilled
at its declared boundaries; the broader dimension/example goal remains
active, with further Lean work following settled mathematics.
