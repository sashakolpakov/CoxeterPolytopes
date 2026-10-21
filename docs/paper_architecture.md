# Six papers, with separation by proof method

21 October 2026. There are **six current papers**, complete at their stated input boundaries, with matching PDFs. The arithmetic20 paper now supplies its two formerly Markdown-only nongolden field proofs and a shared six-root figure: its current PDF is21 pages with zero warnings. Its17-page predecessor is preserved. The other five editions are unchanged. The 33 TeX files under `manuscripts/` include their predecessors,
supporting notes, three generated atlases and input fragments. A third-party
TeX source under `data/external/` is outside that count. The
[catalogue](../manuscripts/catalogue.json) records the individual roles.

The papers should each answer one mathematical question. Shared notation
and a short explanation of a common method belong in each paper that uses
them. A new result obtained by strengthening a constant in the same argument
belongs in the existing paper, rather than in a successor publication.

| Current paper | Common mechanism | Where its cases separate | Present exposition |
| --- | --- | --- | --- |
| [Finite-volume dimension 57](../manuscripts/finite_volume57.tex) | Pay for triangles and quadrilaterals at their corners, bound vertex charges, then contradict global incidence curvature. | Finite spherical vertex forests versus affine cusp factors; dimension58 and the complete higher-dimensional tail. | Complete source at its declared finite-input boundary, with the earlier [full audit](../results/review_finite_volume57_full_resource_20261020.md), [classical bridge repair](../results/finite57_full_budget_bridge_repair_20261020.json) and accepted [triangle/coarse-screen bridges](../results/finite57_publication_bridge_completion_20261020.json). Current [48-page PDF](../manuscripts/build/finite_volume57.pdf) matches its source: main exposition pp.1–12 and essential appendices/table documentation pp.13–48. The8-,41- and44-page editions are preserved. |
| [Compact arithmetic dimension 20](../manuscripts/compact_arithmetic_dimension_bound.tex) | The same angle-payment principle, with arithmetic restrictions and a strict compact counting threshold. | First the actual-ground-field reduction; then the remaining golden-field branch, supporting-star coverage, incidence transfers and the exact forest maximum. | Universal-only source complete; the [exposition audit](../results/review_compact_arithmetic20_exposition_finite57_20261020.md) and [field-closure audit](../results/review_compact_arithmetic20_field_closure_batch_guard_20261020.md) pass at their declared written-proof/finite-data boundaries. No parked19 input or new execution is used. [21-page current PDF](../manuscripts/build/compact_arithmetic_dimension_bound.pdf), with the two nongolden field arguments and shared six-root figure now proved here; the17-page and combined20-page predecessors and parked appendices remain preserved. |
| [Compact dimension 23](../manuscripts/compact_farthest_vertex_bounds.tex) | Choose one farthest vertex. All orthogonal component projections spend the squared norm of the same spatial vector, which is less than one. | Ordinary components first; large A components use coupled opposite-wall pressures, while B/D components use coordinate prefixes. Two terminal forest patterns need stronger A20 and coupled D23/A1 arguments. | Complete source, two vector figures and [fresh full audit PASS](../results/review_compact_farthest23_full_batch_guard_20261020.md), after the middle-index repair. [36-page PDF](../manuscripts/build/compact_farthest_vertex_bounds.pdf) now compiled. |
| [Sharp right-angled dimensions](../manuscripts/right_angled_sharp_dimension_bounds_publication.tex) | One finite-point projection and a height comparison. | The universal bound 8 uses a cubic spherical-code estimate; the entirely ideal bound 4 uses its own explicitly identified published inputs. | Current5-page PDF has a fresh zero-warning build after terminology, face-inheritance and height-notation clarifications; the frozen predecessor and conditional Lean evidence are preserved. A [fresh full mathematical source audit](../results/review_right_angled_sharp_dimension_bounds_full_resource_20261020.md) passes. Keep these two short branches together. |
| [Ideal right-angled four-volume and structure](../manuscripts/ideal_ra4_volume_and_structure_publication.tex) | Cusp/defect counting, followed by actual cutting and recognition. | The volume gap and equality-spectrum reduction precede the restricted Oct/D chain theorem. Only then can that theorem identify the second-volume double. | Complete publication successor and [fresh full source audit PASS](../results/review_ideal_ra4_publication_full_finite57_20261020.md). Exact finite-model algorithms and completeness arguments added, with the same theorem scope. [17-page current PDF](../manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf), with the verified16-page predecessor preserved. |
| [Compact right-angled four-volume constraints](../manuscripts/compact_ra4_volume_constraints.tex) | Perpendicular-facet capacity and neighboring excess. | Zero-excess recognition versus positive-excess and facet-link configurations. Compact metric facets require their own arguments. | Complete14-page source/PDF; later C26/C28 supporting refinements are not silently included. Minimum and second volume remain open. |

The two angle-payment papers form a methodological pair. They should each
explain the same short local-to-global principle, then mark where compact
arithmeticity or finite-volume cusp geometry changes the proof. We should
not require readers of either paper to reconstruct that principle from
an earlier research report. Their substantially different geometric
reductions justify two papers for now. The compact farthest-vertex method
is a different proof and deserves its own paper.

The current written compact bounds are now20 in general and18 in the
arithmetic first-type case over all admissible fields. The
[general20 composition](../results/compact20_and_all_field_arithmetic19_global_composition_candidate_20261021.md)
and [arithmetic18 composition](../results/compact_arithmetic_dimension18_global_composition_20261021.md),
with their full root/outside-author reviews, do not change the source/PDF stages
listed above:23/general and20/arithmetic. Integrate into the SAME two
paper families. The prior [compact21 integration plan](../results/compact21_existing_manuscript_integration_plan_20261021.md)
remains useful historical source mapping, but its terminal target is
superseded; do not append separate chapters for each intermediate bound.

The shared new approach begins with one actual farthest vector and
orthogonal component energies. State actual-point construction, positive
dual interpolation and the summed-foot inequality once, with explicit
scope. The general proof then separates into ordinary components,
large-A paired pressures and B/D prefix/terminal pressures. Its final
rank ledger must include every spectator debit and the complete whole
D21 closure. The arithmetic proof adds a visibly marked hypothesis at
the point where cyclic integrality and nonphysical conjugate positivity
restrict the possible caps. Norm/trace bounds and time-line square
classes then give the stronger packing and A16/A17 results, followed by
whole D20 and D19 exclusions. The D19 actual-field reduction is
rank-specific; it supplies no quadraticity statement at rank18.
A shared geometric theorem may be cited with a precise companion-proof
reference; it must not silently become an external established fact.

The arithmetic source's existing angle-payment20 proof is a different
method. Its complete reviewed edition must be preserved when preparing
the stronger18 successor. Do not retain its full old proof and then
append an unrelated long research dump. The new arithmetic exposition
should explain the first-type criterion, cap inventory, packing and
same-vector final rank argument in that order. All small algebraic
lists need stated completeness proofs; the new handwritten rational
and quadratic-field tables are not unexplained script outputs. Current
Lean accounting for23/20 must never be advertised as verification of20/18.

## Order of exposition

Each paper should begin with its main theorem, examples and comparison
with the published result, followed by a proof roadmap. The roadmap must
name the common reduction and the exact proposition after which separate
cases begin. It should tell the reader why each branch is needed before
introducing its notation.

Next come the common geometric lemmas, then the case-specific arguments.
Computations appear when the proof reaches a genuinely finite obligation.
Each such obligation has the mathematical statement and algorithm specified
in the [computational proof contract](computational_proof_contract.md).
The final composition then cites those statements, rather than filenames
or unexplained execution counts.

The model here is the composition of Kedlaya, Kolpakov, Poonen and
Rubinstein's [*Space vectors forming rational angles*](https://math.mit.edu/~poonen/papers/space_vectors.pdf):
the strategy appears early (§1.3), Algorithm 7.5 is followed by a
termination/correctness theorem, and the finite computation in §8 separates
exact acceptance from the proof that the filters lose no solutions.
These are editorial lessons; our mathematical reductions are different.

## What belongs outside the main paper

Keep raw rows, source hashes, execution logs and historical failed attempts
in the reproducibility package. A mathematical supplement may contain long
classification arguments, but it needs its own definitions, propositions
and proofs. A chronological collection of research notes is not such a
supplement. The 1292-, 1092- and 79-page generated atlases are data archives,
not proposed reading assignments.

The nine-field dimension-19 arithmetic refinement remains preserved.
The user has excluded it from the current universal dimension-20 paper. Its existing verified 20-page combined version remains intact as a preserved predecessor, together with both appendices and their source-bound verification records.

The exceptional compact rank-29 triangle note is supporting history,
not a current paper. Likewise, dimension-24/26/27 predecessors support
the single compact farthest-vertex paper. Earlier right-angled dimension
10/11 and ideal dimension 5/6 notes are preserved without further expansion.

## Source completion and the audit phase

All six current sources now contain their intended arguments, at their explicitly declared finite-data and inherited-theorem boundaries, and all six have internal source-audit passes. Finite57's first repair made the full-budget classical unit-charge argument explicit before the first-failing-factor split; the residual 3/2 and 4/3 arguments alone did not supply that bridge. The latest [bridge receipt](../results/finite57_publication_bridge_completion_20261020.json) records two further accepted insertions: current triangle exclusions and the complete coarse forest-screen domination. Their independent written reviews add no dimension bound or finite execution claim. The 41-page pre-bridge edition remains archived.

The six preserved pre-translation source/PDF pairs are bound by the
[audit-edition receipt](../results/six_papers_audit_edition_20261021.json),
with pages48/17/36/5/17/14. The arithmetic successor now has21 pages and
its own [zero-warning build](../results/compact_arithmetic20_internal_field_completion_final_build_20261021.json); the other five pairs are unchanged. Finite57 is unchanged at source `c300db69…`
and PDF `1f3c020c…`: 12 main pages plus36 essential proof appendices,
with21 underfull spacing warnings only. Its earlier editions survive.
The other five PDFs have zero warnings. Compact23's current
`530d964e…` / `1e6b7155…` pair replaces stale workflow prose with precise
verification scope. SharpRA's `0c29a786…` / `609a9fc0…` pair adds two short
explanations. Both changes have written delta reviews and source-bound
builds. Their predecessor bytes are preserved. Earlier mathematical,
script and selected Lean records retain their original input scopes;
these editorial builds do not add formal geometric coverage or a new
combined replay. Selected earlier figure inspections are not all-page
visual audits. Completion at a stated proof boundary is distinct from
external refereeing or submission readiness.

Replay status is separate:30 original inherited finite57 higher-tail commands pass (65.182s,256,278,528 B maximum RSS), and the remaining tail108 obligation now has a complete56-batch replacement replay. Its old59-second timeout remains preserved. All19 [D13 parts](finite_volume57_lemma_map.md#memory-bounded-d13-replay-contract) PASS (24.408s,184,041,472 B maximum RSS). The original44-job full run failed at the old eager loader; its history does not become a combined runtime PASS when the replacement parts succeed.

After the [distance9 count repair](../results/distance9_symbolic_boolean_count_repair_20261020.json), worker126d5fe3 passed its scanner fixture and metadata, then failed domains on expected tuple/list representation. The [three-comparison correction](../results/review_distance9_expected_key_json_representation_fix_finite57_20261020.md) retains all mathematical domains. Current workerf5f0b470 and fixturef6fd1267 have an actual [scanner-fixture PASS](../results/distance9_json_key_lists_scanner_fixture_execution_20261020.json) (1.732s,89,915,392 B), separate from mathematical interval coverage. All eight current phases now PASS, including384 overlap owners and949 isolation cases:2196/2196 units in56 capped batches,388.008s total and159,432,704 B peak RSS. The [receipt validator](../results/partitioned_replay_coverage_validation_20261020.json) now confirms complete combined coverage19/19 D13 parts and2196/2196 current distance9 units without rerunning mathematics. Old-source passes earn no credit. Resource caps are unchanged.

The [current verification plan](../results/publication_shape_counting_verification_plan_20261020.json)
has123 finite57 full jobs (113 script, nine Lean, one PDF), with11 core jobs. The
[reviewed catalogue change](../results/finite57_distance9_catalogue_integration_20261020.json)
adds the scanner contract and56 actual whole-interval commands in place
of the old108 selection, preserving its registry entry. The sharpRA
selection is unchanged at six jobs and now has an actual
[combined execution PASS](../results/right_angled_sharp_projection_height_core_verification_20261020.json).
These plans and scoped receipts do not claim end-to-end geometric Lean.
The two added [ancestor-forest batches](../results/finite57_ancestor_forest_replay_execution_summary_20261020.json)
actually pass:8,114 tables and10,195,241 Bellman inequalities in5.106s,
with318,799,872 B maximum RSS. Full is the selected scope; some ambient59–65
forest and local-census inputs are not transitively regenerated. No combined
finite57 full-runtime PASS follows.

The seven compiled logic kernels now have62 actual axiom prints.
The [residual Bellman certificate kernel](../formal/CoxeterProofs/ResidualBellmanCertificate.lean)
(`cd2a05a6…`) has capped compilation and eight actual axiom prints, with
only `propext` and `Quot.sound`. Its [assumptions map](residual_bellman_certificate_logic_20261020.md)
separates the abstract integer recurrence from geometry and concrete data,
and describes the local certificate conditions; concrete alphabet/context
coverage and actual table validation remain inputs. This is an abstract
exactness implication, not a new table replay or geometric Lean proof.
The [projection/height kernel](../formal/CoxeterProofs/RightAngledProjectionHeightKernel.lean)
(`2fac8028…`) adds eight actual prints for coordinate projection and strict
wall signs, maximal squared-height recurrence, finite descent and the
scalar radial gap. Actual facet/incidence, future-sheet interpretation,
face inheritance and radial realization remain inputs. SharpRA now
selects three Lean roots and six core jobs with PDF, with the combined
execution PASS linked above.
The [finite-volume tail kernel](../formal/CoxeterProofs/FiniteVolumeTailKernel.lean)
(`b8f03fd4…`) now has actual compile and nine-print PASS records. Its
[assumption map](finite_volume_tail_kernel_assumptions_20261020.md) separates
the105 scalar cap comparisons, payment subtraction and narrow-wall arithmetic
from their geometric and finite-data bridges. The new
[shape kernel](../formal/CoxeterProofs/FiniteVolume58ShapeKernel.lean)
(`40d12b57…`) has0.59s compile and0.55s nine-print audit PASS records.
It derives the uncut-edge/index bound, two841 extrema contributions
reserved once, and the finite754/global8410 doubled shape inequality.
Its [map](finite_volume58_shape_kernel_assumptions_20261020.md) retains
primitive incidence and interpretation inputs; it does not instantiate
the legacy natural-number `shapeLower` field. Replacing the primary old
ambient59 shape slot preserves job counts and all historical evidence.
The seven-kernel62-print split is10/10/8/8/8/9/9 for
cubic/concentration/compact23/Bellman/projection/tail/shape.
Core sizes with PDF are11/10/2/6/7/1 for57/arithmetic20/compact23/sharpRA/
ideal4/compactRA4.

The additional [sharpRA adversarial review](../results/review_right_angled_sharp_fresh_adversarial_finite57_20261020.md) rechecks the whole sharp-dimension argument at its external-input boundary. The [idealRA4 structural review](../results/review_ideal_ra4_structure_batch_guard_20261020.md) focuses on equality forcing the metric facet spectrum and on the actual section/cutting/chain argument. Both are written reviews, not new computations or external referee reports.
Independent [sharpRA foundations](../results/review_right_angled_sharp_foundations_adversarial_20261020.md)
and [finite57 composition](../results/review_finite57_theorem_composition_adversarial_20261020.md)
audits also support the proofs under their explicit inputs.
Primary-source checks cover [Musin/Kolpakov](../results/review_sharp_ra_published_inputs_ideal_branch_20261020.md),
[Vinberg and Alexandrov's mixed27 statement](../results/review_sharp_ra_published_inputs_dufour_20261020.md),
and [Khovanskii's shape identities](../results/review_finite57_published_shape_inputs_20261020.md).
Dufour's original full text was not obtained; the paper's stated route
through Alexandrov was checked directly.

For the two angle-payment papers, strengthen the written connection from
the geometric coverage theorem to the actual enumerator. The residual
forest recurrence already has an inductive justification. The inherited
native-support and queue reductions need equally visible finite bounds,
case coverage and source-to-algorithm correspondence. The exact current
audit is recorded in the
[inventory review](../results/editorial_manuscript_inventory_and_algorithm_audit_20261020.md)
and the
[checker-domain review](../results/review_finite57_compact20_checker_domain_contract_finite57_20261020.md).
These are exposition obligations, not new executed proof checks.

Figures should explain the proof: the common vector and its projections,
an actual cut, or the configuration responsible for a case split. They
should not decorate a list of certificate files. Page budgets are guidance;
an indispensable derivation takes precedence over a target length.
