# Manuscript verification

The [three submissions](arxiv_submission_editions_20261028.md) have complete
written arguments at their stated external and finite-data premises.
The [result-to-manuscript map](result_manuscript_correspondence.md) identifies
their proof homes. Written proof, source review, exact computation, conditional
Lean verification and typesetting are separate kinds of evidence.

## Evidence for the results

| Result | Recorded evidence | Remaining boundary |
|---|---|---|
| General finite volume, dimension ≤57 | [Current-rule completion](../results/finite57_current_rule_release_completion_20261025.json): 259 scoped units, all 64 selected shards, 28,028 maps, 12,518 nonunit applications and 83 prefixes. The [proof/reproduction guide](finite_volume57_reproduce.md) connects the rules, finite data and geometric arguments. | Inherited classification and numerical inputs remain explicit. The separate application conjunction is not a combined fresh 146-job replay or complete transitive regeneration. |
| General compact, dimension ≤17 | [Complete written acceptance](../results/compact_dimension17_global_written_acceptance_20261024.json), [proof map](compact17_consolidated_proof_map_20261025.md), [targeted Lean completion](../results/compact17_targeted_lean_completion_20261025/manifest.json). | Geometry, analytic estimates and complete necessary-cover statements remain premises of the selected Lean kernels. |
| Compact arithmetic, dimension ≤15 | [Global composition acceptance](../results/review_arithmetic16_whole_d16_manual_dyadic_exclusion_and_global15_root_acceptance_20261023.md), [proof map](compact_arithmetic_consolidation_proof_map_20261025.md), [dyadic-kernel compilation and audit](../results/lean_quadratic_dyadic_triple_verification_20261023/manifest.json). | Four computed propositions, field/geometry inputs and the [authenticated external lattice theorem](bugaenko_golden_lattice_cutoff_authentication_20261026.md) remain explicit. The whole-D16 finish is a written argument. |
| Sharp right-angled dimensions and ideal RA4 second volume | [Package acceptance](../results/right_angled_current_release_root_acceptance_20261025.json), [dimension proof map](right_angled_sharp_dimensions_publication_proof_map_20261021.md), [volume proof map](ideal_ra4_volume_and_structure_proof_map.md). | Exact coordinate, tensor and finite-model checks support the written proofs. Geometric realization and published classifications remain explicit inputs. Earlier combined executions retain their original source scopes. |

## What Lean verifies

Selected Lean modules prove internal implications with geometric, literature,
field and named computation statements supplied as parameters. A successful
compilation is recorded separately from an actual axiom audit. Printed
foundational axioms do not enumerate ordinary theorem parameters.

The [formal guide](../formal/README.md) identifies the modules and replay
procedure. The [verification register](formal_verification_records_20261025.md)
and each logic map record actual compilations, separate declaration audits,
legacy cache provenance and failures. The
[selected two-module completion](../results/lean_publication_completion_20261021/manifest.json)
records two actual compilations and two separate audits, covering nine
declarations; its earlier numerical conclusion scopes do not expand to later
bounds. No end-to-end geometric formalization is claimed.

## Plans and scoped replay

Inspect the catalogue without launching proof jobs:

```sh
python3 -B scripts/verify_manuscript.py --list
python3 -B scripts/verify_manuscript.py finite57 --mode core --plan
```

Use the [reproduction guide](../repro/README.md) for a necessary replay.
The catalogue's `scripts`, `logic`, `core` and `full` selections retain their
declared input and execution scopes. Its optional `--pdf` checks the selected
source in a temporary directory; the current arXiv sources have separate
typesetting entries. A PDF build does not verify mathematics.

Scientific, Lean and TeX jobs run serially under the shared lock and measured
percentage resource gates. Do not rerun successful unchanged checks for a
documentation update. Preserve failed/refused attempts and exact receipts;
never count a plan or matching hash as an executed proof check.
