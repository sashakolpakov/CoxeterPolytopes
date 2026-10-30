# Current bounds and verification

Updated 30 October 2026. The repository has three current submission manuscripts.
Their results are internally reviewed written proofs at the stated geometric,
literature and finite-data inputs. No new polytope example is claimed.

## Dimension bounds

| Class | Largest published example | Published ceiling | Repository ceiling |
|---|---:|---:|---:|
| General Coxeter, finite volume | 21 | 995, Prokhorov/Khovanskii | **57** |
| General compact Coxeter | 8 | 29, Vinberg | **17** |
| Compact arithmetic Coxeter | 8 | 29, Vinberg | **15** |
| Noncompact arithmetic Coxeter | 21 | 21, Esselmann | **21**, published and sharp |
| Arithmetic Coxeter, compact or noncompact | 21 | 29 overall | **21**, sharp overall |
| Right-angled, finite volume, all vertex types | 8 | 12, Dufour | **8**, internally sharp |
| Right-angled, entirely ideal | 4 | 6, Kolpakov | **4**, internally sharp |
| Right-angled, compact | 4 | 4, published and sharp | **4**, unchanged |

The general finite-volume, general compact and right-angled arguments assume
no arithmeticity. Arithmetic reflection lattices automatically have the
quadratic-form model; “first type” adds no restriction to these arithmetic rows.
Compact and noncompact arithmetic cases have different proofs. Exact published
scopes and references are in the [assumption register](manuscript_assumptions.md).

## Current manuscripts

| Manuscript | PDF | TeX source |
|---|---|---|
| Dimension bounds 57 / 17 / 15 | [551 pages, 9 figures](../manuscripts/build/arxiv_hyperbolic_coxeter_dimension_bounds_20261028_current_20261030.pdf) | [main.tex](../arxiv/hyperbolic_coxeter_dimension_bounds/main.tex) |
| Sharp right-angled dimensions 8 / 4 | [6 pages, 2 figures](../manuscripts/build/arxiv_sharp_dimension_bounds_8_and_4_20261028_current_20261030.pdf) | [main.tex](../arxiv/sharp_dimension_bounds_8_and_4/main.tex) |
| On the second smallest ideal right-angled hyperbolic four-polytope | [18 pages, 5 figures](../manuscripts/build/arxiv_ideal_four_dimensional_volume_and_structure_20261028_current_20261030.pdf) | [main.tex](../arxiv/ideal_four_dimensional_volume_and_structure/main.tex) |

The combined dimension article has 71 pages of main exposition and complete
proof appendices from page 72. Parts I/IV prove finite volume 57; Parts II/V
prove compact 17; Parts III/VI prove arithmetic compact 15, using the compact
theorem within the same document. The appendices supply necessary case proofs.
The [contents guide](../arxiv/DIMENSION_BOUNDS.md) identifies their locations.

Each source folder contains `main.tex`, `references.bib` and the generated
`main.bbl`, with figures inline. Standard pdfLaTeX and BibTeX, or
`latexmk -pdf main.tex`, build the paper without Python. The
[current manifest](../results/arxiv_submission_current_manifest_20261030.json)
binds the exact sources, bibliographies, PDFs and build receipts. Compilation and source-assembly
review add no new mathematical or Lean execution.
All three four-step pdfLaTeX/BibTeX builds passed with zero TeX or BibTeX
warnings, no undefined references or citations, and Type 1 fonts.

## Verification boundaries

“Internally reviewed” means a complete written argument, a source-bound review
by another author and root acceptance. It does not mean external refereeing,
complete computational regeneration or a formalization of the whole geometry.

| Result | Written and computational evidence | Conditional Lean scope and remaining inputs |
|---|---|---|
| Finite-volume dimension 57 | Complete argument; [259-unit current-rule conjunction](../results/finite57_current_rule_release_completion_20261025.json), including all 64 selected application shards, covers 28,028 maps, 12,518 nonunit applications and 83 prefixes. [Reproduction and proof map](finite_volume57_reproduce.md). | Selected kernels check internal deductions. Geometric statements, native classifications and specified inherited finite data remain inputs. This is not a fresh combined replay of the historical 146-job inventory or full transitive regeneration. |
| Compact dimension 17 | Complete geometric case proofs and global composition; [proof map](compact17_consolidated_proof_map_20261025.md). | [Two modules compiled and separately axiom-audited](../results/compact17_targeted_lean_completion_20261025/manifest.json): seven scalar and three composition declarations. Geometry, analytic estimates, coefficient funding and complete inventories remain explicit premises. |
| Compact arithmetic dimension 15 | Complete field reduction, four computed propositions, final dyadic exclusion and global composition; [dependency map](compact_arithmetic_consolidation_proof_map_20261025.md). | [Fourteen dyadic declarations compiled and separately audited](../results/lean_quadratic_dyadic_triple_verification_20261023/manifest.json). Number fields, geometric realization and the named computed outputs remain external inputs. |
| Right-angled dimensions 8 / 4 | Complete projection, tensor and ideal-equality arguments; two exact script checks. [Proof map](right_angled_sharp_dimensions_publication_proof_map_20261021.md). | Selected projection/height/tensor/count kernels, with 30 audited declarations and a separately recorded five-declaration regular-cell supplement. Actual geometry and the cited code, facet and rigidity theorems remain inputs. |
| Second ideal four-volume and unique double | Complete gap, equality-spectrum and recognition arguments; four exact model/incidence/root checks. [Proof map](ideal_ra4_volume_and_structure_proof_map.md). | The selected ideal kernel has eight audited declarations. Actual cuts, topology, recognition, gluing and published classifications remain written or literature inputs. The chain theorem has its stated facet-spectrum hypothesis. |

Finite57's 25 bounded replacement PASS units and the other 111 inherited
script-command records retain their original scopes; the core eleven is
unchanged. The [current-rule reproduction guide](finite_volume57_reproduce.md)
distinguishes these records from the 259-unit application conjunction. A source
hash, a saved output, a written review and an actual checker execution establish
different things. Failed and deferred attempts remain in their dated records.

The arithmetic proof uses the published golden-lattice theorem at its literal
lattice, full reflection-group and dimension scope. Its
[primary-source authentication](bugaenko_golden_lattice_cutoff_authentication_20261026.md)
is complete; that theorem remains an external mathematical input. The full
invariant-lattice argument and four computed propositions remain necessary.

Lean verifies the internal implications from the hypotheses listed in each
logic map. Printed foundational axiom lists do not enumerate ordinary theorem
parameters. The [formal guide](../formal/README.md) identifies the separate
compilation and actual axiom-audit records; no complete geometric Lean proof
is claimed.

## Volume results and reproduction

The ideal-volume article proves that the second distinct volume is 8π²/3,
uniquely attained by the reflected double of the regular ideal 24-cell. Its
separate chain theorem applies to octahedral and doubled-octahedral facets;
it is not a decomposition theorem for arbitrary ideal four-polytopes.

The [compact RA4 source](../manuscripts/compact_ra4_volume_constraints.tex)
and [proof map](compact_ra4_volume_constraints_reproduce_20261019.md) provide
supporting constraints. The [Guo–Zheng preprint](https://arxiv.org/abs/2610.33053)
establishes the compact 120-cell facet and volume minima. A compact
second-volume classification remains open here.

The [submission guide](arxiv_submission_editions_20261028.md),
[reproduction entry](../repro/README.md) and
[computational proof contract](computational_proof_contract.md) explain the
algorithms, required outputs and serial execution commands. The
[public-data guide](../repro/PUBLIC_DATA.md) states which inputs are available;
no asset from the separate 936-file/33-archive bulk collection has been uploaded.
Preparing the current sources and documentation did not rerun mathematical
workers or Lean modules.
