# arXiv submission editions

Updated 19 November 2026. Four prepared manuscripts belong to the existing
three paper families. One combined article joins the finite-volume, compact
and arithmetic compact dimension proofs; its LP companion recovers those
same bounds through a common feasibility and certificate framework.
The six accepted verification editions retain their existing selections.
Preparing a source folder does not submit it to arXiv.

## Sources, PDFs and results

| Paper | Complete TeX | Matching PDF | Figures | Results |
|---|---|---|---:|---|
| Sharp dimension bounds for right-angled hyperbolic polytopes | [main.tex](../arxiv/sharp_dimension_bounds_8_and_4/main.tex) | [6 pages](../arxiv/sharp_dimension_bounds_8_and_4/main.pdf) | 2 | Finite-volume right-angled dimension at most 8; entirely ideal dimension at most 4. Both bounds are attained by known examples. |
| On the second smallest ideal right-angled hyperbolic four-polytope | [main.tex](../arxiv/ideal_four_dimensional_volume_and_structure/main.tex) | [18 pages](../arxiv/ideal_four_dimensional_volume_and_structure/main.pdf) | 5 | The second ideal volume is 8π²/3, uniquely attained by the reflected 24-cell double. A separate restricted facet-spectrum theorem gives chains of 24-cells. |
| Dimension bounds for hyperbolic Coxeter polytopes | [main.tex](../arxiv/hyperbolic_coxeter_dimension_bounds/main.tex) | [41 pages](../arxiv/hyperbolic_coxeter_dimension_bounds/main.pdf) | 0 | Finite-volume dimension at most 57, compact dimension at most 17, compact arithmetic dimension at most 15. |
| Linear programming certificates for dimension bounds of hyperbolic Coxeter polytopes | [main.tex](../arxiv/linear_programming_dimension_bounds/main.tex) | [36 pages](../arxiv/linear_programming_dimension_bounds/main.pdf) | 0 | LP recovery of the existing bounds 57/17/15 from proved local inequalities and complete domains. |

The [current source/PDF manifest](../release/manifest.json)
records the exact editions. Each upload folder contains `main.tex`,
`references.bib` and the generated `main.bbl`. Figures are inline TikZ;
no separate image or custom class file is needed. All four papers name
**Alexander Kolpakov**, include an AI-use acknowledgment and cite the
repository in the bibliography. The three earlier manuscripts acknowledge
OpenAI's GPT-6-astra and GPT-6.1-sol; the LP manuscript acknowledges Codex.
For each submission, zip only `main.tex`, `references.bib` and
`main.bbl`. Exclude PDFs and temporary files produced by your TeX editor.

## The LP companion submission

The [canonical LP submission](../arxiv/linear_programming_dimension_bounds/README.md)
contains a single `main.tex`, `references.bib`, the generated `main.bbl`
and the matching PDF. The four reviewed TeX files and PDF in
[`manuscripts/linear_programming_dimension_bounds/`](../manuscripts/linear_programming_dimension_bounds/README.md)
remain the prior working edition. The conversion inlines the TeX sections
and supplies a BibTeX bibliography; it does not change the mathematical
models or rerun their checkers.

The local archive is
`dist/CoxeterPolytopes-linear-programming-dimension-bounds-arxiv.zip`,
with its `.sha256` checksum beside it. It contains exactly `main.tex`,
`references.bib` and `main.bbl` at the archive root. Build and package
receipts are under `release/linear_programming_dimension_bounds/`; the
[release manifest](../release/manifest.json) identifies the current files
and their hashes. This is another prepared manuscript in the dimension-bound
family, leaving the three earlier submission editions unchanged.

## The combined dimension article

The concise 41-page submission presents the three theorems, main geometric
arguments and exact interfaces to the [GitHub proof catalogue](../arxiv/hyperbolic_coxeter_dimension_bounds/ADDENDA.md).
Complete case and computational proofs remain in the repository. The existing
[551-page full edition](../manuscripts/addenda/dimension_bounds_20261030/main.pdf)
is preserved byte for byte as the proof addendum. Its contents are:

| Part | Argument | Begins on page |
|---|---|---:|
| I | Finite-volume dimension at most 57 | 5 |
| II | Compact dimension at most 17 | 17 |
| III | Compact arithmetic dimension at most 15 | 35 |
| IV | Finite-volume proof appendices | 72 |
| V | Compact proof appendices | 107 |
| VI | Arithmetic proof appendices | 372 |

The first two results assume no arithmeticity. The arithmetic argument uses
the independently proved compact theorem within the same document. Arithmetic
reflection lattices automatically have the quadratic-form model used there.
The shared introduction explains the methodological split: local angle and
curvature budgets for finite volume, farthest-vertex component energies for
compactness, then field restrictions in the arithmetic case.

Assembly review compared all three retained main texts and appendix groups
with their reviewed illustrated sources. Labels were made distinct, companion
references became internal references, and the bibliography was merged. The
geometric arguments, finite algorithms, completeness proofs and required outputs
remain. The [dimension-paper guide](../arxiv/DIMENSION_BOUNDS.md) gives the
underlying source and proof maps. This is source-preservation and presentation
review, not a new mathematical replay.

## The right-angled articles

The dimension paper begins with the published ceilings and known examples.
An outline explains the shared maximum-height projection and the different
contradictions in the general and ideal cases. The complete cubic tensor
argument and two diagrams remain in the paper.

The volume paper presents a single sequence: the volume gap, the equality facet
spectrum, and the cutting and recognition theorem. Its five diagrams illustrate
the geometry and incidences. The chain theorem assumes every facet is an
octahedron or a doubled octahedron; it does not classify arbitrary ideal
four-polytopes as chains.

The historical introduction relates the second ideal volume to earlier extremal
questions. Potyagailo and Vinberg explicitly posed the compact 120-facet
question. The [Guo–Zheng preprint](https://arxiv.org/abs/2610.33053), Theorem 1.1
and Corollary 1.2, establishes the unique compact 120-cell facet and volume
minima. The earlier ideal theorem gives the 24-facet and volume 4π²/3 minima.
The present paper addresses the next distinct ideal volume and its equality case.

Our compact RA4 paper remains a paper of constraints, without a compact
second-volume classification. The [separate chronology audit](../results/arxiv_compact_prior_publication_audit_20261028.md)
records that this repository was created after the Guo–Zheng preprint submission;
no earlier public priority is claimed for related local results.

## Supporting checks and their scope

| Submission material | Supporting selection or map | Scope |
|---|---|---|
| Sharp dimensions | `right_angled_sharp`; [proof map](right_angled_sharp_dimensions_publication_proof_map_20261021.md) | Two exact checks: displayed 24-cell coordinates and universal rational tensor identities. |
| Ideal volume and structure | `ideal_ra4_volume_and_structure`; [proof map](ideal_ra4_volume_and_structure_proof_map.md) | Four exact checks of the saved double, equality models and specified root/opposition data. |
| Finite-volume 57 | [Reproduction map](finite_volume57_reproduce.md) | The 259-unit current-rule conjunction and earlier bounded replacements retain their explicit premises; no combined 146-job replay is claimed. |
| Compact 17 | [Complete proof map](compact17_consolidated_proof_map_20261025.md) | Full written case proofs, with selected conditional scalar/composition kernels. |
| Compact arithmetic 15 | [Dependency map](compact_arithmetic_consolidation_proof_map_20261025.md) | Four computed premises, the compact theorem and the published golden-lattice theorem remain explicit. |
| LP recovery of 57/17/15 | [LP paradigm and evidence guide](../notes/lp/UNIFORM_LP_PARADIGM.md) | Retained local theorems and complete domains feed the common LP framework. Exact certificates have their own execution records; submission typesetting adds no mathematical replay. |

The two right-angled script-only selections use small tracked inputs and require
none of the unuploaded bulk collection. The broader dimension arguments retain
the availability and execution limits stated in the [public-data guide](../repro/PUBLIC_DATA.md).
The [release evidence guide](right_angled_release_packages_20261025.md) and
[formal guide](../formal/README.md) identify the existing conditional kernels.
Geometry, literature and specified finite computations remain explicit premises.
Printed foundational axioms do not enumerate ordinary theorem parameters.

Inspect selections without launching mathematical jobs:

```sh
python3 -B scripts/run_proof_batch_20261024.py manuscript \
  right_angled_sharp --mode scripts --plan
python3 -B scripts/run_proof_batch_20261024.py manuscript \
  ideal_ra4_volume_and_structure --mode scripts --plan
```

For an actual necessary replay, consult [repro/README.md](../repro/README.md),
remove `--plan`, and supply a fresh ledger path. Execution remains serial and
resource-monitored. The catalogue's `--pdf` option typesets its preserved
5/17-page right-angled editions, not these 6/18-page sources. Historical
`compact20` and `compact_farthest_bounds` selections likewise retain their
20/23-edition scopes; they do not replay all later 15/17 appendices.

## Building the current submission PDFs

The sources are already in their separate `arxiv/` folders. Open `main.tex`
in a TeX editor configured for pdfLaTeX and BibTeX, or run this in the chosen
paper's folder:

```sh
latexmk -pdf main.tex
```

The equivalent manual sequence is:

```sh
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

BibTeX reads `references.bib` and writes `main.bbl`, which is included in the
arXiv source folder. The ordinary workflow needs a standard TeX distribution;
it does not need Python, repository data or mathematical checkers. All figures
are in `main.tex`. Upload `main.tex`, `references.bib` and `main.bbl` together,
with a separate submission for each paper.

The repository's own typesetting run uses the same pdfLaTeX/BibTeX sequence
in a temporary directory and additionally records source hashes, resource
measurements and PDF inspection. It retains the final PDF and `main.bbl`,
then removes temporary build files. It launches no mathematical checker or
Lean module.

For a revised LP submission build, run from the repository root with a fresh
tag and ledger path:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock \
  --ledger results/fresh_lp_arxiv_build_monitor.json -- \
  scripts/build_lp_dimension_arxiv_20261118.py --tag fresh_tag
```

The initial preparation tag is `submission01_20261118`. Its build and
package receipts are
`release/linear_programming_dimension_bounds/build_submission01_20261118.json`
and
`release/linear_programming_dimension_bounds/package_submission01_20261118.json`.

## Verification records

The current [manifest](../release/manifest.json)
binds the reading PDFs to `main.tex`, `references.bib` and `main.bbl`.
Its build receipts record the pdfLaTeX and BibTeX stages, bibliography entry
order, warnings, references and PDF fonts. Source review and sampled layout
inspection have separate scopes; no claim of visual inspection of every page
of the long article is made.

The three earlier four-step builds passed with zero TeX or BibTeX warnings,
no undefined references or citations, and Type 1 fonts. Their PDFs have
6, 18 and 41 pages; the LP submission adds a 36-page PDF with separate
receipts below. The full proof addendum retains its 551 pages. The builds' temporary
products were removed, while `main.bbl` is retained for each submission.

The current build receipts are:

- [Sharp 6-page build](../release/right_angled_sharp/build_20261104.json).
- [Ideal 18-page build](../results/arxiv_ideal_bibtex_build_20261030.json).
- [Combined 41-page build](../release/dimension_bounds/build_ideas_final_20261104.json).
- [LP 36-page build](../release/linear_programming_dimension_bounds/build_language03_20261119.json)
  and [source-package receipt](../release/linear_programming_dimension_bounds/package_language03_20261119.json).
- [Preserved 551-page addendum build](../results/arxiv_dimensions_bibtex_build_20261030.json).

No successful mathematical worker or Lean module was rerun for preparation of
these submission editions. Existing executions retain their original sources
and declared scopes.

The 19 November language revision checks terminology against earlier work by
Kolpakov, Vinberg and Nikulin. The current LP article has 36 pages. Its
mathematical preservation review is
[`lp_manuscript_language_source_review_20261119.json`](../audit/lp_manuscript_language_source_review_20261119.json);
the current build and package use tag `language03_20261119`.
