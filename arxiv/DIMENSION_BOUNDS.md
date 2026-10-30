# Dimension bounds: arXiv submission

For this submission, zip only `main.tex`, `references.bib` and `main.bbl`
from `hyperbolic_coxeter_dimension_bounds/`, separately from the two
right-angled submissions. The source folder contains `main.tex`,
`references.bib` and the generated `main.bbl`. All nine figures are inline TikZ.
It compiles with pdfLaTeX, BibTeX and standard TeX packages; it needs no files
from elsewhere in the repository.

The author is Alexander Kolpakov. The acknowledgment records the use of
OpenAI's GPT-6-astra and GPT-6.1-sol and the author's responsibility for the content.
The repository is cited once in the bibliography using its main URL.

## Contents

The [551-page reading PDF](../manuscripts/build/arxiv_hyperbolic_coxeter_dimension_bounds_20261028_current_20261030.pdf)
collects three existing proof families under a common introduction.

| Part | Result or material | Begins on page |
|---|---|---:|
| Introduction | Three dimension bounds, previous work and the methodological split | 1 |
| I | Finite-volume Coxeter polytopes: dimension at most 57 | 5 |
| II | Compact Coxeter polytopes: dimension at most 17 | 17 |
| III | Compact arithmetic Coxeter polytopes: dimension at most 15 | 35 |
| IV | Complete finite-volume proof appendices | 72 |
| V | Complete compact proof appendices | 107 |
| VI | Complete arithmetic proof appendices | 372 |

The first two bounds assume no arithmeticity. The arithmetic bound uses the
independently proved compact bound. Arithmetic reflection lattices automatically
have the quadratic-form model used in that argument.

The length comes from the complete case proofs. The 71-page main exposition is
not a substitute for those appendices. The preliminary finite-volume payment
rule, the compact higher-dimensional tail and the conditioned arithmetic
chain estimates remain because the final proofs use them. Parked refinements
and exploratory attempts to improve the bounds are excluded.

## Mathematical and computational records

The preparation preserves the mathematical arguments of these retained proof
inputs:

- [Finite volume](../manuscripts/finite_volume57_consolidated.tex), including the reviewed terminal-tail lemma.
- [General compact](../manuscripts/revisions/compact_farthest_vertex_bounds_consolidated.tex), including the complete higher-dimensional tail.
- [Compact arithmetic](../manuscripts/revisions/compact_arithmetic_dimension_bound_consolidated.tex), including the four computed propositions and the direct final dyadic argument.

The detailed proof and computation maps are:

- [Finite-volume proof and reproduction map](../docs/finite_volume57_reproduce.md).
- [Compact-17 proof map](../docs/compact17_consolidated_proof_map_20261025.md).
- [Arithmetic-15 proof map](../docs/compact_arithmetic_consolidation_proof_map_20261025.md).
- [Bugaenko citation check](../docs/bugaenko_golden_lattice_cutoff_authentication_20261026.md).

The exact algorithms, completeness arguments, termination bounds and required
outputs remain in the manuscript. The arithmetic paper's former companion
citations now refer to the compact proof within the same document. The obsolete
Bugaenko source-authentication warning has been replaced by the cited published
theorem; the entire following invariant-lattice argument is retained.

Selected Lean modules verify conditional internal deductions with geometric,
literature and finite-data statements supplied as explicit premises. Preparing
this submission did not rerun mathematical enumerations or Lean modules and does
not claim a new combined mathematical replay. Existing data availability and
execution limits remain documented in the [public-data guide](../repro/PUBLIC_DATA.md).

## Typesetting and editorial checks

The typesetting pipeline uses the paper's `main.tex` and `references.bib` in a
fresh temporary directory, runs pdfLaTeX, BibTeX and two further pdfLaTeX passes,
and retains `main.bbl` alongside the source. The final PDF is kept outside the
upload folder. The manifest and build record distinguish typesetting evidence
from mathematical verification.

The four-step build passed with zero TeX or BibTeX warnings, no undefined
references or citations, and Type 1 fonts. The PDF remains 551 pages, and the
part start pages listed above are unchanged. The build's temporary products
were removed.

The assembly was independently checked against all three original main texts
and appendix groups. Changes concern the common front matter, reference
namespacing, source-provenance cuts, bibliography merge and typography.
All 667 labels and all nine figures are retained. The pinned finite57 replay
source and earlier verification receipts keep their original bindings.

The [current source/PDF manifest](../results/arxiv_submission_current_manifest_20261030.json)
identifies this edition; its current build ledger is
`results/arxiv_dimensions_bibtex_build_20261030.json`. Earlier assembly and
editorial receipts retain their original source bindings.

For ordinary typesetting, open `main.tex` in your TeX editor or run
`latexmk -pdf main.tex` from `hyperbolic_coxeter_dimension_bounds/`.
The equivalent manual sequence is:

```sh
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

No Python program is needed for this sequence. The repository's own builds
add process monitoring and isolated PDF inspection; those checks do not change
the normal TeX workflow.
