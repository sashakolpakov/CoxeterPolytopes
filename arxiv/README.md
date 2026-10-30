# arXiv submission sources

Prepare one submission from each of these folders:

| Paper | Source | Current reading PDF |
|---|---|---|
| Sharp dimension bounds 8 and 4 | [main.tex](sharp_dimension_bounds_8_and_4/main.tex) | [6-page PDF](../manuscripts/build/arxiv_sharp_dimension_bounds_8_and_4_20261028_current_20261030.pdf) |
| On the second smallest ideal right-angled hyperbolic four-polytope | [main.tex](ideal_four_dimensional_volume_and_structure/main.tex) | [18-page PDF](../manuscripts/build/arxiv_ideal_four_dimensional_volume_and_structure_20261028_current_20261030.pdf) |
| Dimension bounds for hyperbolic Coxeter polytopes: 57 / 17 / 15 | [main.tex](hyperbolic_coxeter_dimension_bounds/main.tex) | [551-page PDF](../manuscripts/build/arxiv_hyperbolic_coxeter_dimension_bounds_20261028_current_20261030.pdf) |

Each folder contains `main.tex`, `references.bib` and the generated `main.bbl`.
Figures are inline TikZ. All papers include the author, AI-use acknowledgment
and repository citation. For each submission, zip only `main.tex`,
`references.bib` and `main.bbl`. Exclude PDFs and temporary files created by
your TeX editor.

Open `main.tex` in a TeX editor and select pdfLaTeX with BibTeX, or run
`latexmk -pdf main.tex` from that paper's folder. The equivalent manual sequence is:

```sh
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

The included `main.bbl` supplies the bibliography for the arXiv upload.
`references.bib` remains available for editing and regeneration. No Python
program is needed to typeset these papers.

The large dimension paper has a common introduction, three main arguments and
three groups of complete proof appendices. Its main exposition ends on page 71;
the appendices begin on page 72. See the [dimension-paper guide](DIMENSION_BOUNDS.md)
for the proof maps and verification records.

The verification suites are `right_angled_sharp` and
`ideal_ra4_volume_and_structure`. See the [reproduction guide](../repro/README.md)
for commands and requirements, and the [formal-verification guide](../formal/README.md)
for the conditional Lean scope. Detailed source versions and execution records
remain in the repository documentation. Each paper cites the repository once
in its bibliography.

The [current manifest](../results/arxiv_submission_current_manifest_20261030.json)
records source and bibliography files, PDFs and typesetting evidence. Bibliography
and presentation changes do not claim new mathematical-check or Lean executions.
All three papers passed the four-step pdfLaTeX/BibTeX build, with zero TeX or
BibTeX warnings and no undefined references or citations.
