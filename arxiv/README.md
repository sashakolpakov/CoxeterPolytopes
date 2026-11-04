# arXiv submission sources

The [submission guide](SUBMISSION.md) selects the source and local package
for each of the three paper families. The concise dimension-bound submission
uses the complete earlier edition as its online proof addendum.

| Paper | Source | Current reading PDF |
|---|---|---|
| Sharp dimension bounds 8 and 4 | [main.tex](sharp_dimension_bounds_8_and_4/main.tex) | [6-page PDF](sharp_dimension_bounds_8_and_4/main.pdf) |
| On the second smallest ideal right-angled hyperbolic four-polytope | [main.tex](ideal_four_dimensional_volume_and_structure/main.tex) | [18-page PDF](ideal_four_dimensional_volume_and_structure/main.pdf) |
| Dimension bounds for hyperbolic Coxeter polytopes: 57 / 17 / 15 | [main.tex](hyperbolic_coxeter_dimension_bounds/main.tex) | [41-page PDF](hyperbolic_coxeter_dimension_bounds/main.pdf) |

Each folder contains `main.tex`, `references.bib`, `main.bbl` and the
matching `main.pdf`. The [complete 551-page dimension proof addenda](../manuscripts/addenda/dimension_bounds_20261030/main.pdf)
are retained separately. Current package paths and build receipts are in
[`release/manifest.json`](../release/manifest.json). Generated ZIPs and
checksums live in ignored `dist/`.

Open `main.tex` in a TeX editor and select pdfLaTeX with BibTeX, or run
`latexmk -pdf main.tex` from that paper's folder. The equivalent manual sequence is:

```sh
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

After a successful build, `main.bbl` supplies the bibliography for the upload.
`references.bib` remains available for editing and regeneration. No Python
program is needed to typeset these papers.

The concise dimension article presents all three bounds and their main
geometric arguments. The existing complete proof addendum has 71 pages of
main exposition followed by full case appendices. See the
[GitHub proof catalogue](hyperbolic_coxeter_dimension_bounds/ADDENDA.md) and
[dimension-paper guide](DIMENSION_BOUNDS.md) for the detailed proofs and records.

The verification suites are `right_angled_sharp` and
`ideal_ra4_volume_and_structure`. See the [reproduction guide](../repro/README.md)
for commands and requirements, and the [formal-verification guide](../formal/README.md)
for the conditional Lean scope. Detailed source versions and execution records
remain in the repository documentation. Each paper cites the repository once
in its bibliography.

The [30 October manifest](../results/arxiv_submission_current_manifest_20261030.json)
records the retained compiled editions and their successful typesetting.
The concise edition's preparation is recorded separately. Mathematical
checks and Lean execution retain their own proof maps and receipts.
