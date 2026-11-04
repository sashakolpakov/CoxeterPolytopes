# Manuscripts

Read the **[paper guide](https://sashakolpakov.github.io/CoxeterPolytopes/papers.html)**
for the mathematical methods, assumptions, figures, page counts and evidence.
The three [arXiv submission folders](../arxiv/README.md) contain the current
source folders. Each has `main.tex`, `references.bib` and a generated `main.bbl`.
One combined article joins the finite-volume, compact
and compact arithmetic dimension proofs; the two right-angled articles remain
separate submissions. Each current PDF is kept beside its source in `arxiv/`.

The [concise 57/17/15 submission edition](../arxiv/hyperbolic_coxeter_dimension_bounds/README.md)
is prepared with online proof addenda. Its source archive is indexed in
[`release/manifest.json`](../release/manifest.json); the
[submission guide](../arxiv/SUBMISSION.md) selects it for the dimension-bound
paper. The [551-page complete edition](addenda/dimension_bounds_20261030/main.pdf)
remains the proof addendum, with its [source](addenda/dimension_bounds_20261030/main.tex).

| Submission | PDF | Complete TeX |
|---|---|---|
| Dimension bounds 57 / 17 / 15 | [67 pages](../arxiv/hyperbolic_coxeter_dimension_bounds/main.pdf) | [main.tex](../arxiv/hyperbolic_coxeter_dimension_bounds/main.tex) |
| Sharp right-angled dimensions 8 / 4 | [6 pages, 2 figures](../arxiv/sharp_dimension_bounds_8_and_4/main.pdf) | [main.tex](../arxiv/sharp_dimension_bounds_8_and_4/main.tex) |
| On the second smallest ideal right-angled hyperbolic four-polytope | [18 pages, 5 figures](../arxiv/ideal_four_dimensional_volume_and_structure/main.pdf) | [main.tex](../arxiv/ideal_four_dimensional_volume_and_structure/main.tex) |

The complete proof addendum has 71 pages of main exposition, followed by the complete
case appendices from page 72. The appendices supply necessary proof steps.
The arithmetic bound uses the general compact theorem within the same document.
Its [contents and proof maps](../arxiv/DIMENSION_BOUNDS.md) identify that dependency.

The [compact RA4 source](compact_ra4_volume_constraints.tex) and
[proof map](../docs/compact_ra4_volume_constraints_reproduce_20261019.md)
provide supporting constraints. The external
[Guo–Zheng preprint](https://arxiv.org/abs/2610.33053) establishes the compact
120-cell minimum; the compact second-volume problem remains open here.

The sources use standard pdfLaTeX and BibTeX. Open `main.tex` in a TeX editor,
or run `latexmk -pdf main.tex` in that paper's folder; no Python is needed.
Figures are inline TikZ; all papers name Alexander
Kolpakov and include an AI-use acknowledgment. The
[edition and reproduction guide](../docs/arxiv_submission_editions_20261028.md)
records the manual build sequence, exact evidence and supporting suite IDs.
All three four-step pdfLaTeX/BibTeX builds passed with zero warnings and no
undefined references or citations.
No new mathematical checks or
Lean modules were executed for these presentation revisions.

For serial plans and resource rules, use [repro/README.md](../repro/README.md).
Written proof, source audit, exact computation, conditional Lean and PDF
compilation retain separate verification records.

The editions selected by the [six catalogue selections](catalogue.json)
preserve the exact text consumed by recorded checks. They are verification
inputs within the same paper families. Use the
[submission guide](../arxiv/SUBMISSION.md) for the current upload selection.
