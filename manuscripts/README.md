# Manuscripts

Read the **[paper guide](https://sashakolpakov.github.io/CoxeterPolytopes/papers.html)**
for the mathematical methods, assumptions, figures, page counts and evidence.
The three [arXiv submission folders](../arxiv/README.md) contain the current
single-file sources. One combined article joins the finite-volume, compact
and compact arithmetic dimension proofs; the two right-angled articles remain
separate submissions. Only these three current PDFs are kept in `build/`.

| Submission | PDF | Complete TeX |
|---|---|---|
| Dimension bounds 57 / 17 / 15 | [551 pages, 9 figures](build/arxiv_hyperbolic_coxeter_dimension_bounds_20261028_current_20261030.pdf) | [main.tex](../arxiv/hyperbolic_coxeter_dimension_bounds/main.tex) |
| Sharp right-angled dimensions 8 / 4 | [6 pages, 2 figures](build/arxiv_sharp_dimension_bounds_8_and_4_20261028_current_20261030.pdf) | [main.tex](../arxiv/sharp_dimension_bounds_8_and_4/main.tex) |
| On the second smallest ideal right-angled hyperbolic four-polytope | [18 pages, 5 figures](build/arxiv_ideal_four_dimensional_volume_and_structure_20261028_current_20261030.pdf) | [main.tex](../arxiv/ideal_four_dimensional_volume_and_structure/main.tex) |

The combined article has 71 pages of main exposition, followed by the complete
case appendices from page 72. The appendices supply necessary proof steps.
The arithmetic bound uses the general compact theorem within the same document.
Its [contents and proof maps](../arxiv/DIMENSION_BOUNDS.md) identify that dependency.

The [compact RA4 source](compact_ra4_volume_constraints.tex) and
[proof map](../docs/compact_ra4_volume_constraints_reproduce_20261019.md)
provide supporting constraints. The external
[Guo–Zheng preprint](https://arxiv.org/abs/2610.33053) establishes the compact
120-cell minimum; the compact second-volume problem remains open here.

All three submission PDFs compile from their isolated single-file sources with
zero warnings. Their figures and bibliographies are inline; all name Alexander
Kolpakov and include an AI-use acknowledgment. The
[edition and reproduction guide](../docs/arxiv_submission_editions_20261028.md)
records exact evidence and the supporting suite IDs. No new mathematical checks or
Lean modules were executed for these presentation revisions.

For serial plans and resource rules, use [repro/README.md](../repro/README.md).
Written proof, source audit, exact computation, conditional Lean and PDF
compilation retain separate verification records.

The remaining root TeX files and `revisions/` sources are verification inputs
for the [six catalogue selections](catalogue.json) and their proof maps.
They preserve the exact text consumed by recorded checks; they are not
additional current papers. Use the three `/arxiv` sources above for submission.
