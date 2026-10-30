# arXiv submission sources

Prepare one submission from each of these folders:

| Paper | Source | Current reading PDF |
|---|---|---|
| Sharp dimension bounds 8 and 4 | [main.tex](sharp_dimension_bounds_8_and_4/main.tex) | [6-page PDF](../manuscripts/build/arxiv_sharp_dimension_bounds_8_and_4_20261028_current_20261030.pdf) |
| On the second smallest ideal right-angled hyperbolic four-polytope | [main.tex](ideal_four_dimensional_volume_and_structure/main.tex) | [18-page PDF](../manuscripts/build/arxiv_ideal_four_dimensional_volume_and_structure_20261028_current_20261030.pdf) |
| Dimension bounds for hyperbolic Coxeter polytopes: 57 / 17 / 15 | [main.tex](hyperbolic_coxeter_dimension_bounds/main.tex) | [551-page PDF](../manuscripts/build/arxiv_hyperbolic_coxeter_dimension_bounds_20261028_current_20261030.pdf) |

Each folder contains only `main.tex`. Figures and references are
inline, and all papers include the author, AI-use acknowledgment and repository
links. Zip each folder's contents separately. Keep the reading PDFs and build
records outside the upload folders.

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

All three submission copies were compiled from their isolated `main.tex` files with
zero warnings. These are presentation revisions of the existing papers; this
batch does not claim new mathematical-check or Lean executions.
