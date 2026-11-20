# Submission sources and packages

**The consolidated dimension-bound package is prepared.** The proofread
41-page paper contains all three bounds, their main geometric arguments and
precise references to the GitHub proof catalogue. The article explains the proof strategies and representative geometric and
arithmetic arguments, states the technical propositions, and proves the
global deductions. Complete technical proofs,
computational proofs and reproduction records remain in the repository.
The existing 551-page complete edition is preserved byte for byte.

There are four prepared manuscripts within the existing three paper families.
Local preparation has not submitted anything
to arXiv. Upload packages belong in the ignored `dist/` directory; the
[release manifest](../release/manifest.json) is the current artifact index.

| Paper | Submission source | Detailed proofs and reading copy |
|---|---|---|
| Dimension bounds **57 / 17 / 15** | [Concise main.tex](hyperbolic_coxeter_dimension_bounds/main.tex) and [bibliography](hyperbolic_coxeter_dimension_bounds/references.bib) | [Addendum map](hyperbolic_coxeter_dimension_bounds/ADDENDA.md), [complete source](../manuscripts/addenda/dimension_bounds_20261030/main.tex), [551-page complete proof PDF](../manuscripts/addenda/dimension_bounds_20261030/main.pdf) |
| Sharp right-angled dimensions **8 / 4** | [main.tex](sharp_dimension_bounds_8_and_4/main.tex), with the bibliography in that folder | [Retained 6-page PDF](sharp_dimension_bounds_8_and_4/main.pdf) |
| Second ideal right-angled four-volume | [main.tex](ideal_four_dimensional_volume_and_structure/main.tex), with the bibliography in that folder | [Retained 18-page PDF](ideal_four_dimensional_volume_and_structure/main.pdf) |
| LP certificates for dimension bounds **57 / 17 / 15** | [main.tex](linear_programming_dimension_bounds/main.tex) and [bibliography](linear_programming_dimension_bounds/references.bib) | [36-page PDF](linear_programming_dimension_bounds/main.pdf), [source and package guide](linear_programming_dimension_bounds/README.md) |

## Combined dimension-bound submission

The concise paper includes the general finite-volume improvement of
Prokhorov, the general compact bound, and the compact arithmetic bound.
Complete case proofs and the reproduction package are referenced at a fixed
repository version. The 551-page edition serves as the full proof addendum.

The latest local source archive is:

```text
dist/CoxeterPolytopes-dimension-bounds-arxiv.zip
dist/CoxeterPolytopes-dimension-bounds-arxiv.zip.sha256
```

The archive contains exactly `main.tex`, `references.bib` and `main.bbl`.
The 41-page PDF and prepared bibliography are beside the source. The final
four-pass build passed with zero TeX/BibTeX warnings and no undefined
references or citations. The full 551-page proof addendum is retained
separately, with its original bytes.

**Title:** Dimension bounds for hyperbolic Coxeter polytopes
**Author:** Alexander Kolpakov

The [concise source guide](hyperbolic_coxeter_dimension_bounds/README.md)
contains the standard TeX commands and the monitored repository command.

## Linear programming companion

**Title:** Linear programming certificates for dimension bounds of hyperbolic Coxeter polytopes
**Author:** Alexander Kolpakov

The 36-page companion recovers the existing bounds 57/17/15 from the proved
local inequalities and complete domains. One feasibility and certificate
framework performs the global combination. The arithmetic pruning and other
local theorem inputs retain their stated roles; this preparation runs no
mathematical checkers.

The canonical source folder is
[`arxiv/linear_programming_dimension_bounds/`](linear_programming_dimension_bounds/README.md).
It contains a single `main.tex`, `references.bib`, the generated `main.bbl`
and the matching `main.pdf`. The earlier reviewed four-file working edition
and its PDF remain in
[`manuscripts/linear_programming_dimension_bounds/`](../manuscripts/linear_programming_dimension_bounds/README.md).

The local upload archive and checksum are:

```text
dist/CoxeterPolytopes-linear-programming-dimension-bounds-arxiv.zip
dist/CoxeterPolytopes-linear-programming-dimension-bounds-arxiv.zip.sha256
```

The ZIP contains exactly `main.tex`, `references.bib` and `main.bbl` at
its root. Select `main.tex` and PDFLaTeX on upload. PDFs, editor output,
research data and Lean products are excluded. Build and package receipts
are under `release/linear_programming_dimension_bounds/`, with their current
paths and hashes in the [release manifest](../release/manifest.json).

**Comments:** 36 pages; exact LP certificates, local theorem sources and
reproduction records at https://github.com/sashakolpakov/CoxeterPolytopes.

## Bibliography and compiled editions

Each submission folder contains `main.tex`, `references.bib`, `main.bbl`
and its matching PDF. The
[release manifest](../release/manifest.json) binds the current source/PDF
pairs and build receipts. The [30 October manifest](../results/arxiv_submission_current_manifest_20261030.json)
retains the historical compiled editions; it is not the upload selection.

Standard typesetting uses PDFLaTeX, BibTeX, and two further PDFLaTeX passes.
Figures are inline. Source archives exclude PDFs, TeX intermediates, research
data and Lean products. Mathematical and Lean verification retain their
separate recorded scopes in the [reproduction guide](../repro/README.md).

## Upload metadata for the dimension-bound article

**Comments:** 41 pages; complete proof addenda, computational proofs and
reproduction package at https://github.com/sashakolpakov/CoxeterPolytopes.

Select `main.tex` as the primary source and PDFLaTeX as the compiler.
The upload ZIP contains exactly the three source/bibliography files listed
above. The [proof catalogue](hyperbolic_coxeter_dimension_bounds/ADDENDA.md#what-is-in-the-article-and-what-remains-online)
identifies the central arguments and the technical proofs retained online.
The [current reproduction and Lean companion](../repro/dimension_bounds/README.md)
provides the proposition-to-checker map and public-checkout instructions.
