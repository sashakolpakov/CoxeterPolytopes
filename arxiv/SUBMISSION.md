# Submission sources and packages

**The consolidated dimension-bound package is prepared.** The proofread
15-page paper contains all three bounds, their main geometric arguments and
precise references to the GitHub proof catalogue. Complete case proofs,
computational proofs and reproduction records remain in the repository.
The existing 551-page complete edition is preserved byte for byte.

There are three paper families. Local preparation has not submitted anything
to arXiv. Upload packages belong in the ignored `dist/` directory; the
[release manifest](../release/manifest.json) is the current artifact index.

| Paper | Submission source | Detailed proofs and reading copy |
|---|---|---|
| Dimension bounds **57 / 17 / 15** | [Concise main.tex](hyperbolic_coxeter_dimension_bounds/main.tex) and [bibliography](hyperbolic_coxeter_dimension_bounds/references.bib) | [Addendum map](hyperbolic_coxeter_dimension_bounds/ADDENDA.md), [complete source](../manuscripts/addenda/dimension_bounds_20261030/main.tex), [551-page complete proof PDF](../manuscripts/addenda/dimension_bounds_20261030/main.pdf) |
| Sharp right-angled dimensions **8 / 4** | [main.tex](sharp_dimension_bounds_8_and_4/main.tex), with the bibliography in that folder | [Retained 6-page PDF](sharp_dimension_bounds_8_and_4/main.pdf) |
| Second ideal right-angled four-volume | [main.tex](ideal_four_dimensional_volume_and_structure/main.tex), with the bibliography in that folder | [Retained 18-page PDF](ideal_four_dimensional_volume_and_structure/main.pdf) |

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
The 15-page PDF and prepared bibliography are beside the source. The final
four-pass build passed with zero TeX/BibTeX warnings and no undefined
references or citations. The full 551-page proof addendum is retained
separately, with its original bytes.

**Title:** Dimension bounds for hyperbolic Coxeter polytopes
**Author:** Alexander Kolpakov

The [concise source guide](hyperbolic_coxeter_dimension_bounds/README.md)
contains the standard TeX commands and the monitored repository command.

## Bibliography and compiled editions

The two right-angled paper folders contain `main.tex`, `references.bib`,
`main.bbl` and their matching PDFs. The
[release manifest](../release/manifest.json) binds the current source/PDF
pairs and build receipts. The [30 October manifest](../results/arxiv_submission_current_manifest_20261030.json)
retains the historical compiled editions; it is not the upload selection.

Standard typesetting uses PDFLaTeX, BibTeX, and two further PDFLaTeX passes.
Figures are inline. Source archives exclude PDFs, TeX intermediates, research
data and Lean products. Mathematical and Lean verification retain their
separate recorded scopes in the [reproduction guide](../repro/README.md).
