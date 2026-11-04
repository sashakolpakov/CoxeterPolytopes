# Concise edition of the dimension-bound article

This edition presents the bounds **57 / 17 / 15** through their proof
strategies, central arguments and worked examples. It states the technical
propositions precisely and proves the global deductions. Complete technical
proofs, inventories and computations remain in the existing supplement;
the [addendum map](ADDENDA.md) identifies their exact proof locations.

The [current reproduction and Lean companion](../../repro/dimension_bounds/README.md)
connects those propositions to checkers, exact inputs, execution records and
selected conditional Lean kernels. It distinguishes checks available from
a public checkout from bulk-dependent regeneration.

The version cited throughout is
[`36eeaddf646eff25d4990c9ab730a20ada21f80f`](https://github.com/sashakolpakov/CoxeterPolytopes/tree/36eeaddf646eff25d4990c9ab730a20ada21f80f).
The full manuscript and the six accepted earlier editions are preserved.
This is the concise submission edition within that same paper family.

Build in this directory with standard pdfLaTeX and BibTeX:

```sh
pdflatex -no-shell-escape main.tex
bibtex main
pdflatex -no-shell-escape main.tex
pdflatex -no-shell-escape main.tex
```

The only source inputs are `main.tex` and `references.bib`; a successful
monitored build retains `main.bbl` and `main.pdf` beside the source.
Mathematical proof addenda are online
references and are not required by TeX. Repository typesetting uses the
existing serial resource monitor:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock \
  --ledger results/dimension_bounds_concise_bounded_build_attempt1_20261104.json \
  -- scripts/build_dimension_bounds_concise_20261104.py --tag attempt1
```

Run that command from the repository root. Use a fresh tag and ledger for
a revised build. The latest local upload package is
`dist/CoxeterPolytopes-dimension-bounds-arxiv.zip`, with its checksum beside
it. The builder adds `main.bbl` after successful typesetting. Build receipts
are stored under `release/dimension_bounds/`; the current artifact index is
[`release/manifest.json`](../../release/manifest.json). See the
[submission guide](../../arxiv/SUBMISSION.md).

[Read the 41-page PDF](main.pdf). The final monitored build has zero TeX or
BibTeX warnings, no undefined references or citations, and Type 1 fonts. The full technical proofs are retained in
the unchanged [complete proof addenda](../../manuscripts/addenda/dimension_bounds_20261030/main.pdf).
The release manifest records the matching PDF's page count and typesetting
receipt, including reference and font checks.

The package refers to the existing exact checker and evidence collection;
typesetting this edition does not run mathematical checks or Lean modules.
The [public-data guide](../../repro/PUBLIC_DATA.md) records which larger
inputs are not currently distributed.
