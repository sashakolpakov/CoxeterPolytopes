# Linear programming certificates for dimension bounds of hyperbolic Coxeter polytopes

Reviewed working edition of the LP companion, 18 November 2026.

The current submission source and PDF are in
[`arxiv/linear_programming_dimension_bounds/`](../../arxiv/linear_programming_dimension_bounds/README.md),
following the same layout as the right-angled submissions. This directory
preserves the four source files and 35-page PDF consumed by the reviews below.
Edit the canonical arXiv source for subsequent submission revisions.

- [Read the manuscript — 35 pages](main.pdf)
- [Main TeX source](main.tex)
- [Mathematical and exposition review of the compression](../../audit/lp_manuscript_uniform_compression_review_20261118.json)
- [Citation and AI acknowledgment review](../../audit/lp_manuscript_citation_ack_source_review_20261118.json)
- [Latest acceptance and preservation record](../../audit/lp_manuscript_citation_ack_completion_20261118.json)
- [Final typesetting report](../../results/lp_dimension_manuscript_build07_20261118.json)

The manuscript develops one homogeneous feasibility and certificate
principle for the finite-volume, compact and compact arithmetic cases.
Every actual object supplies local data with a positive excess or deficit;
normalizing that quantity gives a feasible point. A common linear
functional, expressed as corner charges or energy capacities, contradicts
the normalization. One support-restriction proposition proves dimension
monotonicity for a fixed dimension-labelled family.

The constraints and coverage remain specific to each family. Finite volume
uses an endpoint certificate and analytic tails. Compactness uses a complete
catalogue in dimensions 18–29 and the published ceiling 29. Arithmeticity
adds its retained rank-16/17 covers and rows to that same compact catalogue.
Arithmetic pruning remains local theorem input. The global LP combines the
valid inputs; it does not newly derive them from primitive Gram positivity.

The treatment includes the corner features and network formulation; strict
energy normalization and the capacity dual; the complete compact
configuration construction; both arithmetic cover tables; monotone dimension
thresholds; and exact verification records recovering 57, 17 and 15.

The four source files are:

| File | Contents |
| --- | --- |
| `main.tex` | Introduction, geometric conventions, one homogeneous lift/certificate theorem and threshold proposition, case coverage, verification, reproduction and bibliography |
| `finite_volume.tex` | Incidence feasibility, all 26 features, signed charges, finite networks, ideal and analytic-tail blocks, vertex-budget optimization and its computed finite-volume certificate |
| `energy_framework.tex` | One fixed-budget lift/duality theorem for both compact cases, capacity ownership, an overlap example and the threshold energy objective |
| `compact_cases.tex` | Complete compact component estimates, profiles and domains, generated configuration cover, arithmetic inputs and both complete cover tables |

This revision reduces the previous 45-page draft to 35 pages at the same
11-point type size and one-inch margins. It removes repeated LP proofs,
imported geometric proof retellings and duplicated tables. Thirteen ideal
rows are specified by two exact family formulas. The arithmetic covers state
their needed floors directly; additional valid electronic-model rows may
receive weight zero. No separate data appendix was created to move pages
out of the count. Complete covers, actual-wall and rank guards, arithmetic
predicates, support ownership and worked certificates remain in the paper.

The introduction cites Guo and Zheng's *Dimension bounds for finite-volume
hyperbolic Coxeter polytopes*, arXiv:2611.01795, for their bound 81. The
acknowledgment follows the ideal right-angled manuscript's disclosure and
author-responsibility wording, adapted to the LP formulations and Python code.

The final PDF has no typesetting warnings or undefined references and uses
Type1 fonts. Five representative pages were visually inspected. TeX
intermediates and rendered review images were removed; the final PDF sits
beside its source. Typesetting ran under the existing serial 59-second,
1-GiB watchdog. No mathematical solver, Lean module or successful historical
checker was rerun to prepare this manuscript.

To rebuild from the repository root, use the existing bounded runner with
fresh report and monitor paths:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py \
  --parent-lock \
  --ledger results/NEW_lp_manuscript_monitor.json -- \
  scripts/build_lp_dimension_manuscript_20261118.py \
  --replace --report results/NEW_lp_manuscript_build.json
```

The builder copies only these four TeX inputs into an owned temporary
directory and runs up to four sequential pdfLaTeX passes. The bibliography
is included in `main.tex`; BibTeX is unnecessary. The final PDF is replaced
only after successful compilation and source-hash checks.

The [LP paradigm notes](../../notes/lp/UNIFORM_LP_PARADIGM.md) provide a
separate research guide. The [current dimension-bound submission](../../arxiv/hyperbolic_coxeter_dimension_bounds/README.md)
and its frozen verification inputs remain unchanged by this companion draft.
