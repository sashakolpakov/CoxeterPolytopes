# Coxeter Polytopes

Hyperbolic Coxeter dimension bounds, right-angled volume problems, exact
certificates and selected Lean verification. There are **six paper families**.

**[Read the detailed documentation](https://sashakolpakov.github.io/CoxeterPolytopes/)**
· [PDFs and editions](manuscripts/README.md)
· [Reproduction](repro/README.md)
· [Formal logic](formal/README.md)

The **28 October arXiv submission editions** are
[sharp dimensions 8/4, 7 pages](manuscripts/build/arxiv_sharp_dimension_bounds_8_and_4_20261028_layout_final.pdf)
and [ideal volume and structure, 19 pages](manuscripts/build/arxiv_ideal_four_dimensional_volume_and_structure_20261028_historical_context.pdf).
[Sources, figures and verification scope](docs/arxiv_submission_editions_20261028.md)
and the [paper release](https://github.com/sashakolpakov/CoxeterPolytopes/releases/tag/arxiv-papers-2026-10-28)
distinguish these copies from the preserved earlier editions.

| Class | Largest published example | Published ceiling | Repository ceiling |
|---|---:|---:|---:|
| General finite volume | 21 | 995 | **57** |
| General compact | 8 | 29 | **17** |
| Compact arithmetic | 8 | 29 | **15** |
| Noncompact arithmetic | 21 | 21 | **21**, published and sharp |
| Right-angled finite volume | 8 | 12 | **8**, internally sharp |
| Entirely ideal right-angled | 4 | 6 | **4**, internally sharp |
| Compact right-angled | 4 | 4 | **4**, published and sharp |

Repository improvements are internally reviewed written arguments at their
[stated inputs](docs/manuscript_assumptions.md). Selected Lean kernels verify
conditional internal logic; the entire geometry is not formalized. There is
no new polytope example. The published golden-lattice theorem used by
arithmetic 15 now has a [primary-source check](docs/bugaenko_golden_lattice_cutoff_authentication_20261026.md).
Arithmetic reflection lattices automatically admit the quadratic-form model.

The ideal RA4 paper treats the second volume and its unique double, and a
restricted 24-cell chain theorem. The compact RA4 paper proves constraints;
the compact 120-cell facet and volume minima are established in the
[Guo–Zheng preprint](https://arxiv.org/abs/2610.33053). The compact second-volume
problem remains outside this repository's proved results.

Inspect verification plans without launching proof jobs:

```sh
python3 -B scripts/verify_manuscript.py --list
python3 -B scripts/verify_manuscript.py finite57 --mode core --plan
```

Actual execution is serial and monitored. Bulk research assets are not
uploaded; [availability and exact replay limits](repro/PUBLIC_DATA.md) are
explicit. The documentation separates source review, computation, Lean,
PDF compilation and publication evidence.
