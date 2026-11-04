# Coxeter Polytopes

Hyperbolic Coxeter dimension bounds, right-angled volume problems, exact
certificates and selected Lean verification.

**[Read the detailed documentation](https://sashakolpakov.github.io/CoxeterPolytopes/)**
· [PDFs and editions](manuscripts/README.md)
· [Submission packages](arxiv/SUBMISSION.md)
· [Reproduction](repro/README.md)
· [Formal logic](formal/README.md)

The **[consolidated dimension-bound article](arxiv/hyperbolic_coxeter_dimension_bounds/main.tex)**
combines **57 / 17 / 15** in a [67-page PDF](arxiv/hyperbolic_coxeter_dimension_bounds/main.pdf).
It develops the geometric reductions, component estimates and arithmetic
exclusions in detail, and cites the existing
[GitHub proof catalogue](arxiv/hyperbolic_coxeter_dimension_bounds/ADDENDA.md)
for complete case proofs, computational proofs and reproduction records.
The [release index](release/manifest.json) identifies the current source
package, retained compiled editions and proof/reproduction guides.

Three **[arXiv submission folders](arxiv/README.md)** contain the current sources:
[sharp dimensions 8/4, 6 pages](arxiv/sharp_dimension_bounds_8_and_4/main.pdf),
[second ideal four-volume, 18 pages](arxiv/ideal_four_dimensional_volume_and_structure/main.pdf),
and [dimension bounds 57/17/15, 67 pages](arxiv/hyperbolic_coxeter_dimension_bounds/main.pdf).
Each PDF sits beside its source. The [551-page complete dimension proof
addenda](manuscripts/addenda/dimension_bounds_20261030/main.pdf) preserve all
case proofs separately. [Sources and verification scope](docs/arxiv_submission_editions_20261028.md)
distinguish these copies from preserved earlier editions.

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

The [repository map](docs/repository_layout.md) separates manuscript sources,
local upload archives (`dist/`), release metadata, proof records (`results/`),
and audits. Historical publication staging is compressed with verified byte
restoration; recorded scientific inputs retain their original paths.
