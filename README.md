# Coxeter Polytopes

Hyperbolic Coxeter dimension bounds, right-angled volume problems, exact
certificates and selected Lean verification. There are **six paper families**.

**[Read the detailed documentation](https://sashakolpakov.github.io/CoxeterPolytopes/)**
· [PDFs and editions](manuscripts/README.md)
· [Reproduction](repro/README.md)
· [Formal logic](formal/README.md)

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
no new polytope example. Arithmetic 15 retains an explicit external
golden-lattice cutoff assumption; its primary text is not yet authenticated.
Arithmetic reflection lattices automatically admit the quadratic-form model.

The ideal RA4 paper treats the second volume and its unique double, and a
restricted 24-cell chain theorem. The compact RA4 paper proves constraints;
it does not settle the minimum or second volume of compact RA4 polytopes.

Inspect verification plans without launching proof jobs:

```sh
python3 -B scripts/verify_manuscript.py --list
python3 -B scripts/verify_manuscript.py finite57 --mode core --plan
```

Actual execution is serial and monitored. Bulk research assets are not
uploaded; [availability and exact replay limits](repro/PUBLIC_DATA.md) are
explicit. The documentation separates source review, computation, Lean,
PDF compilation and publication evidence.
