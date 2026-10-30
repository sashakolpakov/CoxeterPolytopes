# Results and manuscripts

Updated 30 October 2026. There are **three arXiv submission manuscripts** drawn
from six research families. Only their three current PDFs are kept in
`manuscripts/build/`. The [reader’s index](../manuscripts/README.md) lists
these submissions. Mathematical sources, proof maps and execution receipts
retain their stated scopes.

## Current submissions

| Manuscript | Results and proof organization | Source and PDF |
|---|---|---|
| Dimension bounds for hyperbolic Coxeter polytopes | Finite volume ≤57, compact ≤17, compact arithmetic ≤15. Corner payments and incidence curvature give the first bound. One farthest vector and component energies give the compact bound; field restrictions and dyadic arguments sharpen the arithmetic case. The common introduction identifies the change of method; complete case appendices follow the main exposition. | [main.tex](../arxiv/hyperbolic_coxeter_dimension_bounds/main.tex), [551-page PDF](../manuscripts/build/arxiv_hyperbolic_coxeter_dimension_bounds_20261028_current_20261030.pdf), with 71 pages before the appendices |
| Sharp dimension bounds for right-angled hyperbolic polytopes | Finite-volume dimension ≤8 and entirely ideal dimension ≤4, both attained. A shared projection leads to separate spherical-code and ideal-equality arguments. | [main.tex](../arxiv/sharp_dimension_bounds_8_and_4/main.tex), [6-page PDF](../manuscripts/build/arxiv_sharp_dimension_bounds_8_and_4_20261028_current_20261030.pdf) |
| On the second smallest ideal right-angled hyperbolic four-polytope | Second volume 8π²/3 and its unique reflected double. The volume gap and equality facet spectrum precede the restricted 24-cell chain theorem. | [main.tex](../arxiv/ideal_four_dimensional_volume_and_structure/main.tex), [18-page PDF](../manuscripts/build/arxiv_ideal_four_dimensional_volume_and_structure_20261028_current_20261030.pdf) |

Each submission contains its figures and bibliography in `main.tex`. The
[submission guide](arxiv_submission_editions_20261028.md) gives source/PDF
bindings and separate typesetting commands.

## Research families and proof maps

| Family | Current proof home | Mathematical and verification map |
|---|---|---|
| General finite-volume dimension 57 | Combined dimension article, Parts I and IV | [Proof and reproduction map](finite_volume57_reproduce.md), [lemma map](finite_volume57_lemma_map.md) |
| General compact dimension 17 | Combined dimension article, Parts II and V | [Compact-17 proof map](compact17_consolidated_proof_map_20261025.md) |
| Compact arithmetic dimension 15 | Combined dimension article, Parts III and VI; uses the general compact theorem internally | [Arithmetic-15 proof map](compact_arithmetic_consolidation_proof_map_20261025.md), [Bugaenko primary-source authentication](bugaenko_golden_lattice_cutoff_authentication_20261026.md) |
| Sharp right-angled dimensions 8 / 4 | Sharp-dimension submission | [Projection, tensor and ideal-equality map](right_angled_sharp_dimensions_publication_proof_map_20261021.md) |
| Ideal RA4 volume and structure | Second-volume submission | [Defect, equality spectrum, cuts and recognition](ideal_ra4_volume_and_structure_proof_map.md) |
| Compact RA4 volume constraints | [Supporting source](../manuscripts/compact_ra4_volume_constraints.tex); not a fourth submission | [Constraint proof and verification map](compact_ra4_publication_proof_and_verification_map_20261024.md) |

The source labels in the dimension maps remain useful: the combined submission
prefixes them with `fv:`, `compact:` and `arith:`. Dated receipts continue to
identify their exact source editions.

## Verification boundaries

The results are internally reviewed written arguments at their stated
geometric, literature and finite-data premises. Source review, exact checks,
conditional Lean verification and PDF compilation have separate records.

- **Finite volume 57:** the separate current-rule conjunction has 259 actual
  scoped PASS units, all 64 selected shards, 28,028 maps, 12,518 nonunit
  applications and 83 prefixes. It does not regenerate the native classification
  or inherited numerical data, and is not a combined fresh 146-job replay.
- **Compact 17 and arithmetic 15:** the complete case proofs and global
  compositions have source reviews and root acceptance. Arithmetic 15 retains
  four computed propositions and the authenticated external lattice theorem.
  The final whole-D16 argument is written, rather than an additional census.
- **Right-angled papers:** exact coordinate, tensor and finite-model checks
  support the written proofs. Published classifications and universal geometric
  arguments retain their explicit roles.
- **Lean:** selected modules check internal implications with geometric,
  literature and specified enumeration statements supplied as premises.
  Separate actual axiom audits do not enumerate ordinary theorem parameters.
  There is no end-to-end geometric formalization of these manuscripts.
- **Current PDFs:** the builds and source/PDF bindings establish typesetting
  and document identity. They add no new mathematical replay.

The [catalogue](../manuscripts/catalogue.json) retains the named verification
selections. Its `--pdf` option checks the selected source in a temporary
directory; it does not build a current submission. The
[reproduction guide](../repro/README.md) explains the separate entry points.
Intermediate case exclusions and failed attempts remain in their dated
scientific records; they are not a second list of current open cases.
