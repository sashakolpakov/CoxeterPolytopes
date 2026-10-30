# arXiv submission editions

Updated 30 October 2026. Three independent upload folders contain the current
submission sources. They present existing results from six paper families;
one combined article joins the finite-volume, compact and arithmetic compact
dimension proofs. Preparing a source folder does not submit it to arXiv.

## Sources, PDFs and results

| Paper | Complete TeX | Matching PDF | Figures | Results |
|---|---|---|---:|---|
| Sharp dimension bounds for right-angled hyperbolic polytopes | [main.tex](../arxiv/sharp_dimension_bounds_8_and_4/main.tex) | [6 pages](../manuscripts/build/arxiv_sharp_dimension_bounds_8_and_4_20261028_current_20261030.pdf) | 2 | Finite-volume right-angled dimension at most 8; entirely ideal dimension at most 4. Both bounds are attained by known examples. |
| On the second smallest ideal right-angled hyperbolic four-polytope | [main.tex](../arxiv/ideal_four_dimensional_volume_and_structure/main.tex) | [18 pages](../manuscripts/build/arxiv_ideal_four_dimensional_volume_and_structure_20261028_current_20261030.pdf) | 5 | The second ideal volume is 8π²/3, uniquely attained by the reflected 24-cell double. A separate restricted facet-spectrum theorem gives chains of 24-cells. |
| Dimension bounds for hyperbolic Coxeter polytopes | [main.tex](../arxiv/hyperbolic_coxeter_dimension_bounds/main.tex) | [551 pages](../manuscripts/build/arxiv_hyperbolic_coxeter_dimension_bounds_20261028_current_20261030.pdf) | 9 | Finite-volume dimension at most 57, compact dimension at most 17, compact arithmetic dimension at most 15. |

The [current source/PDF manifest](../results/arxiv_submission_current_manifest_20261030.json)
records the exact editions. Each upload folder contains only `main.tex`.
Figures and references are inline; no separate image, bibliography or custom
class file is needed. All three papers name **Alexander Kolpakov**, acknowledge
the use of OpenAI's GPT-6-astra and GPT-6.1-sol, and cite the repository in the
bibliography. Zip each folder's contents separately. Keep the reading PDFs and
verification receipts outside the upload folders.

## The combined dimension article

The paper has 71 pages of main exposition, followed by complete proof appendices
from page 72. Those appendices supply necessary arguments, not parked refinements.

| Part | Argument | Begins on page |
|---|---|---:|
| I | Finite-volume dimension at most 57 | 5 |
| II | Compact dimension at most 17 | 17 |
| III | Compact arithmetic dimension at most 15 | 35 |
| IV | Finite-volume proof appendices | 72 |
| V | Compact proof appendices | 107 |
| VI | Arithmetic proof appendices | 372 |

The first two results assume no arithmeticity. The arithmetic argument uses
the independently proved compact theorem within the same document. Arithmetic
reflection lattices automatically have the quadratic-form model used there.
The shared introduction explains the methodological split: local angle and
curvature budgets for finite volume, farthest-vertex component energies for
compactness, then field restrictions in the arithmetic case.

Assembly review compared all three retained main texts and appendix groups
with their reviewed illustrated sources. Labels were made distinct, companion
references became internal references, and the bibliography was merged. The
geometric arguments, finite algorithms, completeness proofs and required outputs
remain. The [dimension-paper guide](../arxiv/DIMENSION_BOUNDS.md) gives the
underlying source and proof maps. This is source-preservation and presentation
review, not a new mathematical replay.

## The right-angled articles

The dimension paper begins with the published ceilings and known examples.
An outline explains the shared maximum-height projection and the different
contradictions in the general and ideal cases. The complete cubic tensor
argument and two diagrams remain in the paper.

The volume paper presents a single sequence: the volume gap, the equality facet
spectrum, and the cutting and recognition theorem. Its five diagrams illustrate
the geometry and incidences. The chain theorem assumes every facet is an
octahedron or a doubled octahedron; it does not classify arbitrary ideal
four-polytopes as chains.

The historical introduction relates the second ideal volume to earlier extremal
questions. Potyagailo and Vinberg explicitly posed the compact 120-facet
question. The [Guo–Zheng preprint](https://arxiv.org/abs/2610.33053), Theorem 1.1
and Corollary 1.2, establishes the unique compact 120-cell facet and volume
minima. The earlier ideal theorem gives the 24-facet and volume 4π²/3 minima.
The present paper addresses the next distinct ideal volume and its equality case.

Our compact RA4 paper remains a paper of constraints, without a compact
second-volume classification. The [separate chronology audit](../results/arxiv_compact_prior_publication_audit_20261028.md)
records that this repository was created after the Guo–Zheng preprint submission;
no earlier public priority is claimed for related local results.

## Supporting checks and their scope

| Submission material | Supporting selection or map | Scope |
|---|---|---|
| Sharp dimensions | `right_angled_sharp`; [proof map](right_angled_sharp_dimensions_publication_proof_map_20261021.md) | Two exact checks: displayed 24-cell coordinates and universal rational tensor identities. |
| Ideal volume and structure | `ideal_ra4_volume_and_structure`; [proof map](ideal_ra4_volume_and_structure_proof_map.md) | Four exact checks of the saved double, equality models and specified root/opposition data. |
| Finite-volume 57 | [Reproduction map](finite_volume57_reproduce.md) | The 259-unit current-rule conjunction and earlier bounded replacements retain their explicit premises; no combined 146-job replay is claimed. |
| Compact 17 | [Complete proof map](compact17_consolidated_proof_map_20261025.md) | Full written case proofs, with selected conditional scalar/composition kernels. |
| Compact arithmetic 15 | [Dependency map](compact_arithmetic_consolidation_proof_map_20261025.md) | Four computed premises, the compact theorem and the published golden-lattice theorem remain explicit. |

The two right-angled script-only selections use small tracked inputs and require
none of the unuploaded bulk collection. The broader dimension arguments retain
the availability and execution limits stated in the [public-data guide](../repro/PUBLIC_DATA.md).
The [release evidence guide](right_angled_release_packages_20261025.md) and
[formal guide](../formal/README.md) identify the existing conditional kernels.
Geometry, literature and specified finite computations remain explicit premises.
Printed foundational axioms do not enumerate ordinary theorem parameters.

Inspect selections without launching mathematical jobs:

```sh
python3 -B scripts/run_proof_batch_20261024.py manuscript \
  right_angled_sharp --mode scripts --plan
python3 -B scripts/run_proof_batch_20261024.py manuscript \
  ideal_ra4_volume_and_structure --mode scripts --plan
```

For an actual necessary replay, consult [repro/README.md](../repro/README.md),
remove `--plan`, and supply a fresh ledger path. Execution remains serial and
resource-monitored. The catalogue's `--pdf` option typesets its preserved
5/17-page right-angled editions, not these 6/18-page sources. Historical
`compact20` and `compact_farthest_bounds` selections likewise retain their
20/23-edition scopes; they do not replay all later 15/17 appendices.

## Building the current submission PDFs

The sources are already in their separate `arxiv/` folders. Standard pdfLaTeX
packages suffice to typeset each `main.tex`; the guarded builders additionally
use `pypdf` for PDF inspection. The existing monitored launcher uses macOS
resource probes. Choose fresh tags and ledger paths: prior outputs are preserved.

```sh
python3 -B scripts/run_proof_batch_20261024.py check \
  --ledger results/sharp_submission_build_new.json -- \
  scripts/build_arxiv_submission_20261028.py \
  sharp_dimension_bounds_8_and_4 --tag reviewer1
python3 -B scripts/run_proof_batch_20261024.py check \
  --ledger results/ideal_submission_build_new.json -- \
  scripts/build_arxiv_submission_20261028.py \
  ideal_four_dimensional_volume_and_structure --tag reviewer1
python3 -B scripts/run_proof_batch_20261024.py check \
  --ledger results/dimension_submission_build_new.json -- \
  scripts/build_dimension_bounds_arxiv_20261028.py --tag reviewer1
```

Each build copies only the selected `main.tex` to an isolated temporary directory,
disables shell escape and writes the final PDF outside the upload folder.
It launches no mathematical checker or Lean module. The combined article uses
larger source/PDF/page envelopes while retaining the same serial, monitored
resource policy.

## Verification records

The current [manifest](../results/arxiv_submission_current_manifest_20261030.json)
binds the reading PDFs to the source files. Successful isolated builds have zero
warnings, no undefined references or citations, no other repository inputs and
Type 1 fonts. TeX temporary products were removed. Source review and sampled
layout inspection have separate scopes; no claim of visual inspection of every
page of the long article is made.

The current build receipts are:

- [Sharp 6-page build](../results/arxiv_sharp_current_build_20261030.json).
- [Ideal 18-page build](../results/arxiv_ideal_current_build_20261030.json).
- [Combined 551-page build](../results/arxiv_dimensions_current_build_20261030.json).

No successful mathematical worker or Lean module was rerun for preparation of
these submission editions. Existing executions retain their original sources
and declared scopes.
