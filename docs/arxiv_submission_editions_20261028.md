# The two arXiv submission editions

28 October 2026. These are new reading and submission editions of two existing
paper families. They retain the mathematical conclusions and make the proof
strategy, external inputs and reproduction interfaces easier to follow. The
six accepted earlier editions remain preserved. Preparing these files is not
an arXiv submission, external refereeing, or a new mathematical verification.

## Sources, PDFs and results

| Paper | Complete TeX | Matching PDF | Figures | Results |
|---|---|---|---:|---|
| Sharp dimension bounds for right-angled hyperbolic polytopes | [Submission edition](../manuscripts/revisions/right_angled_sharp_dimension_bounds_arxiv_20261028.tex) | [7 pages](../manuscripts/build/arxiv_sharp_dimension_bounds_8_and_4_20261028_layout_final.pdf) | 2 | Finite-volume right-angled dimension at most 8; entirely ideal dimension at most 4. Both bounds are attained by known examples. |
| Volumes and octahedral chains of ideal right-angled hyperbolic four-polytopes | [Submission edition](../manuscripts/revisions/ideal_ra4_volume_and_structure_arxiv_20261028.tex) | [19 pages](../manuscripts/build/arxiv_ideal_four_dimensional_volume_and_structure_20261028_historical_context.pdf) | 5 | The second ideal volume is 8π²/3, uniquely attained by the reflected 24-cell double. A separate restricted facet-spectrum theorem gives chains of 24-cells. |

The [edition manifest](../results/arxiv_publication_release_manifest_20261028.json)
records exact source/PDF identities. Each source is a single self-contained
TeX file: figures are inline TikZ and the bibliography is inline. No separate
image, bibliography or custom class file is needed. Standard TeX packages
remain normal typesetting dependencies. Both papers name **Alexander Kolpakov**
and acknowledge the use of OpenAI's GPT-6-astra and GPT-6.1-sol.

The [28 October paper release](https://github.com/sashakolpakov/CoxeterPolytopes/releases/tag/arxiv-papers-2026-10-28)
collects these two editions. It is separate from the preserved research-data
releases: the 936-file bulk payload remains unuploaded.

For separate arXiv uploads, save each source as `main.tex` in its own folder
and package that folder's contents independently. The repository's local
`arxiv/` folders are upload staging, kept outside the tracked publication tree.
The public sources linked above are exact edition copies, not links into that
local staging area. Reading PDFs and verification receipts do not belong in
the TeX upload folders.

## What changed in the presentation

The dimension paper now opens with the published 12/6 ceilings and the known
8/4 examples. Its proof roadmap identifies one shared maximum-height
projection, followed by two different contradictions. The existing descent
diagram is retained, and a second diagram explains the ideal metric-equality
branch. The complete cubic tensor argument stays in the paper.

The ideal-volume paper makes the sequence of arguments explicit: first the
volume gap, then the equality facet spectrum, then the independent cutting
and recognition theorem. The compressed equality argument is expanded into
complete paragraphs. All five geometric and incidence diagrams remain.
The chain theorem still assumes every facet is an octahedron or a doubled
octahedron; it does not classify arbitrary ideal four-polytopes as chains.

The historical introduction connects the two extremal questions, least facet
number and least volume. Potyagailo and Vinberg explicitly posed the compact
120-facet question. The [Guo–Zheng preprint](https://arxiv.org/abs/2610.33053),
Theorem 1.1 and Corollary 1.2, establishes the unique compact 120-cell facet
and volume minima. Kolpakov's earlier ideal theorem establishes the analogous
24-facet and volume 4π²/3 minima. The present ideal paper continues that volume
problem with the second distinct value and its unique equality case. The
original Potyagailo–Vinberg citation is attributed to its explicit facet
question, rather than to a separate printed volume question.

Our compact RA4 manuscript remains a paper of constraints. It supplies neither
the Guo–Zheng unrestricted compact minimum proof nor a compact second-volume
classification. The [separate chronology audit](../results/arxiv_compact_prior_publication_audit_20261028.md)
records that this repository was created on 21 October, after their 19 October
preprint submission. Dated local notes do not establish earlier public priority.

## Reproducing the supporting checks

The submission sources point to the public
[companion snapshot `04c9bab7`](https://github.com/sashakolpakov/CoxeterPolytopes/tree/04c9bab7).
That pin identifies the code, finite inputs and recorded verification scopes;
it does not imply that the new submission sources existed in that older
snapshot. The [right-angled release guide](right_angled_release_packages_20261025.md)
and the following maps identify the mathematical statements checked.

| Submission paper | Catalogue ID | Script-only selection | Proof map |
|---|---|---:|---|
| Sharp dimensions | `right_angled_sharp` | 2 exact checks | [Projection, tensors and ideal equality](right_angled_sharp_dimensions_publication_proof_map_20261021.md) |
| Ideal volume and structure | `ideal_ra4_volume_and_structure` | 4 exact checks | [Defect, equality spectrum, cuts and recognition](ideal_ra4_volume_and_structure_proof_map.md) |

The dimension selection checks the displayed 24-cell coordinates and the
universal rational tensor identities. The ideal selection checks the saved
double face lattice, the small equality-model data, an independently written
saved-data check, and the specified D4 root/opposition model. These finite
checks support the displayed mathematics; they do not replace the universal
geometric arguments or the published classifications.

All inputs selected by these script-only jobs are small files already tracked
in the companion repository. They require none of the unuploaded collection
of 936 bulk research files. The workers use Python's standard library;
the guarded public entry requires Python 3.11 or later and, at this snapshot,
uses macOS-specific resource probes. The following planning commands execute
no mathematical jobs:

```sh
python3 -B scripts/run_proof_batch_20261024.py manuscript \
  right_angled_sharp --mode scripts --plan
python3 -B scripts/run_proof_batch_20261024.py manuscript \
  ideal_ra4_volume_and_structure --mode scripts --plan
```

For a necessary actual replay, remove `--plan` and supply a previously unused
ledger path. For example:

```sh
python3 -B scripts/run_proof_batch_20261024.py manuscript \
  right_angled_sharp --mode scripts \
  --ledger results/sharp_replay_new.json
python3 -B scripts/run_proof_batch_20261024.py manuscript \
  ideal_ra4_volume_and_structure --mode scripts \
  --ledger results/ideal4_replay_new.json
```

Execution is serial under the shared lock and the percentage resource policy.
The `core` selection additionally requests selected conditional Lean kernels;
`full` includes the separately marked historical model-producer job where
declared. The formal stages preserve separate compilation and declared actual
axiom audits. Reused legacy dependency objects retain the narrower provenance
recorded in the release guide. Published facts, geometric realization and the
named computational outputs remain explicit premises. Printed foundational
axioms do not enumerate ordinary theorem parameters, and no end-to-end geometric
Lean proof is claimed.

The catalogue's `--pdf` option still builds its preserved **5-page and
17-page publication editions**. It does not automatically build either
28 October submission source. The current 7/19-page reading PDFs have their
own exact source bindings and typesetting records below.

## Building the submission PDFs from a source checkout

The guarded [submission builder](../scripts/build_arxiv_submission_20261028.py)
requires `pdflatex` and the Python package `pypdf`. Its monitored execution
uses the same macOS resource interface. The successful build and PDF-inspection
environment used `pypdf==6.13.0`. Prepare the local upload folders by
copying the versioned public sources without changing their bytes:

```sh
mkdir -p arxiv/sharp_dimension_bounds_8_and_4 \
  arxiv/ideal_four_dimensional_volume_and_structure
cp manuscripts/revisions/right_angled_sharp_dimension_bounds_arxiv_20261028.tex \
  arxiv/sharp_dimension_bounds_8_and_4/main.tex
cp manuscripts/revisions/ideal_ra4_volume_and_structure_arxiv_20261028.tex \
  arxiv/ideal_four_dimensional_volume_and_structure/main.tex
```

Build the two papers separately through the existing guarded entry. Use unused
ledger filenames and tags; the builder refuses to replace a prior output:

```sh
python3 -B scripts/run_proof_batch_20261024.py check \
  --ledger results/sharp_submission_build_new.json -- \
  scripts/build_arxiv_submission_20261028.py \
  sharp_dimension_bounds_8_and_4 --tag reviewer1
python3 -B scripts/run_proof_batch_20261024.py check \
  --ledger results/ideal_submission_build_new.json -- \
  scripts/build_arxiv_submission_20261028.py \
  ideal_four_dimensional_volume_and_structure --tag reviewer1
```

This builds only the selected PDF in an isolated temporary directory. The
resulting PDF is retained under `manuscripts/build/`, outside the upload
folder. It neither runs the mathematical checks nor changes any preserved
publication source or PDF. The archived build reports retain the original
local staging-source names; the edition manifest binds their source hashes
to the byte-identical public revision files.

## Build and review evidence

Both current submission copies were compiled from only their exact `main.tex`
bytes in a fresh temporary directory, with shell escape disabled. TeX read no
other repository inputs. Both builds have zero warnings, no undefined
references or citations, and Type 1 fonts. The figures, introductions,
reproduction pages and references were visually inspected. Temporary TeX
products and reviewed preview PNGs were removed; the PDFs and receipts remain.

| Evidence | Scope |
|---|---|
| [Two-paper presentation acceptance](../results/arxiv_two_manuscripts_presentation_audit_20261028.json) | Source review, initial current source/PDF bindings, figure and reproduction mapping, visual review and cleanup. Its sharp-paper binding remains current. |
| [Sharp final typesetting execution](../results/arxiv_sharp_presentation_tex_layout_final_20261028.json) | Isolated 7-page, 2-figure build at its exact source hash; no mathematical verification claim. |
| [Sharp source/PDF build report](../results/arxiv_sharp_submission_build_20261028.json) | Original exact build JSON, with the staging-source path retained; the edition manifest binds the public source copy. |
| [Ideal historical-context acceptance](../results/arxiv_ideal_volume_historical_context_20261028.json) | Supersedes only the earlier ideal source/PDF binding; records the precise historical citations and current 19-page edition. |
| [Ideal final typesetting execution](../results/arxiv_ideal_historical_context_tex_20261028.json) | Isolated 19-page, 5-figure build at the current source hash; no mathematical verification claim. |
| [Ideal source/PDF build report](../results/arxiv_ideal_submission_build_20261028.json) | Original exact build JSON; the edition manifest joins it to the public 19-page source. |

No successful mathematical worker or Lean module was rerun for this
presentation update. Existing executions retain their original sources and
declared scopes. The preserved [5-page sharp edition](../manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf)
and [17-page ideal edition](../manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf)
remain available with their original records.
