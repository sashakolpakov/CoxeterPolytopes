# Coxeter Polytopes

Dimension bounds, geometry and exact proof tools for hyperbolic Coxeter
polytopes. This repository contains six mathematical manuscripts, their
PDFs, supporting proofs, reproducibility records and selected Lean kernels.
The main results concern dimension bounds and right-angled polytopes;
the four-dimensional papers also study volume and combinatorial structure.

Start with the [documentation website](https://sashakolpakov.github.io/CoxeterPolytopes/)
([Sphinx source](docs/site/index.rst)), the
[six-paper reader’s index](manuscripts/README.md), or the
[current bounds and verification guide](docs/current_bounds_and_verification.md).
The repository is [sashakolpakov/CoxeterPolytopes](https://github.com/sashakolpakov/CoxeterPolytopes).

## Papers

All six editions are complete at their stated published-input,
companion-theorem and finite-data boundaries, with matching PDFs.

| Paper | Result in the PDF | Read | Source |
| --- | --- | --- | --- |
| General finite volume | Dimension at most **57** | [PDF, 48 pages](manuscripts/build/finite_volume57.pdf) | [TeX](manuscripts/finite_volume57.tex) |
| Sharp right-angled dimensions | At most **8** in finite volume and **4** if entirely ideal; both sharp | [PDF, 5 pages](manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf) | [TeX](manuscripts/right_angled_sharp_dimension_bounds_publication.tex) |
| General compact | Dimension at most **23** | [PDF, 36 pages](manuscripts/build/compact_farthest_vertex_bounds.pdf) | [TeX](manuscripts/compact_farthest_vertex_bounds.tex) |
| Compact arithmetic, first type | Dimension at most **20** | [PDF, 21 pages](manuscripts/build/compact_arithmetic_dimension_bound.pdf) | [TeX](manuscripts/compact_arithmetic_dimension_bound.tex) |
| Ideal right-angled four-volume and structure | Published minimum, second volume and its unique double, restricted 24-cell chains | [PDF, 17 pages](manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf) | [TeX](manuscripts/ideal_ra4_volume_and_structure_publication.tex) |
| Compact right-angled four-volume constraints | Facet capacity, excess inequalities and recognition results | [PDF, 14 pages](manuscripts/build/compact_ra4_volume_constraints.pdf) | [TeX](manuscripts/compact_ra4_volume_constraints.tex) |

The finite-volume paper has 12 pages of main exposition and 36 pages
of essential proof and algorithm appendices. Each paper’s proof map is
linked from the [reader’s index](manuscripts/README.md).

Later internally reviewed written proofs lower the general compact
ceiling to **19** and the compact arithmetic first-type ceiling to **17**
over every actual admissible field. Those stronger proofs are supporting
research in the same two paper families; they are **not yet theorems of
the PDFs stating 23 and 20**. The [global composition](results/compact_dimension19_and_arithmetic_dimension17_global_composition_20261021.md)
and [comparison guide](docs/current_bounds_and_verification.md) give their
precise hypotheses and review records. General finite volume remains 57.

The [latest accepted supporting batch](results/b18_a11_and_restricted_continuations_acceptance_20261021.json)
also excludes general B18+actual-I and arithmetic A11+actual-I+(A5 or D5).
Five general nineteen and seven arithmetic seventeen groups remain
necessary possibilities; these are not realized examples. The batch's
[result and evidence map](results/publication_b18_a11_batch_evidence_map_20261021.md)
separates these exclusions from three restricted geometric continuations.

The right-angled 8 and entirely ideal 4 bounds match published examples.
No new polytope example is claimed. The unrestricted compact right-angled
120-cell minimum and second-volume problems remain open. Preserved
field-specific arithmetic refinements are excluded from the current
universal paper; they are not additional papers.

The literal golden norm79 construction host is now proved nonreflective:
two full actual D9 germs give an original-lattice integral chamber period
of infinite order. The [proof](results/compact_h9_prime79_two_actual_d9_germ_period_candidate_20261021.md)
has full internal reviews and separate exact execution records. Its actual-germ
and arithmetic-covolume inputs are explicit; no field, form-class, other
lattice or infinite-family exclusion follows.

## What “verified” means here

The repository keeps four kinds of evidence separate:

- **Written proofs:** the manuscripts and supporting notes state their
  geometric, published and finite-data inputs.
- **Internal mathematical reviews:** source-bound reviews assess complete
  arguments at those boundaries. They are not external refereeing.
- **Exact computation:** scripts check specified identities, finite domains
  and certificates. Recorded replay does not imply regeneration of every
  transitive input.
- **Partial Lean formalization:** selected compiled kernels prove internal
  implications under explicit premises. Actual hyperbolic geometry,
  literature inputs and concrete application bridges remain partly outside
  Lean; prepared modules earn no compilation credit.

For finite57, matching bounded records cover **111 of 113 selected script
commands**. Two legacy eager loaders still need replacement in the full
selection; their reviewed 25-unit replacement has **2 of 25 actual PASS
records**. There is no combined full-manuscript runtime PASS. See the
[finite57 reproduction guide](docs/finite_volume57_reproduce.md).

The [verification guide](docs/manuscript_verification.md),
[assumption register](docs/manuscript_assumptions.md) and per-paper proof
maps explain what each recorded check establishes. The
[Lean directory](formal/README.md) contains the selected logic sources and
coverage information. Written proof completion does not mean that Lean
formalization is complete.

The Regular24 and classical-composition kernels have bounded compilation and separate
actual axiom-audit PASS records: Regular24 has **five actual declarations**,
and ClassicalDimensionComposition has **four**. The
[completion receipt](results/lean_publication_completion_20261021/manifest.json)
binds unchanged sources, individual compile/audit ledgers and actual output.
Only standard axioms occur; no successful dependencies were rerun. These
are conditional coordinate and composition results, with geometry, actual
field interpretation and script premises still explicit. They discharge
the two selected compilation/audit prerequisites, not complete geometric
formalization or remote publication.

The additional [endpoint triangle kernel](docs/endpoint_triangle_norm_kernel_logic_20261021.md)
now has bounded compilation and a separate **ten-declaration** axiom audit
PASS. It derives cycle cancellation, an integral Bezout witness and finite
norm contradictions. Its [receipt](results/lean_endpoint_triangle_verification_20261021/manifest.json)
records the final source, actual output and preserved failed attempts;
geometric and number-field application inputs remain explicit.

The [Lorentz obstruction kernel](docs/lorentz_obstruction_kernel_logic_20261021.md)
also has bounded compilation and a separate **seven-declaration** axiom audit
[PASS](results/lean_lorentz_obstruction_verification_20261021/manifest.json),
using only standard axioms and no `sorryAx`. It verifies the conditional
second-orbit cutoff and restricted Schur/rank algebra. Actual finite-order
averaging, chamber/lattice geometry and arithmetic interpretation remain
written inputs; it is not a complete geometric formalization.

## Large research inputs

Large scientific inputs are packaged separately from Git. The
[data guide](repro/PUBLIC_DATA.md) describes the versioned checksum manifest,
sequential downloads and restoration to the original paths. Inspect the
required assets without network requests or file writes:

```sh
python3 -B scripts/fetch_public_data.py --plan
```

Release assets have not yet been uploaded. Their preparation is not a
successful public restore or a mathematical replay. The guide explains
selective restoration and checksum validation once the assets are available.

## Inspect a verification plan

From the repository root, Python 3.11 or later can list the available
paper IDs and print plans without starting proof processes:

```sh
python3 -B scripts/verify_manuscript.py --list
python3 -B scripts/verify_manuscript.py right_angled_sharp --mode core --pdf --plan
python3 -B scripts/verify_manuscript.py finite57 --mode full --pdf --plan
```

**Keep `--plan` on the finite57 full command.** Its current selection
still contains the two eager loaders described above. A plan is an
inventory, not a proof run or a successful verification record.

For actual replay, follow the relevant proof map and
[bounded reproduction instructions](docs/manuscript_verification.md).
Checks run serially through the repository’s resource-controlled runners,
with 40 GiB free disk and 6 GiB available RAM required for admission.
Do not substitute uncapped or parallel runs. Lean uses the
[pinned toolchain](formal/lean-toolchain); dependencies and additional
checker requirements are documented with their respective packages.

## Repository guide

- [manuscripts/](manuscripts/README.md): the six papers, PDFs and preserved supporting editions.
- [docs/](docs/site/index.rst): reading guides, methods, assumptions and reproduction contracts.
- [formal/](formal/README.md): conditional Lean kernels and their assumption maps.
- [scripts/](scripts/): exact checkers and bounded orchestration tools.
- [results/](results/): supporting written proofs, certificates, reviews and execution receipts.
