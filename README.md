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

All six editions are **written, compiled and internally reviewed**, with
matching PDFs. This describes manuscript preparation at the stated
published-input, companion-theorem and finite-data boundaries. It does
not mean complete computational replay, external refereeing or full
geometric verification in Lean; the remaining obligations are listed below.

| Paper | Result in the PDF | Read | Source |
| --- | --- | --- | --- |
| General finite volume | Dimension at most **57** | [PDF, 48 pages](manuscripts/build/finite_volume57.pdf) | [TeX](manuscripts/finite_volume57.tex) |
| Sharp right-angled dimensions | At most **8** in finite volume and **4** if entirely ideal; both sharp | [PDF, 5 pages](manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf) | [TeX](manuscripts/right_angled_sharp_dimension_bounds_publication.tex) |
| General compact | Dimension at most **23** | [PDF, 36 pages](manuscripts/build/compact_farthest_vertex_bounds.pdf) | [TeX](manuscripts/compact_farthest_vertex_bounds.tex) |
| Compact arithmetic Coxeter | Dimension at most **20** | [PDF, 21 pages](manuscripts/build/compact_arithmetic_dimension_bound.pdf) | [TeX](manuscripts/compact_arithmetic_dimension_bound.tex) |
| Ideal right-angled four-volume and structure | Published minimum, second volume and its unique double, restricted 24-cell chains | [PDF, 17 pages](manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf) | [TeX](manuscripts/ideal_ra4_volume_and_structure_publication.tex) |
| Compact right-angled four-volume constraints | Facet capacity, excess inequalities and recognition results | [PDF, 14 pages](manuscripts/build/compact_ra4_volume_constraints.pdf) | [TeX](manuscripts/compact_ra4_volume_constraints.tex) |

The finite-volume paper has 12 pages of main exposition and 36 pages
of essential proof and algorithm appendices. Each paper’s proof map is
linked from the [reader’s index](manuscripts/README.md).

Later internally reviewed written proofs lower the general compact
ceiling to **18** and the compact arithmetic Coxeter ceiling to **17**
over every actual admissible field. Those stronger proofs are supporting
research in the same two paper families; they are **not yet theorems of
the PDFs stating 23 and 20**. The [compact eighteen composition](results/compact_dimension18_global_composition_and_whole_d19_exclusion_candidate_20261022.md),
[arithmetic seventeen composition](results/compact_dimension19_and_arithmetic_dimension17_global_composition_20261021.md)
and [comparison guide](docs/current_bounds_and_verification.md) give their
precise hypotheses and review records. General finite volume remains 57.

The admissible quadratic-form model, also called arithmetic of first or
simplest type, is automatic for arithmetic hyperbolic reflection lattices.
It is not an additional restriction on the arithmetic Coxeter results.
The [external-input register](docs/manuscript_assumptions.md#arithmetic-reflection-groups-and-the-quadratic-form-model)
records this published reduction and its references.

The [latest compact acceptance](results/compact_dimension18_global_written_acceptance_20261022.json)
excludes all six remaining D19 long-head classes, including the entire
continuous dotted interval. Together with the complete dimension-nineteen
forest exclusions and the independent higher-dimensional tail, this
proves the **global compact ceiling of eighteen**, without arithmeticity.
It has full source-bound internal branch and composition reviews. The next
general frontier is dimension eighteen. Three arithmetic seventeen groups
remain: whole A17, whole B17 and whole D17. These are necessary
possibilities, not realized examples.

The [complete arithmetic D16+I exclusion](results/arithmetic17_d16_actual_i_complete_exclusion_composition_20261022.md)
now has [full outside review and root acceptance](results/d16_actual_i_exclusion_and_finite57_hexagon_root_acceptance_20261022.md).
The actual-field reduction and separate quadratic, cubic and quartic
contradictions exhaust that entire forest. Its removal leaves the three
whole families above; the arithmetic ceiling remains seventeen. The
[compact eighteen logic kernel](docs/compact_dimension_eighteen_composition_logic_20261022.md)
now has [bounded compilation and a separate five-declaration axiom audit PASS](results/lean_compact18_composition_verification_20261022/manifest.json),
with no `sorryAx`. It checks the final composition under23 named written
theorems and supplies no geometric instance. The first failed elaboration
and its repaired, independently reviewed successor are preserved.

A new [finite-volume financing theorem](results/finite57_heavy_hexagon_nonadjacent_affine_cusp_global_discharge_candidate_20261022.md)
uses two actual nonadjacent cusps to cover the heavy vertices on every
hexagon in its stated all-heavy domain. Its complete geometry, receiving
multiplicities and prior debts have full internal review and root acceptance.
Other faces and costly forests remain open; no dimension56 bound follows.

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

The distinct golden norm239 lattice has reviewed actual edges through
E8+A1, E7+I2(5), E7+A2 and E7+2A1. The
[fourth-edge proof](results/compact_h9_prime239_fourth_edge_and_e7_two_a1_germ_candidate_20261022.md)
states its complete finite algorithm and actual twenty-unit output.
The [fifth-edge proof](results/compact_h9_prime239_fifth_edge_and_e8_long_a1_germ_candidate_20261022.md)
reaches E8+A1 with an isolated square2b root; all fifteen units and the
separate combine passed, followed by full mathematical/runtime review.
The subsequent [original-module period proof](results/compact_h9_prime239_original_base_to_fifth_e8_long_a1_period_candidate_20261022.md)
and its [root acceptance](results/review_prime239_original_base_fifth_germ_period_root_acceptance_20261022.md)
now prove this literal lattice **nonreflective**: an integral future
chamber symmetry has infinite order. This conclusion concerns that
module only and is a written geometric proof, with the five earlier
edge computations as inherited inputs; no new checker or Lean run is claimed.

A different norm379 host now has a
[reviewed first transverse edge and full E7+I2(5) germ](results/review_prime379_first_transverse_wall_full_germ_root_acceptance_20261022.md).
All five profiles and the separate combination passed:55 selectors and
21,835 actual scalar rows. Its subsequent
[released-m edge now has a fully reviewed E7+2A1 germ](results/review_prime379_third_edge_and_e7_two_a1_germ_root_acceptance_20261022.md),
with one isolated norm2b wall. Five new profile runs and a separate combine
passed; the earlier RAM refusal is preserved with zero execution credit.
Reflectivity and finite volume remain open; no new polytope example is claimed.

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

For finite57, **all 25 bounded replacement units passed** and now replace
the two eager loaders in the primary catalogue. Together with inherited
records for 111 unchanged commands, they provide separate matching
evidence for the current 136 script selections. There is no combined
146-job full-manuscript runtime PASS. See the
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

Two further small kernels now have bounded compilation and separate actual
axiom audits: [D18 centre scalars](results/lean_d18_centre_verification_20261022/manifest.json)
has **five declarations**, and [cubic B16 logic](results/lean_cubic_b16_short_cap_verification_20261022/manifest.json)
has **six**. Only standard axioms occur. Their field, geometric and actual
reference-applicability premises remain explicit; the preparation headers
are preserved as historical source bytes. Each audit retains one harmless
unused-variable warning.

The [D15 typed-seven centre kernel](docs/d15_typed_seven_centre_kernel_logic_20261022.md)
also has bounded compilation and a separate **six-declaration** axiom audit
[PASS](results/lean_d15_typed_seven_centre_verification_20261022/manifest.json).
It derives the continuous radical estimate and conditional funded margins
from primitive scalar premises. Actual geometry, classification and the
reference-applicability bridge remain explicit inputs. The failed tactic
predecessor is preserved with zero verification credit.

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

A plan is an inventory, not a proof run or a successful verification
record. Finite57's [current integration receipt](results/finite57_complete_loader_catalogue_integration_20261022.json)
binds the completed bounded replacements and the actual146-job plan.
Full replay still requires the exact inherited inputs described in its guide.

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
