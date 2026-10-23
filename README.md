# Coxeter Polytopes

Dimension bounds, explicit hyperbolic reflection chambers and right-angled
volume problems. The repository has **six paper families**, their frozen
TeX/PDF editions, supporting proofs, reproducibility records and selected
Lean verification.

Start with the [documentation website](https://sashakolpakov.github.io/CoxeterPolytopes/),
the [current bounds and verification guide](docs/current_bounds_and_verification.md)
or the [paper index](manuscripts/README.md).

## Current dimension results

The repository improvements are **internally reviewed written proofs**:
complete arguments at their declared input boundaries, another author's
source review and coordinating acceptance. Publication, external refereeing,
script replay and geometric formalization have separate statuses.

| Class | Largest published example | Published ceiling | Repository ceiling |
|---|---:|---:|---:|
| General Coxeter, finite volume | 21 | 995 | **57** |
| General compact Coxeter | 8 | 29 | **18** |
| Compact arithmetic Coxeter | 8 | 29 | **16** |
| Noncompact arithmetic Coxeter | 21 | 21 | **21**, published and sharp |
| Right-angled, finite volume, all vertex types | 8 | 12 | **8**, internally sharp |
| Right-angled, entirely ideal | 4 | 6 | **4**, internally sharp |
| Right-angled, compact | 4 | 4 | **4**, published and sharp |

Arithmetic reflection lattices automatically admit the quadratic-form
model; “first type” adds no restriction to these arithmetic rows. Exact
published inputs and scopes are in the [assumption register](docs/manuscript_assumptions.md).
There is **no new polytope example** from the active construction search.

## The six papers

All six editions are **written, compiled and internally reviewed**. Later
compact18/arithmetic16 proofs belong in the existing compact families;
the frozen PDFs still state23/20. Earlier notes, fragments and atlases
are supporting material, rather than additional proposed papers.

| Paper | Result in the PDF | PDF |
|---|---|---|
| General finite volume | Dimension at most57 | [48 pages](manuscripts/build/finite_volume57.pdf) |
| Sharp right-angled dimensions | Finite volume at most8; entirely ideal at most4 | [5 pages](manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf) |
| General compact | Dimension at most23; supporting proof18 | [36 pages](manuscripts/build/compact_farthest_vertex_bounds.pdf) |
| Compact arithmetic | Dimension at most20; supporting proof16 | [21 pages](manuscripts/build/compact_arithmetic_dimension_bound.pdf) |
| Ideal right-angled four-volume and structure | Published minimum; second volume and unique double; restricted24-cell chains | [17 pages](manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf) |
| Compact right-angled four-volume constraints | Facet capacity, excess and recognition | [14 pages](manuscripts/build/compact_ra4_volume_constraints.pdf) |

The57 paper has12 pages of main exposition and36 pages of essential proof
and algorithm appendices. The [paper index](manuscripts/README.md) and
[result correspondence](docs/result_manuscript_correspondence.md) identify
each technique, source, proof map and verification boundary. The ideal
four-volume paper does not classify all ideal polytopes as24-cell chains.
The compact four-volume paper does not establish a120-cell minimum or
a unique second-volume double.

## What the proofs do

The finite-volume argument assigns weights to actual polygon corners.
Spherical and affine vertex diagrams bound the available payment; global
face counting requires more. It is one simultaneous rule, including every
shared debt. The compact arguments use one farthest vertex and one spatial
vector: orthogonal component projections have additive squared norms and
force a cost above its budget. Arithmeticity adds field and integrality
restrictions. The sharp right-angled proof uses actual facet projections
followed by spherical-code bounds.

At exact compact dimension18, only **wholeB18 or wholeD18** remains in
the necessary farthest-vertex list. These are hypothetical configurations,
not examples. The [complete B18 terminal inventory](results/review_compact18_b18_uniform_terminal_lengths_complete_seventeen_inventory_root_acceptance_20261023.md)
gives unconditional same-vector sums S17<3.760 and S18<3.830.
The [near3/own6 role](results/review_compact18_b18_near_three_own_six_complete_short_exclusion_root_acceptance_20261023.md)
and [all finite-own far-helper roles](results/review_compact18_b18_all_finite_far_helper17_raw_prefix_root_acceptance_20261023.md)
are excluded. The [complete dotted far-helper exclusion](results/review_compact18_b18_dotted_far_helper17_complete_root_acceptance_20261023.md)
now leaves ONLY the short-wall contact family for the seventeenth opposite.
The [D18 inventory](results/review_compact18_d18_complete_centre16_inventory_two_holdouts_root_acceptance_20261023.md)
and [subsequent reduction](results/review_compact18_d18_nonpure_strong15_finite_fork3_head_root_acceptance_20261023.md)
leave two stated weak-centre profiles. Neither whole family is excluded.

Arithmetic16 retains [nine original groups](results/review_arithmetic16_a13_actual_i_ordinary_dihedral_group_nine_root_acceptance_20261023.md).
Universal A12/A13/A14 energies and actual-field constraints remove their
stated alternatives. The [fresh B16 head theorem](results/review_arithmetic16_whole_b16_head_near_short_three_scale_root_acceptance_20261023.md)
and [complete short-family reduction](results/review_arithmetic16_whole_b16_complete_short_white_isolated_dotted_field_root_acceptance_20261023.md)
now leave both actual opposites isolated, the short one dotted. The actual
field is quadratic or cubic; its two own parameters satisfy strict
conjugate bounds and an additional half-product integrality condition.
Every nonisolated short role is excluded. The isolated pair, wholeB16,
original group15 and the ceiling16 remain open.

Toward finite volume56, the allocation arguments now finance
[all zero-exchanged row7 and near-row8 exits](results/review_finite57_zero_exchanged_all_row7_near_row8_financing_root_acceptance_20261023.md),
retaining every earlier payment. Far-row8 and positive exchanged products
remain among the explicit obstacles. No universal364 cap or56 theorem is
claimed.

The construction lane walks actual edges in original lattices. Both the
[glued host](results/review_glued_e7_a2_two_unit_releases_original_period_and_literal_nonreflectivity_root_acceptance_20261023.md)
and the [squared-unit host](results/review_h9_squared_unit_e7_return_original_period_literal_nonreflectivity_root_acceptance_20261023.md)
are proved nonreflective by original chamber symmetries of infinite order.
These exclude two literal hosts, not their commensurability classes.
A different [prime41 host and first edge](results/review_h9_cubed_unit_prime41_original_host_first_edge_root_acceptance_20261023.md)
and its [second actual unit edge](results/review_h9_prime41_second_actual_unit_edge_original_germ_obstruction_root_acceptance_20261023.md)
are proved. The [third edge and original involution](results/review_h9_prime41_third_actual_edge_original_involution_full_hexagon_root_acceptance_20261023.md)
close its WHOLE compact E7 hexagon, of area4pi/3. The involution has
order two; ambient volume remains unresolved.

## Verification and reproduction

The [computational proof contract](docs/computational_proof_contract.md)
requires each algorithm's domain, completeness and termination argument,
exact arithmetic, execution result and replay boundary. Hash agreement and
PDF compilation establish neither mathematical completeness nor correctness.

Finite57's25 bounded B13/frontier replacement units have all actually
passed and are integrated in the primary catalogue. The146-job full
selection has no combined runtime PASS or complete transitive regeneration;
111 other command matches retain their inherited execution scope. See
the [reproduction guide](docs/finite_volume57_reproduce.md).

Selected Lean kernels have actual compilation and **separate actual axiom
audits**. Their logic maps expose the ordinary theorem parameters that
printed foundational axioms do not enumerate.

| Kernel | Actual audited declarations | Scope and receipt |
|---|---:|---|
| B18 terminal pressure | 10 | [Continuous scalar pressure contradictions](results/lean_b18_terminal_pressure_verification_20261023/manifest.json) |
| A13 endpoint graft | 11 | [Continuous gains and paid-square logic](results/lean_a13_i_endpoint_graft_verification_20261023/manifest.json) |
| A12 adjusted loss | 10 | [Continuous loss and finite-list paid square](results/lean_a12_adjusted_loss_verification_20261023/manifest.json) |
| A18 dotted endpoints | 10 | [Continuous endpoint inequalities](results/lean_a18_two_dotted_endpoint_verification_20261023/manifest.json) |
| Heavy-centre forest | 11 | [Cone, conditional chain and incidence algebra](results/lean_heavy_centre_forest_verification_20261023/manifest.json) |

The new B18 module compiled in1.71seconds; its separate audit took1.16seconds
and reused the object. All ten prints use only standard foundational axioms,
with no `sorryAx`. Two unused-variable warnings are preserved. Actual geometry,
inventories, legal pressure, terminal sums and real instances remain external;
this is not wholeB18 or end-to-end manuscript formalization. Earlier failures
and resource deferrals retain zero PASS credit. The [formal guide](formal/README.md)
records the complete coverage.

For named, serial reproduction batches:

```sh
python3 -B scripts/verify_manuscript.py finite57 --mode core
python3 -B scripts/verify_manuscript.py right_angled_sharp --mode core
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.B18TerminalPressureKernel
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.B18TerminalPressureKernel --audit-axioms
```

Successful modules and dependencies need not be rebuilt without changed inputs.

## Resources and publication

Run one scientific, Lean or TeX child at a time under the shared lock.
Admission requires40GiB free disk and6GiB available RAM. Script/TeX batches
use59seconds/1GiB/one thread; Lean uses2048MiB/one thread and a separate
axiom audit. Remove owned unused products after each batch. Start cleanup
below20GiB free disk, preserving sources, certificates, final PDFs and history.

Before every push, audit the **actual Git index** and build Sphinx with
warnings treated as errors through the publication watchdog. Final PDFs and
verification stdout are evidence. Compiler objects, dependencies, TeX
intermediates, caches and generated sites stay untracked.

The optional [bulk-data manifest](repro/public-data-manifest.json) binds936
scientific files in33 archives. **No bulk asset has been uploaded.** See
[repro/README.md](repro/README.md) and [external input retrieval](repro/EXTERNAL_INPUTS.md)
for availability and the distinction between byte restoration and mathematical
verification. No third-party redistribution licence is implied.
