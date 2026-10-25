# Coxeter Polytopes

Dimension bounds, explicit hyperbolic reflection chambers and right-angled
volume problems. The repository has **six paper families**, their frozen
TeX/PDF editions, supporting proofs, reproducibility records and selected
Lean verification.

Start with the [documentation website](https://sashakolpakov.github.io/CoxeterPolytopes/),
the [current bounds and verification guide](docs/current_bounds_and_verification.md)
or the [paper index](manuscripts/README.md).

Current work is the [six-paper audit and consolidation](docs/audit_and_consolidation_status_20261024.md).
That guide separates fresh mathematical audits, manuscript revisions and the
new bounded execution entry. Its 14 infrastructure tests passed; necessary
proof jobs run serially under the monitored resource policy below.
The [right-angled release records](docs/right_angled_release_packages_20261025.md)
now have [current root acceptance](results/right_angled_current_release_root_acceptance_20261025.json)
for both papers, with external inputs and historical execution boundaries explicit.
The compact17 and arithmetic15 working articles now have complete
[mathematical and editorial source acceptance](results/compact17_arithmetic15_final_working_source_root_acceptance_20261025.json).
Their [compact](docs/compact17_consolidated_proof_map_20261025.md) and
[arithmetic](docs/compact_arithmetic_consolidation_proof_map_20261025.md)
maps expose the common method, case divisions and assumed inputs.
The [completed editorial revision](results/compact17_arithmetic15_editorial_completion_20261025.json)
now has matching PDFs, separate source reviews and sampled rendered-page review.
The two new compact Lean modules retain actual compilation and separate actual
axiom audits. The six frozen editions remain unchanged.

## Current dimension results

The repository improvements are **internally reviewed written proofs**:
complete arguments at their declared input boundaries, another author's
source review and coordinating acceptance. Publication, external refereeing,
script replay and geometric formalization have separate statuses.

| Class | Largest published example | Published ceiling | Repository ceiling |
|---|---:|---:|---:|
| General Coxeter, finite volume | 21 | 995 | **57** |
| General compact Coxeter | 8 | 29 | **17** |
| Compact arithmetic Coxeter | 8 | 29 | **15** |
| Noncompact arithmetic Coxeter | 21 | 21 | **21**, published and sharp |
| Right-angled, finite volume, all vertex types | 8 | 12 | **8**, internally sharp |
| Right-angled, entirely ideal | 4 | 6 | **4**, internally sharp |
| Right-angled, compact | 4 | 4 | **4**, published and sharp |

Arithmetic reflection lattices automatically admit the quadratic-form
model; “first type” adds no restriction to these arithmetic rows. Exact
published inputs and scopes are in the [assumption register](docs/manuscript_assumptions.md).
The general bounds assume no arithmeticity: their restrictions come from
physical geometry, native angles and spherical/affine links. Field,
integrality and nonphysical-conjugate conditions enter only the arithmetic lane.
There is **no new polytope example** from the active construction search.

## The six papers

All six editions are **written, compiled and internally reviewed**. Later
compact17/arithmetic15 proofs belong in the existing compact families;
the frozen PDFs still state23/20. Earlier notes, fragments and atlases
are supporting material, rather than additional proposed papers.

| Paper | Result in the PDF | PDF |
|---|---|---|
| General finite volume | Dimension at most57 | [48 pages](manuscripts/build/finite_volume57.pdf) |
| Sharp right-angled dimensions | Finite volume at most8; entirely ideal at most4 | [5 pages](manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf) |
| General compact | Dimension at most23; supporting proof17 | [36 pages](manuscripts/build/compact_farthest_vertex_bounds.pdf) |
| Compact arithmetic | Dimension at most20; supporting proof15 | [21 pages](manuscripts/build/compact_arithmetic_dimension_bound.pdf) |
| Ideal right-angled four-volume and structure | Published minimum; second volume and unique double; restricted24-cell chains | [17 pages](manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf) |
| Compact right-angled four-volume constraints | Facet capacity, excess and recognition | [14 pages](manuscripts/build/compact_ra4_volume_constraints.pdf) |

The57 paper has 12 pages of main exposition and 36 pages of essential proof
and algorithm appendices. The [paper index](manuscripts/README.md) and
[result correspondence](docs/result_manuscript_correspondence.md) identify
each technique, source, proof map and verification boundary. The ideal
four-volume paper does not classify all ideal polytopes as24-cell chains.
The compact four-volume paper does not establish a120-cell minimum or
a unique second-volume double.

Three later working verification copies now have matching bounded builds:
[finite57,50 pages](manuscripts/build/finite_volume57_consolidated.pdf),
[compact17,285 pages](manuscripts/build/compact_farthest_vertex_bounds_consolidated.pdf)
and [arithmetic15,216 pages](manuscripts/build/compact_arithmetic_dimension_bound_consolidated.pdf).
The compact revisions have no overfull boxes, undefined references or package
warnings; nine/eight underfull spacing warnings remain. Their main expositions
occupy about19/38 pages, followed by complete case appendices. These remain
revisions in the same families; the six frozen editions are preserved. The
[paper index](manuscripts/README.md) records exact builds and sampled visual coverage.

## What the proofs do

The finite-volume proof assigns weights to actual polygon corners. Finite
vertex diagrams and affine cusp links bound the available payment; global
incidence counting requires more. Every transfer uses one simultaneous rule
with the full earlier debt retained.

The compact proofs use one farthest vertex and one spatial vector. Positive
norm certificates and orthogonal component projections force a cost above
its budget. Arithmeticity adds actual-field, integrality and angle restrictions.
The sharp right-angled proof uses actual facet projections and spherical codes.

The [compact17 proof map](docs/compact17_consolidated_proof_map_20261025.md)
shows how the component estimates split into chains, classical prefixes and
actual-rank forest consumers. The [arithmetic15 map](docs/compact_arithmetic_consolidation_proof_map_20261025.md)
adds field restrictions and identifies its four computed premises and the
literal literature assumption. Both complete working articles have
mathematical and editorial source acceptance and matching working PDF builds.
The two new [conditional compact Lean kernels](docs/compact_dimension17_conditional_logic_map_20261024.md)
now have actual compilation and separate actual axiom audits for ten declarations.
Geometry, analytic estimates and complete inventories remain explicit premises.
Older successful kernels retain their scopes and were not rerun.

The active work is verification and publication of these six families.
Research toward56/16/14, new examples and the spectral approach is parked;
its [detailed record](docs/current_bounds_and_verification.md) remains available.
No partial forest list is a realized example or a stronger dimension bound.

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

Finite57's [current-rule checks](results/finite57_current_rule_release_completion_20261025.json)
now pass for all 64 shards: 28,028 full maps, including 12,518 nonunit
obligations and 83 prefixes. The [final 259-unit reader](results/finite57_final_current_rule_conjunction_execution_20261025.json)
also passed at the [explicit source and geometric premises](results/finite57_old_affine_longer_tail_root_acceptance_20261025.json).
The [reproduction guide](docs/finite_volume57_reproduce.md) explains exactly
what was checked, which inputs remain assumed, and how earlier failures are preserved.
This is source-bound current-rule verification; native enumeration and the
complete geometric proof have not been regenerated or formalized end to end.

Selected Lean kernels have actual compilation and **separate actual axiom
audits**. Their logic maps expose the ordinary theorem parameters that
printed foundational axioms do not enumerate.

The packages below have **actual bounded compilation and separate
actual axiom audits**. Their [logic maps](formal/README.md) expose the field, census,
geometry and real-instance premises that remain external. The B16 predecessor's
failed compilation and all unused-variable warnings are preserved.

| Kernel | Actual audited declarations | Scope and receipt |
|---|---:|---|
| Quadratic dyadic triples | 14 | [Discrete valuation/parity and norm deductions](results/lean_quadratic_dyadic_triple_verification_20261023/manifest.json) |
| B16 isolated-pair scalars | 12 | [Cubic and two paid-prefix contradictions](results/lean_b16_isolated_pair_scalar_verification_20261023/manifest.json) |
| B18 terminal pressure | 10 | [Continuous scalar pressure contradictions](results/lean_b18_terminal_pressure_verification_20261023/manifest.json) |
| A13 endpoint graft | 11 | [Continuous gains and paid-square logic](results/lean_a13_i_endpoint_graft_verification_20261023/manifest.json) |
| A12 adjusted loss | 10 | [Continuous loss and finite-list paid square](results/lean_a12_adjusted_loss_verification_20261023/manifest.json) |
| A18 dotted endpoints | 10 | [Continuous endpoint inequalities](results/lean_a18_two_dotted_endpoint_verification_20261023/manifest.json) |
| Heavy-centre forest | 11 | [Cone, conditional chain and incidence algebra](results/lean_heavy_centre_forest_verification_20261023/manifest.json) |
| D18 dotted-head scalars and compact17 composition | 10 across two modules | [Conditional scalar and discrete deductions](results/compact17_targeted_lean_completion_20261025/manifest.json) |

These are conditional internal kernels. Actual geometry, complete inventories,
number fields, enumerations and primitive pressures remain explicit inputs.
Earlier failures and resource deferrals earn no PASS credit. The
[formal guide](formal/README.md) records each package's exact coverage.

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
The user's25October policy admits new batches at **2% free disk and 15%
available RAM**, measured against total filesystem capacity and physical RAM.
Stop active jobs below 1% disk or6% RAM; start cleanup below 2% disk. Gate on
exact byte ratios and record totals, percentages and the policy source hash.
Script/TeX batches use59 seconds/1 GiB/one thread; Lean uses2048 MiB/one thread
and a separate axiom audit. Remove only owned unused products after each
batch, preserving sources, certificates, final PDFs, evidence and history.

Before every push, audit the **actual Git index** and build Sphinx with
warnings treated as errors through the publication watchdog. Final PDFs and
verification stdout are evidence. Compiler objects, dependencies, TeX
intermediates, caches and generated sites stay untracked.

The optional [bulk-data manifest](repro/public-data-manifest.json) binds936
scientific files in 33 archives. **No bulk asset has been uploaded.** See
[repro/README.md](repro/README.md) and [external input retrieval](repro/EXTERNAL_INPUTS.md)
for availability and the distinction between byte restoration and mathematical
verification. No third-party redistribution licence is implied.
