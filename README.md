# Coxeter Polytopes

Dimension bounds, hyperbolic reflection geometry, exact proof tools and selected
Lean verification. The repository contains **six paper families**, their TeX/PDF
editions, supporting mathematical arguments and reproducibility records.

Start with the [documentation website](https://sashakolpakov.github.io/CoxeterPolytopes/),
the [current bounds and verification guide](docs/current_bounds_and_verification.md),
or the [paper index](manuscripts/README.md).

## Current dimension results

The new bounds below are **internally reviewed written proofs**. Published
results, actual script replay and conditional Lean verification are distinguished
in the [verification guide](docs/current_bounds_and_verification.md).

| Class | Largest published example | Published ceiling | Repository ceiling |
|---|---:|---:|---:|
| General Coxeter, finite volume | 21 | 995 | **57** |
| General compact Coxeter | 8 | 29 | **18** |
| Compact arithmetic Coxeter | 8 | 29 | **16** |
| Noncompact arithmetic Coxeter | 21 | 21 | **21**, published and sharp |
| Right-angled, finite volume, all vertex types | 8 | 12 | **8**, internally sharp |
| Right-angled, entirely ideal | 4 | 6 | **4**, internally sharp |
| Right-angled, compact | 4 | 4 | **4**, published and sharp |

The compact18 and arithmetic16 improvements are supporting research in the
existing compact paper families. Their frozen PDFs still state23 and20.
Arithmetic reflection lattices automatically admit the quadratic-form model;
“first type” is not an extra restriction here. Exact published inputs and
references are in the [assumption register](docs/manuscript_assumptions.md).

At **exact compact dimension18**, only **whole B18 or whole D18** remains
in the accepted necessary frontier. The new complete
[whole-A18 exclusion](results/review_compact18_whole_a18_and_two_classical_group_frontier_root_acceptance_20261023.md)
follows the complete
[D17+actual-I](results/review_compact18_whole_d17_actual_i_and_three_group_frontier_root_acceptance_20261023.md)
and [B17+actual-I](results/review_compact18_whole_b17_actual_i_and_four_group_frontier_root_acceptance_20261023.md)
exclusions. These are configurations that any hypothetical candidate must have
at one selected vertex, not examples. Arithmetic dimension16 still has thirteen
necessary groups. The global compact ceiling stays18.

## The six papers

All six editions are **written, compiled and internally reviewed**. This does
not mean external refereeing, complete script replay or full geometric Lean
verification. The [paper index](manuscripts/README.md) gives the sources, proof
maps and exact verification boundaries.

| Paper | Result in the PDF | PDF |
|---|---|---|
| General finite volume | Dimension at most57 | [48 pages](manuscripts/build/finite_volume57.pdf) |
| Sharp right-angled dimensions | Finite volume at most8; entirely ideal at most4 | [5 pages](manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf) |
| General compact | Dimension at most23; later written work reaches18 | [36 pages](manuscripts/build/compact_farthest_vertex_bounds.pdf) |
| Compact arithmetic | Dimension at most20; later written work reaches16 | [21 pages](manuscripts/build/compact_arithmetic_dimension_bound.pdf) |
| Ideal right-angled four-volume and structure | Published minimum; second volume and unique double; restricted 24-cell chains | [17 pages](manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf) |
| Compact right-angled four-volume constraints | Facet capacity, excess inequalities and recognition results | [14 pages](manuscripts/build/compact_ra4_volume_constraints.pdf) |

The57 paper has12 pages of main exposition and36 pages of essential proof and
algorithm appendices. The other TeX files are preserved predecessors, supporting
notes, fragments and atlases, not additional proposed papers.

## What the methods do

The general finite-volume proof assigns payments to actual polygon corners.
Spherical and affine vertex diagrams limit the available payments, while global
face counting requires a larger total. The compact proofs instead use one
farthest vertex: projections of the **same vector** have additive squared norms,
and their combined required cost exceeds its norm budget. The sharp right-angled
proof uses actual facet projections followed by spherical-code bounds.

The active construction search walks actual edges in explicit original lattices.
The latest [glued nine-dimensional continuation](results/review_glued_e7_a2_heptagon_and_ambient_e2_edge_full_scope_root_acceptance_20261023.md)
closes a compact heptagon face and then leaves it along an additional actual
ambient edge. Finite volume of that chamber is still unresolved: this is **not a
new example**. A matching local diagram or a closed face is insufficient.

The finite-volume57 improvement lane now has a
[complete finite-centre classification and safe discharge union](results/review_finite57_finite_centre_safe_union_and_heavy_escape_root_acceptance_20261023.md),
following complete [isolated finite-cap](results/review_finite57_all_isolated_finite_short_fork_receivers_root_acceptance_20261023.md)
and [ideal-centre](results/review_finite57_ideal_fork_centre_joint_parallel_financing_root_acceptance_20261023.md)
financing results. Costly and shared-capacity cases remain; there is no56 theorem.

## Verification and reproduction

Each algorithm is documented with its mathematical domain, completeness and
termination argument, and the exact result of execution. Hash agreement and PDF
compilation are not mathematical checks. The [computational proof contract](docs/computational_proof_contract.md),
[finite57 reproduction guide](docs/finite_volume57_reproduce.md) and per-paper maps
separate those obligations.

Finite57 has all25 bounded B13/frontier replacement units actually passed and
integrated. The146-job full selection has no combined runtime PASS or complete
transitive regeneration; inherited and newly replayed checks retain their exact
scopes. The frozen paper's proof is internally reviewed at its declared finite
input boundaries.

The newest [A18 scalar kernel](docs/a18_two_dotted_endpoint_scalar_kernel_logic_20261023.md)
has a [bounded compilation and separate ten-declaration axiom audit PASS](results/lean_a18_two_dotted_endpoint_verification_20261023/manifest.json).
It proves continuous ordered-ring inequalities and paid-square logic for two
dotted endpoint opposites. Its real geometric interpretation, complete actual
inventories and legal payment premises remain explicit. All ten prints use only
`propext`, `Classical.choice` and `Quot.sound`; none uses `sorryAx`.
The earlier five-case frontier kernel keeps its original scope and is not
relabelled as verification of the later two-case written frontier.

For named, serial verification batches:

```sh
python3 -B scripts/verify_manuscript.py finite57 --mode core
python3 -B scripts/verify_manuscript.py right_angled_sharp --mode core
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.A18TwoDottedEndpointScalarKernel
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.A18TwoDottedEndpointScalarKernel --audit-axioms
```

See [formal/README.md](formal/README.md) for other modules and their assumptions.
Successful local objects need not be rebuilt without changed inputs.

## Resource and publication discipline

Run one scientific, Lean or TeX child at a time under the shared lock. Admission
requires40GiB free disk and6GiB available RAM. Script/TeX batches use59seconds,
1GiB and one thread; Lean uses2048MiB and one thread, followed by a separate axiom
audit. Unused owned build products are removed after each batch. Start a cleanup
agent only below20GiB free disk. Preserve sources, certificates, final PDFs and
execution history; never kill unrelated processes.

Before pushing, inspect the **actual Git index** for compiler output, caches,
large blobs and secrets, and build Sphinx with warnings treated as errors through
the publication watchdog. Final PDFs and bounded verification stdout are evidence;
Lean objects, TeX intermediates and generated sites are untracked.

The optional [bulk-data manifest](repro/public-data-manifest.json) describes936
scientific files in33 archives. **No bulk asset has been uploaded.** These files
must not be routed through Git as a substitute. See [repro/README.md](repro/README.md)
and [external input retrieval](repro/EXTERNAL_INPUTS.md) for the distinction between
byte restoration and mathematical verification.
