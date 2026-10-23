# The six papers and their techniques

There are **six current paper families**, with written, compiled and internally
reviewed TeX/PDF editions. The compact18 and arithmetic16 supporting improvements
belong in the existing compact families when integrated. They are not extra
papers, and the frozen PDFs still state23 and20. Earlier TeX notes, fragments,
predecessors and generated atlases are preserved supporting material.

| Paper and TeX | Main technique | Result in this PDF | PDF | Proof/reproduction map | Verification ID |
|---|---|---|---|---|---|
| [General finite volume](finite_volume57.tex) | Corner payments, Coxeter forests and cusp curvature | Dimension at most57 | [48 pages](build/finite_volume57.pdf),12 main +36 essential appendices | [Proof map](../docs/finite_volume57_lemma_map.md) | `finite57` |
| [Sharp right-angled dimensions](right_angled_sharp_dimension_bounds_publication.tex) | Actual projections, descending heights and spherical codes | Finite volume at most8; entirely ideal at most4; both internally sharp | [5 pages](build/right_angled_sharp_dimension_bounds_publication.pdf) | [Proof map](../docs/right_angled_sharp_dimensions_publication_proof_map_20261021.md) | `right_angled_sharp` |
| [General compact](compact_farthest_vertex_bounds.tex) | One farthest vertex and additive component energies | Dimension at most23; later supporting proof18 | [36 pages](build/compact_farthest_vertex_bounds.pdf) | [Proof map](../docs/compact_farthest_vertex_bounds_proof_map_20261020.md) | `compact_farthest_bounds` |
| [Compact arithmetic](compact_arithmetic_dimension_bound.tex) | Field reduction and exact angle rules, with the same geometric budget | Dimension at most20; later supporting proof16 | [21 pages](build/compact_arithmetic_dimension_bound.pdf) | [Proof map](../docs/compact_arithmetic20_proof_map.md) | `compact20` |
| [Ideal RA4 volume and structure](ideal_ra4_volume_and_structure_publication.tex) | Volume gap, equality recognition and restricted decomposition | Published minimum; second volume and unique double; restricted24-cell chains | [17 pages](build/ideal_ra4_volume_and_structure_publication.pdf) | [Proof map](../docs/ideal_ra4_volume_and_structure_proof_map.md) | `ideal_ra4_volume_and_structure` |
| [Compact RA4 volume constraints](compact_ra4_volume_constraints.tex) | Facet capacity, excess and recognition | Necessary volume/structure constraints | [14 pages](build/compact_ra4_volume_constraints.pdf) | [Proof map](../docs/compact_ra4_volume_constraints_reproduce_20261019.md) | `compact_ra4` |

The57 and sharp right-angled papers are the priority entry points. The
[current comparison guide](../docs/current_bounds_and_verification.md) distinguishes
published bounds, supporting improvements and actual verification. The
[detailed result correspondence](../docs/result_manuscript_correspondence.md)
preserves the longer provenance. The
[six-paper delivery receipt](../results/six_paper_final_tex_delivery_20261021.json)
binds the matching frozen sources, PDFs and build records.

## Where the common arguments split

Finite57 first establishes one simultaneous corner-payment rule. Finite vertices
then require exact spherical-forest maxima; ideal vertices require affine-factor
charges and actual cusp multiplicities. Both branches reunite in one incidence
contradiction. The higher-dimensional tail is separate and does not assume a
facet is itself Coxeter.

The compact proofs start with the same original farthest vector. Ordinary,
large-A and large-B/D component estimates have separate geometric arguments,
then reunite through additive squared norms. Rank-specific endpoint and fork
cases close the surviving configurations. Arithmeticity adds field/angle
restrictions rather than a second unrelated geometric framework.

The sharp right-angled paper shares one actual projection lemma between its
general and entirely ideal cases. Spherical-code estimates and the relevant
published facet minima enter at clearly identified later steps. Its internally
sharp8/4 bounds do not assume arithmeticity.

The ideal four-dimensional paper combines the published minimum, second-volume
recognition and restricted24-cell chains as one sequence. It does not classify
arbitrary ideal right-angled polytopes as24-cell garlands. The compact RA4 paper
does not prove that the120-cell is the absolute volume minimum, nor that its
double is the unique second-volume example; those remain open here.

## Supporting research and actual checks

General compact18 and arithmetic16 have complete reviewed supporting compositions
and separately compiled conditional Lean modules. At exact compact18 the latest
accepted necessary frontier is wholeB18/D18 after the
[complete whole-A18 exclusion](../results/review_compact18_whole_a18_and_two_classical_group_frontier_root_acceptance_20261023.md).
Arithmetic16 retains [nine original groups](../results/review_arithmetic16_a13_actual_i_ordinary_dihedral_group_nine_root_acceptance_20261023.md).
Its entire A13 row is now excluded, after universal and cubic7-field arguments
and a simultaneous actual-I endpoint/weak-foot square. The cubic step uses
a direct integral-basis proof, not enumeration. Neither list asserts examples.
The PDFs have not been revised for these later bounds. The new compact
terminal reductions and literal-host nonreflectivity belong to supporting
research, not additional proposed manuscripts.

The new mixed row1/row2 financing theorem belongs to the finite-volume family.
The selected strong17 B18 branch, unconditional isolated17 exclusion and
both finite/dotted weak pure D18 centre16 exclusions belong to the general
compact family. The zero-exchanged costly-exit financing theorem belongs
to the finite-volume family. The squared-unit host's actual edges, compact
E7 face, full D6+B3 corners and original order-two chamber symmetry are
construction research; ambient finite volume remains open.
The universal arithmetic A14>4/5 energy theorem and its two-ordinary-A1
consumer belong to the same arithmetic family; the nine group IDs remain.

Each computational algorithm needs a mathematical specification: its complete
domain, a termination and completeness argument, exact output, and the part
actually replayed. The [computational proof contract](../docs/computational_proof_contract.md)
and individual proof maps give those boundaries. A source hash or a successful
PDF build proves neither coverage nor a mathematical theorem.

Finite57 has all25 bounded B13/frontier replacement units passed and integrated.
The full146-job selection, including136 script selections, nine Lean roots and
one PDF build, has no combined runtime PASS or complete transitive regeneration.
Historical and bounded replacements retain their individual execution records.

The [A12 continuous-loss scalar kernel](../docs/a12_adjusted_loss_kernel_logic_20261023.md)
has [bounded compilation and a separate ten-declaration audit PASS](../results/lean_a12_adjusted_loss_verification_20261023/manifest.json).
It verifies the conditional loss and paid-square argument; actual geometry,
field inventories and the later universal A12/A13 theorems remain external.
Earlier failed/refused attempts are preserved with zero verification credit.

The [A13 endpoint-graft scalar source/map](../docs/a13_i_endpoint_graft_scalar_kernel_logic_20261023.md)
have complete static review and root preclear. A failed admitted compile and
two before-child RAM refusals are preserved. Two explicit identity rewrites
have outside/root clearance and an actual1.33second compile PASS. The separate
audit was initially refused at5.3GiB RAM and now has an
[actual eleven-declaration PASS in1.28seconds](../results/lean_a13_i_endpoint_graft_verification_20261023/manifest.json),
with no successful compilation rerun and no `sorryAx`. Failed diagnostic
prints earn zero credit. This is not full geometric
verification of the A13 result.

The newest [A18 continuous scalar kernel](../docs/a18_two_dotted_endpoint_scalar_kernel_logic_20261023.md)
has an [actual bounded compile and separate ten-declaration audit PASS](../results/lean_a18_two_dotted_endpoint_verification_20261023/manifest.json).
Actual geometry, legal prices and a real instance remain explicit. It verifies
the two-dotted scalar branch, not wholeA18 or the newer two-family written
frontier. The previous five-case kernel retains its original scope. See
[formal/README.md](../formal/README.md) for all conditional logic boundaries.

Run named checks through the existing coordinator/runner, for example:

```sh
python3 -B scripts/verify_manuscript.py finite57 --mode core
python3 -B scripts/verify_manuscript.py right_angled_sharp --mode core
```

One scientific/Lean/TeX child runs at a time. Admission requires40GiB disk and
6GiB available RAM; the existing time, memory and thread caps remain. Remove
owned unused products after every batch, preserving mathematical evidence.
The preceding reader guide is [byte-preserved](../results/source_archive/reader_guides_before_classical_frontier_refresh_20261023/manifest.json).
