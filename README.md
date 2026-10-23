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
ceiling to **18** and the compact arithmetic Coxeter ceiling to **16**
over every actual admissible field. Those stronger proofs are supporting
research in the same two paper families; they are **not yet theorems of
the PDFs stating 23 and 20**. The [compact eighteen composition](results/compact_dimension18_global_composition_and_whole_d19_exclusion_candidate_20261022.md),
[arithmetic sixteen composition](results/arithmetic17_whole_a17_complete_endpoint_census_and_global16_composition_20261022.md)
and [comparison guide](docs/current_bounds_and_verification.md) give their
precise hypotheses and review records. General finite volume remains 57.

At exact compact dimension18, the [universal A12 energy above13/20](results/review_compact_a12_point65_and_general18_eight_groups_root_acceptance_20261022.md)
closes its entire row without arithmeticity and with arbitrary finite
spectators. The later [complete D16 typed-seven exclusion](results/review_compact18_whole_d16_typed_seven_root_acceptance_and_seven_group_frontier_20261022.md)
also closes that whole row. The [complete B16 exclusion](results/review_compact18_whole_b16_typed_seven_root_acceptance_and_six_group_frontier_20261023.md)
uses an actual neighboring A16 endpoint and one jointly funded vector to
close both typed-seven cases. The new [complete A17 acceptance](results/review_compact18_whole_a17_actual_i_and_five_group_frontier_root_acceptance_20261023.md)
removes both rank-one alternatives; the [D17 ordinary exclusion](results/review_compact18_d17_ordinary_rank_one_root_acceptance_20261023.md)
leaves actual I there. **Five necessary groups remain:** whole A18,
B17+actual I, whole B18, D17+actual I and whole D18. A separate
[analytic endpoint theorem](results/review_compact18_a18_paired_finite_analytic_root_acceptance_20261023.md)
excludes the branch with two isolated finite endpoint opposites in A18;
it does not remove that whole group. The [new logic and arithmetic kernel](docs/general_eighteen_endpoint_refinement_logic_20261023.md)
has a [bounded compile and separate eleven-declaration audit PASS](results/lean_general18_endpoint_refinement_verification_20261023/manifest.json),
with geometric and analytic bridges explicit. The earlier
[six-group kernel](results/lean_general18_six_group_frontier_verification_20261023/manifest.json)
retains its actual PASS. At arithmetic dimension16, the [complete D14 exclusion](results/review_arithmetic16_whole_d14_root_acceptance_and_thirteen_group_frontier_20261022.md)
combines ordinary-dihedral and both actual typed-seven cases, leaving
**thirteen groups**. These are reviewed written refinements, not new global
ceilings or examples.

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
general frontier is dimension eighteen. The arithmetic dimension-seventeen
frontier, including whole A17, is now completely excluded below.

The [complete arithmetic D16+I exclusion](results/arithmetic17_d16_actual_i_complete_exclusion_composition_20261022.md)
now has [full outside review and root acceptance](results/d16_actual_i_exclusion_and_finite57_hexagon_root_acceptance_20261022.md).
The actual-field reduction and separate quadratic, cubic and quartic
contradictions exhaust that entire forest. At that stage its removal left
three whole families; the later D17 exclusion below left two. At that
checkpoint the arithmetic ceiling remained seventeen. The
[compact eighteen logic kernel](docs/compact_dimension_eighteen_composition_logic_20261022.md)
now has [bounded compilation and a separate five-declaration axiom audit PASS](results/lean_compact18_composition_verification_20261022/manifest.json),
with no `sorryAx`. It checks the final composition under23 named written
theorems and supplies no geometric instance. The first failed elaboration
and its repaired, independently reviewed successor are preserved.

The [whole D17 exclusion](results/review_arithmetic17_whole_d17_complete_exclusion_root_acceptance_20261022.md)
now combines complete field reductions with an exact36-field census.
Its120 pairs leave two necessary configurations; both violate one
simultaneously funded metric budget, forcing squared norm greater than
1.00184 where it must be below1. The four actual census units and full
written composition have separate source-bound reviews. This removes a
forest family; at that checkpoint the global arithmetic ceiling stayed17. The
[new scalar Lean kernel](results/lean_classical17_census_verification_20261022/manifest.json)
also passed a bounded compilation and separate four-declaration audit.
It checks the numerical implications under explicit geometric/reference
premises; the enumeration and actual geometry remain written inputs.

The later [whole-B17 exclusion](results/review_arithmetic17_whole_b17_complete_exclusion_root_acceptance_20261022.md)
also has full source-bound internal review and root acceptance. Fresh
head/short-edge geometry, all finite-head cases and actual B17 cycle
normalization put its quadratic branch inside the same executed census;
both survivors violate the funded norm budget. Its incomplete predecessor's
finite-head gap and exact correction are preserved. No successful census
job was repeated. This left only whole A17. Its
[accepted endpoint reduction](results/review_arithmetic17_whole_a17_endpoint_reduction_root_acceptance_20261022.md)
forces both endpoint opposites to be isolated, excludes the golden field
and bounds both edge cosh lengths below2.1. The later
[geometric field reduction](results/review_arithmetic17_whole_a17_four_actual_fields_root_acceptance_20261022.md)
forces one finite endpoint and derives four possible fields. The
[complete quadratic exclusion](results/review_arithmetic17_whole_a17_quadratic_exclusion_root_acceptance_20261022.md)
removes two: all25-field/159-pair census units and the separate full replay
passed with zero survivors. The two remaining cubic fields K7 and K9
are now excluded too: two complete field units and a separate all-record
replay checked **2,432 coefficient triples with zero survivors**. The
[actual evidence](results/a17_finite_cubic_endpoint_runtime_manifest_20261022.json)
and [full runtime review](results/review_a17_finite_cubic_endpoint_actual_runtime_and_conditional_exclusion_foundations_20261022.md)
bind that finite computation. The [complete composition](results/arithmetic17_whole_a17_complete_endpoint_census_and_global16_composition_20261022.md)
and its [full review and root acceptance](results/review_arithmetic17_whole_a17_global16_complete_composition_root_acceptance_20261022.md)
close all fourteen original dimension-seventeen forest groups. Together
with the independently established higher-dimensional tail, this proves
the **global compact arithmetic ceiling16**. The final sixteen composition now has a
[bounded Lean compile and separate three-declaration axiom audit PASS](results/lean_arithmetic16_composition_verification_20261022/manifest.json).
Its [logic map](docs/arithmetic_dimension_sixteen_composition_logic_20261022.md)
exposes27 written geometric/arithmetic/data premises; no geometric instance
is supplied. The first pre-child RAM refusal is preserved, and the successful
compile was reused. The six PDFs remain unchanged.

Fresh work at **exact dimension16** now gives a complete necessary frontier
with **thirteen groups**, following six whole-group exclusions and one
A14 subcase from the original nineteen-group list. The
[accepted refinement](results/review_arithmetic16_whole_a16_and_universal_a3_root_acceptance_20261022.md)
combines whole-A16 exclusion, arithmetic A11 energy above13/20 with arbitrary
spectators, the D13 rank-three and B14 typed-seven exclusions, and the selected
A14+actual-I+ordinary-A1 exclusion. The new compact A3 energy above13/50
requires **no arithmeticity**; together with A11 and the actual typed-seven
defect it removes another group. The later complete D14 proof closes all
three of that row's original spectator cases. These are necessary alternatives, not
realized examples, and do not lower the arithmetic ceiling to15. The
[comparison guide](docs/current_bounds_and_verification.md) lists all thirteen.
A separate [conditional frontier kernel](docs/arithmetic_sixteen_frontier_refinement_logic_20261022.md)
now has a [bounded compile and separate six-declaration axiom audit PASS](results/lean_arithmetic16_frontier_thirteen_group_verification_20261022/manifest.json)
under nine explicit written interfaces, with no `sorryAx`. Its first failed
compile and pre-child RAM refusal remain preserved. The existing
final16 three-declaration Lean PASS above is unchanged.

A new [finite-volume financing theorem](results/finite57_heavy_hexagon_nonadjacent_affine_cusp_global_discharge_candidate_20261022.md)
uses two actual nonadjacent cusps to cover the heavy vertices on every
hexagon in its stated all-heavy domain. Its complete geometry, receiving
multiplicities and prior debts have full internal review and root acceptance.
Other faces and costly forests remain open; no dimension56 bound follows.

The [two-/three-heavy pentagon theorem](results/review_finite57_two_or_three_heavy_pentagon_union_root_acceptance_20261022.md)
also now finances its complete selected domain, using a hyperbolic-line
inequality and bounded actual receiver multiplicities. It replaces the
older three-heavy allocation on overlap. One-heavy pentagons and complete
global coverage remain open; the finite-volume ceiling is still57.

Two further selected cases now have reviewed off-face cusp financing:
[isolated fork3 gives affine E8](results/review_finite57_short_fork_three_off_face_e8_financing_root_acceptance_20261022.md),
and [isolated fork4 gives affine B4](results/review_finite57_short_fork_four_off_face_b4_financing_root_acceptance_20261022.md).
Both use complete source multiplicities and a compatible full8 transfer.
Other contacts and costly forests still lack global coverage.

The new [isolated parallel short-fork theorem](results/review_finite57_isolated_parallel_short_fork_financing_root_acceptance_20261022.md)
finances a further selected heavy-vertex domain at actual off-face ideal
cusps. All seven receiving families have positive remaining curvature,
with minimum3,485 after their full loads. Other short-fork contacts and
costly forests remain open; this does not prove dimension56.

The [merged short-fork matching theorem](results/review_finite57_merged_short_fork_matching_root_acceptance_20261023.md)
adds full8 financing for remaining heavy sources with a private eligible
receiver or at least two merging neighbors. It retains the earlier partial
pentagon rule; complete coverage and dimension56 remain open.

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

The literal norm379 host is also now
[proved nonreflective](results/review_prime379_fourth_germ_original_period_root_acceptance_20261022.md).
Its fourth edge reaches a full E8 corner with an isolated norm2b wall.
Five profile runs and a separate combination checked all15 selectors and294
scalar rows. The original integral map between the two full E8 corners has
trace2947015+4768360phi, proving infinite order and infinite chamber volume.
The proof has complete source-bound mathematical/runtime reviews; the
[execution manifest](results/prime379_fourth_released_r_runtime_manifest_20261022.json)
records1.006seconds total and under61MB maximum sampled process-tree RSS.
This excludes that original lattice only. No new example or field-wide
nonreflectivity theorem follows; the earlier RAM refusal remains preserved.

The new [glued golden E7/A2 module](results/review_glued_e7_a2_original_module_first_two_edges_and_triangle_root_acceptance_20261022.md)
has accepted actual first edges and a complete triangle face. Its
[ambient continuation](results/review_glued_e7_a2_second_ambient_edge_and_quadrilateral_root_acceptance_20261022.md)
reaches full D6+three-A1 and D6+A3 corners and identifies a complete compact
quadrilateral, with its intrinsic angles computed by projection. The later
[actual U edge](results/review_glued_e7_a2_released_e2_actual_a5_bridge_root_acceptance_20261023.md)
and [actual F edge](results/review_glued_e7_a2_v_e2_actual_d4_edge_root_acceptance_20261023.md)
reach full H3+A5+A1 and H3+D4+A1+[2b] corners. Complete first-wall
proofs show the proposed later mirrors T and P2 are nonactual. Finite
nine-dimensional chamber volume remains unresolved. The accepted
[complete H3+D4 face](results/review_glued_e7_a2_complete_right_hexagon_root_acceptance_20261023.md)
is a compact right-angled hexagon of area pi. Its finite face symmetries
are not ambient lattice periods; no new example is claimed.

The distinct determinant-norm551=19*29 host is now also
[proved nonreflective](results/review_composite551_binary_release_and_literal_nonreflectivity_root_acceptance_20261022.md).
Its [seven actual profile runs and separate combine](results/composite551_transverse3_actual_runtime_manifest_20261022.json)
reach a full E7+two-A1 germ. The subsequent complete manual next-edge proof
reaches E8 with the norm2b spectator and constructs an original integral
chamber period of trace843+1354phi>10. Infinite order and the arithmetic
index/covolume argument exclude finite reflection covolume for this literal
module. The manual continuation earns no new worker or Lean execution
credit. Other invariant lattices and the commensurability class are not
excluded; no new example is claimed.

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
