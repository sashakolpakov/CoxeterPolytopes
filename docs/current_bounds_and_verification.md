# Current bounds and verification

Updated23 October2026. There are **six paper families**, each with a written,
compiled and internally reviewed TeX/PDF edition. Later supporting proofs improve
the compact bounds without multiplying papers. This page reports the current
mathematical state; dated receipts and preserved predecessor guides record history.

## Dimension comparison

| Class | Largest published example | Published ceiling | Current internally reviewed ceiling | Relation to the PDF |
|---|---:|---:|---:|---|
| General Coxeter, finite volume | 21 | 995, Prokhorov/Khovanskii | **57** | Current48-page paper |
| General compact Coxeter | 8 | 29, Vinberg | **18** | Supporting proof;36-page PDF states23 |
| Compact arithmetic Coxeter | 8 | 29, Vinberg | **16** | Supporting proof;21-page PDF states20 |
| Noncompact arithmetic Coxeter | 21 | 21, Esselmann | **21**, sharp and published | No new claim |
| Arithmetic Coxeter, compact or noncompact | 21 | 29 overall | **21**, sharp overall | New compact proof plus published noncompact theorem |
| Right-angled, finite volume, all vertex types | 8 | 12, Dufour | **8**, internally sharp | Current5-page paper |
| Right-angled, entirely ideal | 4 | 6, Kolpakov | **4**, internally sharp | Same5-page paper |
| Right-angled, compact | 4 | 4, published and sharp | **4**, unchanged | External starting result |

“Internally reviewed” means a complete written argument with the stated external
and finite-data interfaces, another author's source review and root acceptance.
It is distinct from publication, refereeing, complete computational replay and
end-to-end formalization. References and exact published scopes are recorded in
[the assumption register](manuscript_assumptions.md) and the
[right-angled proof map](right_angled_sharp_dimensions_publication_proof_map_20261021.md).
The quadratic-form model is automatic for arithmetic reflection lattices;
“first type” is not an extra restriction on the arithmetic rows.

## Current compact frontiers

The [global compact18 composition](../results/compact_dimension18_global_composition_and_whole_d19_exclusion_candidate_20261022.md)
excludes dimension19 using the complete same-vertex forest list and all its
branches, together with the previous higher-dimensional tail. Its
[written acceptance](../results/compact_dimension18_global_written_acceptance_20261022.json)
and [conditional Lean execution](../results/lean_compact18_composition_verification_20261022/manifest.json)
retain separate scopes.

At **exact dimension18**, the necessary frontier now contains only:

1. wholeB18;
2. wholeD18.

The complete [B17+actual-I exclusion](../results/review_compact18_whole_b17_actual_i_and_four_group_frontier_root_acceptance_20261023.md)
first reduced five groups to four. The complete
[D17+actual-I exclusion](../results/review_compact18_whole_d17_actual_i_and_three_group_frontier_root_acceptance_20261023.md)
then reduced four to three. The latter includes a fresh rank17 centre inventory,
all foreign contacts, both actual fork signs, adaptive13/14 budgets and one
head/helper/I square. An explicitly preserved interval clarification separates
its small-parameter branch from the joint-point branch.

The subsequent complete
[whole-A18 exclusion](../results/review_compact18_whole_a18_and_two_classical_group_frontier_root_acceptance_20261023.md)
reduces three groups to two. Its complete endpoint classification uses actual
edge-length bounds, explicit affine obstructions and all finite/dotted own
contacts. The last case forces two actual walls to be orthogonal by the native
angle gap, then uses their positive H4+A1 projection and one fixed norm budget.
The proof includes all weak-end patterns and internal holes. It incorporates
the earlier both-isolated theorem rather than assuming endpoint isolation.
This complete geometric result is written and internally reviewed, not
Lean-verified. The remaining two frontiers are necessary configurations,
not examples, and no compact17 theorem follows.

The complete [B18 window/head reduction](../results/review_compact18_whole_b18_fresh_window_head_and_finite_helper_root_acceptance_20261023.md)
and [D18 window/strong5/head reduction](../results/review_compact18_whole_d18_fresh_window_strong_five_head_root_acceptance_20261023.md)
are now internally reviewed supporting theorems. They retain the same original
farthest vector and actual endpoints. The B18 finite-helper obstruction uses
finite-angle conjugates without an arithmetic-lattice hypothesis. Terminal
and dotted-cap couplings remain; neither whole family is excluded.

The [global arithmetic16 composition](../results/review_arithmetic17_whole_a17_global16_complete_composition_root_acceptance_20261022.md)
closes every original exact-dimension17 forest group over every actual admissible
field. Its [three-declaration conditional Lean audit](../results/lean_arithmetic16_composition_verification_20261022/manifest.json)
leaves the geometric, field and computed-data interfaces explicit. The next
[arithmetic16 frontier](../results/review_arithmetic16_universal_a12_five_sevenths_root_acceptance_and_ten_group_frontier_20261023.md)
first had ten necessary groups. The new all-field theorem gives every whole
A12 component energy>5/7, with arbitrary finite spectators. A complete helper,
near4 and endpoint argument handles every remaining weak-foot pattern in one
vector. This closes original groups6/7, after the separate two-defect proof
closed group8; no A12-dominant group remains. The older
[thirteen-group frontier kernel](../results/lean_arithmetic16_frontier_thirteen_group_verification_20261022/manifest.json)
retains that earlier conditional scope. It does not verify these later deletions,
realized examples or arithmetic15.

The subsequent [universal A13 energy theorem](../results/review_arithmetic16_a13_fifty_three_seventieths_and_two_ordinary_alternatives_root_acceptance_20261023.md)
gives e_A13>.757198412>53/70 and removes its ordinaryA3 and ordinaryT2+A1
alternatives. The [direct cubic7-field argument](../results/review_arithmetic16_k7_a13_integer_strip_and_three_a13_alternatives_root_acceptance_20261023.md)
gives e_A13>.777074642 in that actual field and removes actualO7+A1.
It bounds three integral-basis coefficients directly; no census ran.
The subsequent [actual-I endpoint-graft theorem](../results/review_arithmetic16_a13_actual_i_ordinary_dihedral_group_nine_root_acceptance_20261023.md)
excludes the last **A13+actual shortI+ordinaryT2** alternative and deletes
original group9 entirely. It pays both endpoint pressures and every weak foot
in one A13 square; its current ordinaryT2 floor is3/19. There are now
**nine necessary groups**, with original IDs1,2,4,10,11,14,15,18,19.
The global arithmetic ceiling remains16; no A13-dominant group survives.
The [universal A14 theorem](../results/review_arithmetic16_a14_point597_four_fifths_root_acceptance_20261023.md)
now gives e_A14>4/5 with arbitrary finite spectators and every actual field.
It excludes A14+two ordinary nonshortA1. Original group10 still retains
ordinaryT2 and actualN7/O7; the nine original group IDs do not change.

The [A13 scalar logic map](a13_i_endpoint_graft_scalar_kernel_logic_20261023.md)
and source have full different-author review and
[root precompile acceptance](../results/review_a13_i_endpoint_graft_scalar_kernel_root_precompile_acceptance_20261023.md).
The [first actual admission was refused](../results/lean_a13_i_endpoint_graft_verification_20261023/attempt_1_resource_deferral.json)
BEFORE Lean at4.7GiB RAM. The next admitted compile
[failed at two tactic scalings](../results/lean_a13_i_endpoint_graft_verification_20261023/attempt_2_failed_compile.json).
Two explicit identity rewrites have [outside/root clearance](../results/review_a13_endpoint_graft_two_identity_rewrites_root_preclear_20261023.md);
the repaired retry was first [refused before Lean at5.8GiB](../results/lean_a13_i_endpoint_graft_verification_20261023/attempt_3_resource_deferral.json).
A later admitted retry has an [actual compile PASS in1.33seconds](../results/lean_a13_i_endpoint_graft_verification_20261023/pre_audit_manifest.json),
at6.49GiB available RAM. The separate audit was initially refused at5.3GiB.
A later [actual eleven-declaration audit PASS](../results/lean_a13_i_endpoint_graft_verification_20261023/manifest.json)
took1.28seconds at6.34GiB, reusing the successful compilation. All eleven
prints use only standard foundational axioms; no `sorryAx` occurs. Failed diagnostic dependencies,
including error-propagated sorryAx, have zero verification credit.
The source derives continuous scalar gains; actual geometry, inventories,
legal prices and full weak-foot absorption remain external. Its source/map
are immutable preparation records; later actual receipts must supply execution.

The new [continuous A12 kernel](a12_adjusted_loss_kernel_logic_20261023.md)
has complete different-author source/map review. It derives full-interval
adjusted losses and the finite-list paid square from primitive ring data,
including the reference's exact rational value. Its first admission was refused
BEFORE Lean at5.7GiB RAM. The second admitted compile failed on two proof tactics;
its output is preserved. The reviewed repairs now have an
[actual bounded compile and separate ten-declaration audit PASS](../results/lean_a12_adjusted_loss_verification_20261023/manifest.json):
2.37seconds and1.51seconds respectively. All prints use only standard
foundational axioms; no `sorryAx` occurs. The6GiB admission gate is unchanged,
and owned temporary work was removed. The full geometric A12/A13 theorems remain
outside this narrower scalar kernel. Its source header and logic map remain
dated preparation records; the execution receipt supplies the later status.

The [latest B18/D18 reductions](../results/review_compact18_b18_strong17_and_d18_point58_root_acceptance_20261023.md)
are complete reviewed necessary theorems. D18 has a fresh .58 window and both
actual fork cosh lengths<12, or<10 if5 is .58-strong. In a selected strong17
B18 pure-helper/isolated-short configuration the two actual walls had to be
orthogonal. The [complete actual-A17-neighbour argument](../results/review_b18_strong_helper_actual_a17_neighbor_root_acceptance_20261023.md)
now excludes that locus too, closing the selected strong17 branch. The new
[unconditional isolated17 theorem](../results/review_compact18_b18_unconditional_isolated_seventeen_root_acceptance_20261023.md)
excludes every isolated17 own label and the entire dotted interval, without
strong17 or a premise on the actual short18 opposite. The earlier selected branch's
auxiliary actual vertex is not chosen farthest; a weakened comparison keeps
the same original vector. The selected [weak pure dotted centre16 branch in D18](../results/review_compact18_d18_weak_pure_dotted_centre_branch_complete_root_acceptance_20261023.md)
is now excluded entirely. A complete actual fork list forces a far-attached
opposite; all three far profiles then contradict the same funded reference,
including every zero/finite/parallel/dotted mutual case. The new
[finite-own pure-helper theorem](../results/review_compact18_d18_finite_own_weak_pure_centre_root_acceptance_20261023.md)
rebuilds the strong15 and all-head references at finite g, then closes all
actual fork roles and continuous mutual cases. Thus no weak PURE far-left3
centre16 helper remains, whether finite-own or dotted-own. Other centre16
profiles remain open. The earlier selected weak
pure-helper17 branch is [fully excluded](../results/review_compact18_b18_selected_native_five_weak_seventeen_complete_branch_root_acceptance_20261023.md).
Other terminal profiles remain; neither whole family is excluded.

## Verification status

| Result or package | Written status | Executed verification and remaining boundary |
|---|---|---|
| General finite volume57 | Complete paper, internal composition review | All25 bounded B13/frontier replacement units passed and are integrated. Matching inherited evidence covers111 other script selections. No combined146-job full runtime PASS or complete transitive regeneration. Selected counting/table/scalar kernels compiled and separately audited; actual geometry remains external. |
| General compact18 | Complete supporting composition, outside review and root acceptance | Its five-declaration composition kernel compiled and separately audited under23 written inputs. The later two-family frontier is written mathematics; older five-/six-case modules are not relabelled as verifying it. |
| Compact arithmetic16 | Complete supporting composition, outside review and root acceptance | Three actual audited composition declarations under27 explicit written inputs; separate thirteen-group kernel also passed. No full number-field or geometric instance. |
| Right-angled8 and ideal4 | One complete5-page paper, internal reviews and exact checks | Projection/height, tensor-polynomial and count logic have selected script/Lean checks. Actual geometry and the cited facet/code/classification results remain external. |
| Ideal RA4 volume and structure | Complete17-page paper and internal reviews | Exact combinatorial/model checks and selected conditional logic; no complete geometric Lean formalization or external referee approval. |
| Compact RA4 volume constraints | Complete14-page paper and internal reviews | Exact capacity/excess checks and conditional logic. It does not prove a global minimum or a unique second-volume double. |

### Heavy-centre forest logic

The [new kernel map](heavy_centre_forest_scalar_kernel_logic_20261023.md)
has an [actual bounded compilation and separate eleven-declaration audit PASS](../results/lean_heavy_centre_forest_verification_20261023/manifest.json).
Primitive square equations imply the cone inequalities; strong finite-chain
induction derives no-return from two local geometric bridges. Natural variables
index the chain or count incidences; metric quantities remain in an ordered ring.
Actual hyperbolic instances, graph-cycle realization,
forest Euler and reverse ports remain external. All eleven prints use only
standard foundational axioms, with no `sorryAx`. Compilation took2.4seconds
and the separate audit1.93seconds. The earlier5.8GiB admission refusal and
one unused-variable warning are retained. The531912-byte object is ignored;
owned temporary proof work is absent. No older successful module was rerun.

### A18 continuous scalar verification

[A18TwoDottedEndpointScalarKernel](a18_two_dotted_endpoint_scalar_kernel_logic_20261023.md)
now has [actual bounded compilation and a separate ten-declaration audit PASS](../results/lean_a18_two_dotted_endpoint_verification_20261023/manifest.json).
The unchanged12998-byte source derives radical chords from square identities,
the full continuous core gap, endpoint gains, one paid completion and all four
endpoint budget patterns over a genuine ordered ring. Actual geometric pressures,
all-hole references, individual legal prices, orthogonality and a real instance
remain explicit premises. It supports the two-dotted branch, not a wholeA18
exclusion or the later mixed theorem by itself.

Compilation took3.25seconds and the admitted audit1.77seconds. All ten printed
axiom sets contain exactly`propext`, `Classical.choice` and `Quot.sound`; no
`sorryAx` occurs. One unused-variable warning is retained, and no cosmetic rerun
was made. The earlier audit admission refusal is preserved alongside the later
PASS. Owned temporary proof work is absent; its691656-byte object is an ignored
cache, not a publication asset. The linked map and preclear are immutable
preparation records; the execution receipt supplies the later status.

The earlier [eleven-declaration endpoint kernel](general_eighteen_endpoint_refinement_logic_20261023.md)
retains its [actual PASS](../results/lean_general18_endpoint_refinement_verification_20261023/manifest.json)
for the preceding FIVE-case deduction and finite analytic integer seals. Its
scope is not expanded by later source reviews. Printed foundational axioms do
not list theorem parameters; consult each logic map for external assumptions.

## Work toward finite volume56

The global ceiling remains57. New accepted selected-domain results finance:

* [all isolated finite short-fork caps](../results/review_finite57_all_isolated_finite_short_fork_receivers_root_acceptance_20261023.md), with unique actual cheap finite centre receivers;
* [every actual ideal fork-centre neighbour](../results/review_finite57_ideal_fork_centre_joint_parallel_financing_root_acceptance_20261023.md), with all parent families, cusp ports and previous debts;
* the [safe finite-centre union](../results/review_finite57_finite_centre_safe_union_and_heavy_escape_root_acceptance_20261023.md), from the complete twelve-class finite endpoint table.

The last result also proves that a component of heavy vertices cannot be closed
under all three centre moves: its actual three-face would be compact and
right-angled, and three facet colours would force every polygon to have an even
number of sides at least six, contradicting Euler's formula.

The [new metric forest theorem](../results/review_finite57_heavy_centre_forest_and_linear_exit_expansion_root_acceptance_20261023.md)
strengthens this: heavy-to-heavy centre edges have cosh length at least16,
and successive directions are orthogonal, so no cycle is possible. Every
heavy set S has exactly |S|+2c(S) centre-edge incidences leaving it. A full
component with H vertices has H+2 nonheavy exits. This supplies linear
expansion. The later [mixed row1/row2 theorem](../results/review_finite57_row_one_two_variable_cone_component_financing_root_acceptance_20261023.md)
extends the forest to both receiving types, including every finite dihedral
label and arbitrary dotted exchanged products. A variable-angle cone argument
handles the nonright row2 turns. Full-component counting supplies enough safe
eight-unit tokens when all nonheavy centre neighbours have these types,
retaining the old K/Q/L1 payments. Other exits remain explicit in its criterion;
there is no universal364 cap or dimension56 theorem.

The new [zero-exchanged costly-exit theorem](../results/review_finite57_zero_exchanged_bad_centre_affine_financing_root_acceptance_20261023.md)
finances full-table rows4/5/6/11. Positive actual affine kernels recover
old-p edges and all fifteen whole cusp-parent families. Global multiplicities
and canonical budgets retain every old shared load; the minimum joint
remaining balance is5161/2. Among still-unpaid nonheavy finite centre exits
with zero exchanged product only rows7/8 remain. Positive exchanged products
remain open. This written theorem does not prove a universal364 cap or56.

## Construction search

The glued E7/ternary nine-dimensional host is now [proved nonreflective](../results/review_glued_e7_a2_two_unit_releases_original_period_and_literal_nonreflectivity_root_acceptance_20261023.md).
After the actual norm-four corner, two complete23-profile first-wall proofs
give two actual finite vertices and their full outward germs. Each proof keeps
all primitive norm classes, both embeddings, original-module membership and
unknown-wall firstness. Five displayed projection-sign errors in excluded
rows of the first edge were [corrected with byte-preserved history](../results/review_glued_e7_a2_e2_release_five_sign_successor_root_acceptance_20261023.md);
the repaired exclusions and actual first wall survive.

Two individually nonintegral metric involutions have a product G that preserves
the ENTIRE original module, proved by a universal parity identity. It maps the
two actual full germs and preserves the full chamber. Its exact trace5+4phi>10
proves infinite order. Distinct powers give infinitely many cosets of the full
reflection group in the arithmetic lattice, hence infinite chamber volume.
This closes precisely this literal host. It gives no new compact9 example,
field-wide theorem or commensurability-class exclusion.

A different [squared-unit coupling host](../results/review_h9_squared_unit_glued_host_first_two_actual_edges_root_acceptance_20261023.md)
uses the original E7 summand and ternary Gram
`[[2,-1,phi²],[-1,2,0],[phi²,0,-2phi]]`. Complete first-wall proofs retain
all four primitive norm classes. Its [whole E7 face](../results/review_h9_squared_unit_third_edge_full_e7_quadrilateral_root_acceptance_20261023.md)
is a compact quadrilateral with angles(pi/3,pi/2,pi/2,pi/2), area pi/6.
The [first ambient D6 edge](../results/review_h9_squared_unit_first_d6_ambient_edge_root_acceptance_20261023.md)
leaves that face and reaches the full actual germ D6+[2]+[4b]+[2], b=2+phi.
Two [further actual edges](../results/review_h9_squared_unit_d6_b2_to_b3_actual_edges_root_acceptance_20261023.md)
reach D6+B2+[4b] and then a full D6+B3 germ, determinant16. A complete
nine-profile proof handles a tied norm-four mirror by proving it redundant;
the actual new wall has original norm2. This advances the construction; ambient finite volume,
compactness and any new H9 example remain unproved. The subsequent
[long D6-leaf edge](../results/review_h9_squared_unit_d6_long_leaf_edge_original_involution_root_acceptance_20261023.md)
has a complete nine-profile firstness proof and reaches another full D6+B3
vertex. A product of two individually nonintegral reflections is original
on the entire lattice and preserves the full chamber. Its order is exactly
two, with trace6; it proves neither infinite chamber volume nor finite volume.

## Papers, algorithms and resource limits

The [six-paper index](../manuscripts/README.md) maps results to techniques,
source/PDF editions and proof maps. The
[result-manuscript correspondence](result_manuscript_correspondence.md) preserves
more detailed history. Supporting improvements remain in the same paper families;
parked refinements and earlier editions are preserved rather than added as papers.

The [computational proof contract](computational_proof_contract.md) requires each
algorithm's domain, completeness and termination proof, exact arithmetic, output
and replay boundary. Use [finite57's reproduction guide](finite_volume57_reproduce.md)
and [the formal-project guide](../formal/README.md) for named batches. Hashes,
source review, arithmetic replay and geometric formalization are different checks.

One scientific/Lean/TeX child runs at a time under the shared lock. Admission is
40GiB disk/6GiB available RAM; runtime floors are8GiB/3GiB. Script/TeX batches
retain59seconds/1GiB/one thread; Lean2048MiB/one thread, then a SEPARATE audit.
Remove only owned unused products after each batch. Cleanup starts below20GiB
free disk; sources, certificates, final PDFs and execution history are preserved.

Publication requires the actual staged-index debris/secret audit and strict
Sphinx build through the256MiB sampled-RSS publication watchdog. No compiler
objects, caches or generated sites belong in Git. The optional936-file/33-archive
bulk payload remains unuploaded and excluded from Git; see
[repro/README.md](../repro/README.md).

The exact preceding reader guides are [byte-preserved](../results/source_archive/reader_guides_before_classical_frontier_refresh_20261023/manifest.json).
Historical execution receipts retain their original source bindings and are not
silently rebound to this rewritten current guide.
