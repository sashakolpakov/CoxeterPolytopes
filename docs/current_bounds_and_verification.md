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

The [global arithmetic16 composition](../results/review_arithmetic17_whole_a17_global16_complete_composition_root_acceptance_20261022.md)
closes every original exact-dimension17 forest group over every actual admissible
field. Its [three-declaration conditional Lean audit](../results/lean_arithmetic16_composition_verification_20261022/manifest.json)
leaves the geometric, field and computed-data interfaces explicit. The next
[arithmetic16 frontier](../results/review_arithmetic16_whole_d14_root_acceptance_and_thirteen_group_frontier_20261022.md)
still has **thirteen necessary groups**; its separate
[frontier kernel](../results/lean_arithmetic16_frontier_thirteen_group_verification_20261022/manifest.json)
checks that conditional refinement, not realized examples or arithmetic15.

## Verification status

| Result or package | Written status | Executed verification and remaining boundary |
|---|---|---|
| General finite volume57 | Complete paper, internal composition review | All25 bounded B13/frontier replacement units passed and are integrated. Matching inherited evidence covers111 other script selections. No combined146-job full runtime PASS or complete transitive regeneration. Selected counting/table/scalar kernels compiled and separately audited; actual geometry remains external. |
| General compact18 | Complete supporting composition, outside review and root acceptance | Its five-declaration composition kernel compiled and separately audited under23 written inputs. The later two-family frontier is written mathematics; older five-/six-case modules are not relabelled as verifying it. |
| Compact arithmetic16 | Complete supporting composition, outside review and root acceptance | Three actual audited composition declarations under27 explicit written inputs; separate thirteen-group kernel also passed. No full number-field or geometric instance. |
| Right-angled8 and ideal4 | One complete5-page paper, internal reviews and exact checks | Projection/height, tensor-polynomial and count logic have selected script/Lean checks. Actual geometry and the cited facet/code/classification results remain external. |
| Ideal RA4 volume and structure | Complete17-page paper and internal reviews | Exact combinatorial/model checks and selected conditional logic; no complete geometric Lean formalization or external referee approval. |
| Compact RA4 volume constraints | Complete14-page paper and internal reviews | Exact capacity/excess checks and conditional logic. It does not prove a global minimum or a unique second-volume double. |

### New continuous scalar verification

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
number of sides at least six, contradicting Euler's formula. This guarantees an
exit, **not** sufficient receiving capacity. Costly and multiply used finite
exits remain. There is no unproved Hall completion or universal364 cap.

## Construction search

The current glued E7/ternary nine-dimensional host has a sequence of proved
actual ambient edges. The [three-edge B3 theorem](../results/review_glued_e7_a2_three_actual_edges_and_original_b3_break_root_acceptance_20261023.md)
uses complete first-wall profiles, both embeddings and original-module tests.
The [next accepted continuation](../results/review_glued_e7_a2_heptagon_and_ambient_e2_edge_full_scope_root_acceptance_20261023.md)
recognizes that B3 facet as the old-f3 wall, closes a whole compact heptagon of
area3pi/2, and proves a further ambient edge to the old-e2 wall. Its full new
germ is two H3 factors, two norm2 A1s and the retained norm2b line.

The immediate matching metric maps fail original integrality. A closed face or
matching germ does not prove an ambient chamber period. Finite volume,
reflectivity/nonreflectivity and a new compact9 example remain unresolved for
this host. Other literal hosts have separate nonreflectivity results; these do
not exclude whole commensurability classes or arbitrary invariant lattices.

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
