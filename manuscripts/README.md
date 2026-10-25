# The six papers and their techniques

There are **six current paper families**, with written, compiled and internally
reviewed TeX/PDF editions. The compact17 and arithmetic15 supporting improvements
belong in the existing compact families when integrated. They are not extra
papers, and the frozen PDFs still state23 and 20. Earlier TeX notes, fragments,
predecessors and generated atlases are preserved supporting material.

The [current audit and consolidation guide](../docs/audit_and_consolidation_status_20261024.md)
tracks fresh reviews and same-family revisions. Exploratory work toward
56/16/14 stays outside these TeX arguments. The delivered pairs below remain
frozen; later working builds do not automatically replace accepted editions.

The compact17 and arithmetic15 working articles now have complete
[mathematical and editorial source acceptance](../results/compact17_arithmetic15_final_working_source_root_acceptance_20261025.json).
Their [general compact map](../docs/compact17_consolidated_proof_map_20261025.md)
and [arithmetic map](../docs/compact_arithmetic_consolidation_proof_map_20261025.md)
show the common foundation, actual branch points and external/computed
premises. Their [editorial completion record](../results/compact17_arithmetic15_editorial_completion_20261025.json)
binds matching PDF builds, separate source reviews and sampled rendered-page
inspection. Complete case appendices remain; the frozen editions below stay23 and20.

| Paper and TeX | Main technique | Result in this PDF | PDF | Proof/reproduction map | Verification ID |
|---|---|---|---|---|---|
| [General finite volume](finite_volume57.tex) | Corner payments, Coxeter forests and cusp curvature | Dimension at most57 | [48 pages](build/finite_volume57.pdf),12 main +36 essential appendices | [Proof map](../docs/finite_volume57_lemma_map.md) | `finite57` |
| [Sharp right-angled dimensions](right_angled_sharp_dimension_bounds_publication.tex) | Actual projections, descending heights and spherical codes | Finite volume at most8; entirely ideal at most4; both internally sharp | [5 pages](build/right_angled_sharp_dimension_bounds_publication.pdf) | [Proof map](../docs/right_angled_sharp_dimensions_publication_proof_map_20261021.md) | `right_angled_sharp` |
| [General compact](compact_farthest_vertex_bounds.tex) | One farthest vertex and additive component energies | Dimension at most23; later supporting proof17 | [36 pages](build/compact_farthest_vertex_bounds.pdf) | [Proof map](../docs/compact_farthest_vertex_bounds_proof_map_20261020.md) | `compact_farthest_bounds` |
| [Compact arithmetic](compact_arithmetic_dimension_bound.tex) | Field reduction and exact angle rules, with the same geometric budget | Dimension at most20; later supporting proof15 | [21 pages](build/compact_arithmetic_dimension_bound.pdf) | [Proof map](../docs/compact_arithmetic20_proof_map.md) | `compact20` |
| [Ideal RA4 volume and structure](ideal_ra4_volume_and_structure_publication.tex) | Volume gap, equality recognition and restricted decomposition | Published minimum; second volume and unique double; restricted24-cell chains | [17 pages](build/ideal_ra4_volume_and_structure_publication.pdf) | [Proof map](../docs/ideal_ra4_volume_and_structure_proof_map.md) | `ideal_ra4_volume_and_structure` |
| [Compact RA4 volume constraints](compact_ra4_volume_constraints.tex) | Facet capacity, excess and recognition | Necessary volume/structure constraints | [14 pages](build/compact_ra4_volume_constraints.pdf) | [Proof map](../docs/compact_ra4_volume_constraints_reproduce_20261019.md) | `compact_ra4` |

The57 and sharp right-angled papers are the priority entry points. The
[current comparison guide](../docs/current_bounds_and_verification.md) distinguishes
published bounds, supporting improvements and actual verification. The
[detailed result correspondence](../docs/result_manuscript_correspondence.md)
preserves the longer provenance. The
[six-paper delivery receipt](../results/six_paper_final_tex_delivery_20261021.json)
binds the matching frozen sources, PDFs and build records.

## Working verification copies

These are revisions within the existing three families. Compilation confirms
a matching PDF; mathematical source review and conditional formal verification
have their own records. The compact revisions now have reviewed presentation
changes and sampled visual acceptance. Their long actual-rank case appendices
retain complete arguments. No frozen paper is silently replaced.

| Working source | Result | Matching PDF | Bounded build receipt and remaining warnings |
|---|---|---|---|
| [General finite volume](finite_volume57_consolidated.tex) | Dimension at most57 | [50 pages](build/finite_volume57_consolidated.pdf) | [Actual build](../results/finite57_terminal_tail_working_tex_build_20261025.json);21 underfull-box warnings; no overfull boxes or undefined references |
| [General compact](revisions/compact_farthest_vertex_bounds_consolidated.tex) | Dimension at most17 | [285 pages](build/compact_farthest_vertex_bounds_consolidated.pdf) | [Actual build](../results/compact17_editorial_layout_final_build_20261025.json); nine underfull spacing warnings; no overfull boxes, undefined references or package warnings |
| [Compact arithmetic](revisions/compact_arithmetic_dimension_bound_consolidated.tex) | Dimension at most15 | [216 pages](build/compact_arithmetic_dimension_bound_consolidated.pdf) | [Actual build](../results/arithmetic15_editorial_publication_style_build_20261025.json); eight underfull spacing warnings; no overfull boxes, undefined references or package warnings |

All three builds report no undefined references and removal of their owned
temporary build files. The compact17/arithmetic15 completion record lists the
12/6 final-PDF pages visually inspected; this is sampled layout review. Their
common expositions end around pages19/38, with genuine case-specific proofs
following. Earlier receipts retain their original scope and dated filenames.
Earlier failed attempts and byte-preserved TeX repair preimages remain evidence.
The current finite57 source includes the reviewed terminal-tail lemma. Its
[pinned replay predecessor](../results/finite57_terminal_tail_working_source_replay_separation_20261025.json),
`revisions/finite_volume57_consolidated.tex` at `0a01a796…`, stays unchanged
because existing ordinary-premise and S8 receipts bind that source. These
are versions within one paper family.

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
does not prove that the 120-cell is the absolute volume minimum, nor that its
double is the unique second-volume example; those remain open here.

## Supporting research and actual checks

The later compact17 and arithmetic15 arguments are complete reviewed written
supporting proofs. Their frozen paper editions remain23 and 20. Consolidate
future improvements within those families rather than create a paper per
numerical increment. The [result correspondence](../docs/result_manuscript_correspondence.md)
links the individual sources, reviews, accepted interfaces and verification.

The [general compact17 theorem](../results/compact_dimension17_global_written_acceptance_20261024.json)
now excludes wholeD18, the only necessary alternative at exact 18. The prior
global18 tail excludes higher dimensions. The final dotted-head argument
uses an actual head point, four complete initial references, BOTH individual
anchor prices and a fresh local whole-A18 lemma. All use the SAME original
farthest vector and every actual wall; no arithmeticity is assumed.
Its [full root acceptance](../results/review_compact_dimension17_whole_d18_and_global_composition_root_acceptance_20261024.md)
records the complete different-author review and 21 bound inputs. The
[two conditional Lean sources](../docs/compact_dimension17_conditional_logic_map_20261024.md)
have complete static reviews and now have
[actual bounded compilation and separate actual axiom audits](../results/compact17_targeted_lean_completion_20261025/manifest.json)
for ten declarations. These verify conditional scalar and discrete logic;
geometry, calculus, classification and enumeration remain explicit inputs.
The first scalar compilation failure, caused by a reserved binder name,
and its reviewed three-token rename are preserved. The earlier resource
deferrals remain historical records. A fresh exact 17
forest analysis now has a [complete reviewed necessary cover](../results/review_compact_dimension17_fresh_complete_necessary_cover_root_acceptance_20261024.md)
with thirteen conservative dominant rows, retaining H3/H4 under BOTH
typed-seven roles. The later whole-family
[A11/A12](../results/review_compact17_whole_a11_a12_fresh_four_foot_root_acceptance_20261024.md)
and [A14](../results/review_compact17_whole_a14_point79_all_feet_root_acceptance_20261024.md)
exclusions remove three rows. The
[fresh classical-tail proof](../results/review_compact17_whole_d14_and_b15_d15_rank_two_tails_root_acceptance_20261024.md)
also closes wholeD14 and every B15/D15 rank-two spectator with
energy>3/19. The later
[complete actual-typed-seven proofs](../results/review_compact17_whole_b15_d15_typed_seven_and_complete_rows_root_acceptance_20261024.md)
remove wholeB15/D15. The
[complete A15 non-typed theorem](../results/review_compact17_a15_all_non_typed_rank_two_root_acceptance_20261024.md)
and [complete actual-typed continuation](../results/review_compact17_whole_a15_both_actual_typed_seven_and_six_row_cover_root_acceptance_20261024.md)
now remove ENTIRE A15. These are complete original-vector
arguments, including all isolated feet and paid helper corrections.
The [A16 rank-one theorem](../results/review_compact17_a16_rank_one_point9_ordinary_a1_exclusion_root_acceptance_20261024.md)
first excluded ordinary A1. The
[corrected whole-A16 theorem](../results/review_compact17_whole_a16_actual_i_corrected_and_five_row_cover_root_acceptance_20261024.md)
now excludes ACTUAL I too, leaving FIVE dominant rows:
A17,B16,B17,D16,D17. Its complete628-line proof and 337-line outside
review preserve the predecessor's explicit own3/WHITE4 coverage omission
and sole repair, with every helper/head cross and individual price checked.
The [whole-A17 dotted-pair theorem](../results/review_compact17_whole_a17_two_isolated_dotted_one_third_root_acceptance_20261024.md)
closes its two isolated dotted endpoints with one raw reference and
a uniform scalar certificate. The
[finite-pair theorem](../results/review_compact17_a17_both_isolated_finite_analytic_root_acceptance_20261024.md)
also closes both isolated finite endpoints by complete monotone native
gaps and rational angle certificates. Mixed finite/dotted and nonisolated
endpoints remain; wholeA17 is not excluded.
No16 ceiling or new example follows.
Arithmetic15 now retains ten necessary groups. The
[complete A13+ordinary-A2 proof](../results/review_arithmetic15_a13_ordinary_a2_whole_original12_root_acceptance_20261024.md)
closes original12 when combined with the accepted old-label-at-least4
and two-ordinary-A1 theorems. It retains all foreign contacts and uses
one ALL-foot reference, a common actual triangle and Cartan determinants.
These field and integrality arguments use arithmeticity; the general
compact exclusions do not. No arithmetic14 bound follows.
The [new arithmetic common-cap alphabet](../results/review_arithmetic_large_a_common_end_cap_trace_alphabet_root_acceptance_20261024.md)
uses a local Lorentzian principal to identify the ACTUAL field, then an
exact trace bound and native-angle tests. Five fields survive at its
stated cap interface; this is not a global field list or a new exclusion
among the ten remaining groups. The
[uniform A12 theorem](../results/review_arithmetic15_uniform_a12_three_quarters_and_whole_original9_root_acceptance_20261024.md)
proves component energy>3/4 with both endpoints and ALL remaining feet
paid in one square, closing original9. The later
[complete original10 proof](../results/review_arithmetic15_whole_original10_prime_thirteen_two_endpoints_and_eleven_group_root_acceptance_20261024.md)
also closes ACTUAL I+T2 by a prime13 obstruction and a single fully
paid endpoint/foot square. I keeps its joint cost and is never assigned
ordinary A1 cost. This arithmetic argument does not enter the general lane.
The later [complete original13 theorem](../results/review_arithmetic15_original13_prime_seven_and_ten_group_root_acceptance_20261024.md)
uses the same norm-divisibility mechanism at prime7, with a fresh
two-A1 endpoint catalogue. Its complete403-line proof and 279-line
outside review leave ten necessary groups, still at ceiling15.

Toward finite volume56, all selected weighted ideal budgets and all native6
finite budgets are accepted. Actual generation and separate replay cover
6,626 small-factor ideal cases. The new polygon-refund method also finances
complete centre-exit rows5 and 11 for literal2D17+D23. Both positive far-row6
rank cases now have complete refund/receiver proofs too. ALL positive
near-row6 exits are now financed at both ranks: exceptions use full
actual receiving budgets, while regular arrivals pay for themselves with
distinct receiver-corner refunds. The original dual-index error and a
later omitted old-V23 debt are both preserved with their exact repairs.
An all-actual-facet extension restricts finite native labels to 2/3/4/6
at dimension57. ALL positive NEAR end3 row4 exits are now financed at both
ranks, sharing one old-plus-new exception reserve and refunding distinct
regular receiver corners. The corrected raw A15 charge362 is retained.
The complementary [complete FAR-end3 proof](../results/review_finite57_far_end_three_and_whole_row_four_root_acceptance_20261024.md)
now closes ENTIRE row4. Quadrilateral sources transfer8 with no corner
refund; regular cases use actual larger polygons. All 36 ideal parent
templates and old debts fit one shared8N union. The
[complete positive-NEAR8 theorem](../results/review_finite57_positive_near_eight_and_entire_row_eight_root_acceptance_20261024.md)
now closes ENTIRE row8 too, together with accepted zero-NEAR8/all-FAR8.
Four paired actual corner refunds preserve the old B22 receiver unit;
all E8 exceptions share ONE full8N reserve. The corrected
[positive NEAR7 proof](../results/review_finite57_all_near_row_seven_corrected_root_acceptance_20261024.md)
now covers both ranks too. Its full D5/A11 charge includes the two
omitted ordinary cross units; the failed predecessor is preserved.
The [per-port composition theorem](../results/review_finite57_all_centre_port_activation_near_seven_and_conditional_l12_root_acceptance_20261024.md)
allows unused covered ports at already-paid sources. The new
[complete FAR7 proof](../results/review_finite57_all_positive_far_seven_and_full_literal_heavy_l12_root_acceptance_20261024.md)
closes the remaining positive row7 ports and therefore EVERY full L12
component of the LITERAL2D17+D23 domain. Correct short duals exclude
the apparent A14/A20 cycle; full receiving unions retain all old debts.
Other costly forests, zero-arrival A16 receivers and complete face/source
coverage remain open; no56 theorem follows.
The new [selected A16+A1 port theorem](../results/review_finite57_a16_ideal_orthogonal_ports_selected_financing_and_residual_root_acceptance_20261024.md)
finances every such vertex with an ideal next-to-end port or orthogonal
exchange. Actual-edge unions keep the old ideal allowance; six new
whole cusp families have complete strict budgets. Both finite port
endpoints and the twelve necessary finite templates remain explicit.
The later [all-covered-port theorem](../results/review_finite57_a16_all_covered_ports_joint_reserves_eight_templates_root_acceptance_20261024.md)
adds seven finite-safe receiving forests and upgrades EVERY covered port
to 8, keeping old5.5 tags and adding only2.5. Eight finite residual
templates remain. Its aggregate16g>=11H condition is sufficient but
unproved universally; no56 theorem or new census follows.

The compact spectral work proves a full-wall Gram band and a fourth-order
transverse obstruction. Its unknown-return and bounded central-prefix work
has no numerical bound or execution credit. Construction research has one
new actual original H9 facet and E6+B2+A1 endpoint after four authenticated
units, followed by a manually proved J facet and full E6+H3 endpoint.
On a separate branch, complete closed strips identify the new actual R7
facet with D7+2A1 endpoint; an old-U edge reaches D6+3A1. All four direct
matches to the earlier such germ are nonoriginal. This is no period or
volume theorem. A fresh mixed direction proves the 50-entry recorded
outer relaxation still has infinite volume; unknown walls may cut it.
The [next short-leaf edge](../results/review_h9_prime41_p7_short_leaf_actual_e7_return_root_acceptance_20261024.md)
now has unique actual R* and full E7+2A1 endpoint. All eleven necessary
first-wall profiles are covered by original parities and closed strips
at both embeddings, including a redundant same-time root. The metric
return is a nonoriginal order-two reflection. The
[next original norm4 facet B*](../results/review_h9_returned_e7_actual_dyadic_edge_return_words_and_outer_gap_root_acceptance_20261024.md)
has a fresh complete firstness proof and full E7+[4]+[2b41] endpoint.
The enlarged54-entry outer relaxation still has infinite volume; every
nonzero alternating power of the natural returns is nonoriginal.
The [N-release theorem](../results/review_h9_prime41_n_release_order_upper_and_complete_nonunit_firstness_root_acceptance_20261024.md)
proves a finite actual first endpoint, excludes every unit first profile
and reduces it to TWO nonunit conics. The
[bounded norm-one upper worker](../results/review_h9_prime41_n_release_norm_one_upper_worker_root_preexecution_clearance_20261024.md)
has full static code/contract review but is UNEXECUTED below RAM admission.
It supplies no explicit unit or first facet. Full-chamber volume remains open.
The prepared Y-release worker was not run or used. Ambient finite volume and
a new example remain open. These are supporting
methods and research notes, not additional paper editions.

Finite57's25 bounded replacement units have actual PASS and catalogue
integration. Its full136-script/nine-Lean-root/one-PDF selection has no combined
146-job replay or complete transitive regeneration; the 111 other command
matches retain inherited execution scope. The core11 selection is unchanged.
See the [reproduction guide](../docs/finite_volume57_reproduce.md).

Finite57 now has [all 64 actual selected applications](../results/finite57_longer_tail_all64_actual_application_completion_20261025.json)
and a [successful 259-unit current-rule reader](../results/finite57_final_current_rule_conjunction_execution_20261025.json),
at the accepted ordinary and supplemental premises. The
[reproduction guide](../docs/finite_volume57_reproduce.md) records the source
identities, coverage and remaining external inputs. This is separate from
the historical146-job catalogue and complete transitive regeneration.

Selected Lean packages have compilation and SEPARATE actual axiom audits,
including the new [fourteen dyadic declarations](../results/lean_quadratic_dyadic_triple_verification_20261023/manifest.json)
and [twelve repaired B16 scalar declarations](../results/lean_b16_isolated_pair_scalar_verification_20261023/manifest.json),
plus the [ten compact17 declarations in two modules](../results/compact17_targeted_lean_completion_20261025/manifest.json).
Their [logic maps and receipts](../formal/README.md) retain the external geometric,
field, classification and enumeration premises. This is not a complete geometric
formalization. A compiled PDF, a source review and an executed replay remain
different evidence. All six TeX/PDF bindings stay unchanged.
