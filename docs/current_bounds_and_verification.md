# Dimension bounds and verification

Updated 22 October 2026. There are six TeX papers, each **written,
compiled and internally reviewed**, with a matching PDF. This manuscript
status does not imply complete replay or full geometric Lean verification.
Parked arithmetic refinements are excluded and preserved.
Fresh internal mathematical audits found no blocker at the stated
published-input and finite-data boundaries. This is not external
refereeing or complete geometric formalization.

The [reader’s index](../manuscripts/README.md) gives the six PDFs and
reproduction commands. The [result–manuscript map](result_manuscript_correspondence.md)
separates the theorem printed in each PDF from later supporting research.
The [delivery receipt](../results/six_paper_final_tex_delivery_20261021.json)
binds the finished editions; the
[completed-audit checkpoint](../results/six_paper_completed_audits_and_restricted_continuations_20261021.json)
records those reviews. The [global-bound checkpoint](../results/compact19_arithmetic17_written_acceptance_20261021.json)
records compact nineteen and arithmetic seventeen; the
[previous research checkpoint](../results/compact19_arithmetic17_fresh_frontiers_and_paid_endpoint_acceptance_20261021.json)
supplies the fresh necessary forests. The
[previous A-family checkpoint](../results/a_family_and_prime79_frontier_continuations_acceptance_20261021.json)
accepts further exclusions and five construction-model eliminations. The
[previous endpoint checkpoint](../results/endpoint_case_exclusions_and_prime79_interval_acceptance_20261021.json)
adds reviewed A17/A18 and arithmetic A15/A16 exclusions, with full root
and other-author mathematical reviews, plus one actually executed norm79
interval. The [previous acceptance](../results/a14_a18_and_prime79_period_acceptance_20261021.json)
adds the last general A18 companion exclusion, all arithmetic A14 cases,
and literal norm79 nonreflectivity with separate script and Lean records.
The [previous acceptance](../results/b18_a11_and_restricted_continuations_acceptance_20261021.json)
excludes general B18+actual-I and arithmetic A11+actual-I+(A5 or D5),
leaving five and seven necessary groups at that checkpoint. The
[previous acceptance](../results/d18_exclusion_b16_reduction_and_prime239_first_edge_acceptance_20261021.json)
also excludes general D18+actual-I, narrows the arithmetic B16+I short
opposite, and accepts a restricted quadrilateral theorem and the literal
norm239 first binary edge, leaving four and seven necessary groups.
The [previous acceptance](../results/d17_typed7_cubic_b16_and_prime239_transverse_acceptance_20261022.json)
excludes general D17 with either actual typed-seven defect and the
entire cubic short-cap branch of arithmetic B16+I, and accepts the next
norm239 transverse edge with full E7+I2(5) germ. The next
[accepted batch](../results/quadratic_b16_b19_and_rank_transfer_acceptance_20261022.json)
closes the quadratic branch and hence the whole B16+I row. The
[previous acceptance](../results/b15_b19_short4_and_prime239_e7_a2_acceptance_20261022.json)
also excludes B15 with either actual typed-seven defect, leaving three
general and five arithmetic necessary groups at that checkpoint. It narrows B19 to one
last-long wall, accepts the third norm239 edge and records eleven actual
Lean declarations in two new conditional kernels. The
[latest acceptance](../results/whole_b19_exclusion_acceptance_20261022.json)
closes the last B19 wall and hence the whole B19 forest, leaving only
whole A19 and D19 in general dimension nineteen. These batches do not
lower a global ceiling. The six-paper writing milestone is complete;
research has resumed within the existing paper families.

## Dimension comparison

A ceiling is an upper bound on dimension, not an existence claim.
The repository columns report internally reviewed written arguments.

| Class | Largest known example | Published ceiling | Repository written ceiling | Finished PDF |
| --- | ---: | ---: | ---: | --- |
| Arbitrary finite-volume Coxeter | 21 | 995, Prokhorov | 57 | [57](../manuscripts/build/finite_volume57.pdf), 48 pages |
| Arbitrary compact Coxeter | 8 | 29, Vinberg | 19 | [23](../manuscripts/build/compact_farthest_vertex_bounds.pdf), 36 pages |
| Compact arithmetic Coxeter | 8 | 29, Vinberg | 17 | [20](../manuscripts/build/compact_arithmetic_dimension_bound.pdf), 21 pages |
| Noncompact arithmetic Coxeter | 21 | 21, Esselmann | 21, published and sharp | No separate new paper |
| Arithmetic Coxeter, compact or noncompact | 21 | 29 overall | 21, sharp as a union | Compact argument plus published noncompact theorem |
| Compact right-angled | 4 | 4, Potyagailo–Vinberg | 4, published and sharp | Background in the right-angled paper |
| Entirely ideal right-angled | 4 | 6, Kolpakov | 4, sharp | [4](../manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf), same 5-page paper |
| Right-angled with mixed finite/ideal vertices | 8 | 12, Dufour | 8, sharp | [8](../manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf), same 5-page paper |

The universal right-angled theorem permits all vertex types and proves
dimension at most eight without arithmeticity. Its entirely ideal branch
proves four. Compact right-angled four was already sharp in the literature.
Noncompact arithmetic dimension twenty is excluded by the published
classification; its sharp ceiling twenty-one must not be confused with
the compact arithmetic problem.

The admissible quadratic-form model is automatic for arithmetic hyperbolic
reflection lattices: the terms first type and simplest type impose no
additional restriction here. The
[external-input register](manuscript_assumptions.md#arithmetic-reflection-groups-and-the-quadratic-form-model)
states the published reduction used to enter that model. Historical
source labels retain their original wording and hash bindings.

The new [global composition](../results/compact_dimension19_and_arithmetic_dimension17_global_composition_20261021.md)
proves compact nineteen and arithmetic Coxeter seventeen. It has
full root verification and two outside-author composition reviews:
[foundations](../results/review_compact19_arithmetic17_global_composition_foundations_reviewer_20261021.md)
and [second review](../results/review_compact_dimension19_arithmetic17_global_composition_second_full_scope_20261021.md).
Every new branch also has a full source-bound mathematical review.
The arithmetic proof covers every actual admissible ground field,
with no prescribed field list. These are internally reviewed written
proofs, with their geometric and published inputs explicit. The finished
compact PDFs still state twenty-three and twenty; the new nineteen/
seventeen results belong in those same two paper families when integrated.

## What has been checked

| Paper | Written audit | Reproduction and formal boundary |
| --- | --- | --- |
| Finite volume 57 | [Fresh full composition audit](../results/review_finite57_geometric_composition_new_audit_20261021.md); no blocker at named finite-input boundaries | Matching bounded records cover 111/113 selected script commands. The reviewed 25-unit replacement for two eager commands has **3 PASS, 22 without PASS**. There is no combined 123-job full replay or complete transitive regeneration. Conditional Lean checks selected counting, table and scalar implications. |
| Sharp right-angled 8/ideal 4 | Fresh [geometry](../results/review_sharp_ra_geometric_foundations_adversarial_20261021.md) and [tensor/equality](../results/review_sharp_ra_tensor_ideal_equality_adversarial_20261021.md) audits | Exact core checks and selected conditional Lean have actual successful runs. Coordinate projection, height and cubic identities are formalized; published bounds and actual geometric interfaces remain explicit premises. |
| Compact arithmetic PDF 20 | Full exposition/field-closure reviews and [completed internal-field proof review](../results/review_compact_arithmetic20_internal_field_translation_foundations_20261021.md) | Its matching 21-page build has zero warnings. Accounting and finite checks keep their declared scopes; this is not formal verification of the later seventeen theorem. |
| General compact PDF 23 | [Full source review](../results/review_compact_farthest23_full_batch_guard_20261020.md), with source-bound editorial delta review | Matching 36-page build has zero warnings. The compiled accounting kernel retains twenty-three; actual geometric inequalities and later nineteen remain outside that kernel. |
| Ideal RA4 volume/structure | [Fresh full adversarial audit](../results/review_ideal_ra4_volume_structure_fresh_adversarial_20261021.md) | Universal second-volume equality is proved before the restricted chain theorem is applied. Exact model outputs and conditional Lean have explicit scopes; the actual cutting and recognition proofs are written geometry. |
| Compact RA4 constraints | [Fresh full adversarial audit](../results/review_compact_ra4_volume_constraints_fresh_adversarial_20261021.md) | Graph and geometric cases are written out. The catalogue selection checks typesetting only; there is no universal geometric Lean proof. |

The finite57 [reproduction guide](finite_volume57_reproduce.md) records
all completed bounded replacements and the still-unintegrated loader.
Its immutable [partial successor ledger](../results/finite57_loader_replay_quadratic_b16_20261022.json)
preserves the parser and first two B13 pin-group PASS records, including
the second group's earlier resource refusal before its later successful run.
Resume with `b13:pins-2` in a new ledger; never
replace this partial receipt with a completion claim.

Recent actual Lean additions include the
[two-foot payment kernel](two_foot_paid_logic_20261021.md) and
[large-parent arithmetic kernel](large_inert_parent_arithmetic_logic_20261021.md),
with separate compile and actual axiom-audit records. The
[regular 24-cell height kernel](regular_ideal_24cell_height_logic_20261021.md)
now has bounded compilation and a separate audit of **five actual declarations**.
The [classical composition kernel](classical_dimension_composition_logic_20261021.md)
also compiled and passed its separate **four-declaration** audit for the
conditional nineteen/seventeen logic and quartic norm contradiction.
The [completion receipt](../results/lean_publication_completion_20261021/manifest.json)
binds their unchanged sources, four serial execution ledgers and actual output.
Only standard Lean axioms appear; none uses `sorryAx`. Two Regular24
unused-section-variable warnings are preserved; Classical has none.
The earlier [admission refusal](../results/regular24_height_compile_admission_refusal_20261021.json)
remains historical evidence. Each map separates derived logic from actual
geometry, external/script inputs and field-norm bridges. Those premises are
not discharged by an ordinary axiom print or by this conditional compilation.

The new [endpoint triangle kernel](endpoint_triangle_norm_kernel_logic_20261021.md)
also has bounded compilation and a separate actual ten-declaration audit
PASS. Its [receipt](../results/lean_endpoint_triangle_verification_20261021/manifest.json)
records the successful final source and both failed tactic predecessors,
which earn no credit. It derives cycle cancellation, the integral Bezout
expression and norm-product logic; actual fields, cyclic integrality and
endpoint geometry remain explicit application premises.

The [Lorentz obstruction kernel](lorentz_obstruction_kernel_logic_20261021.md)
now has bounded compilation and a separate actual seven-declaration axiom
audit [PASS](../results/lean_lorentz_obstruction_verification_20261021/manifest.json).
Only standard axioms occur, with no `sorryAx`; failed predecessors retain
zero credit. It derives the second-orbit cutoff, the golden comparison,
and restricted large-D Schur/rank implications. Actual averaging/Cauchy
realization, full wall germs, original-lattice preservation, arithmetic
covolume and Gram/rank interpretation remain written application inputs.

Two further small kernels now have actual bounded compilation and separate
axiom-audit PASS records: [D18 centre scalars](../results/lean_d18_centre_verification_20261022/manifest.json)
has five declarations, and [cubic B16 arithmetic/energy logic](../results/lean_cubic_b16_short_cap_verification_20261022/manifest.json)
has six. Only standard axioms occur; neither uses `sorryAx`. Each has one
preserved unused-variable warning. Their unchanged source headers and logic
maps retain their original preparation status; these execution receipts
record the later successful verification. Actual geometry, arithmetic-field
and exhaustive-cap applicability remain explicit premises.

## Recent dimension proofs and next work

The previous general dimension-twenty frontier is now closed. Its whole
A20 and B20 exclusions are joined by a complete whole-D20 proof. Actual
long-head isolation leaves native five, finite labels at least seven,
or a dotted head. The corrected [native-five proof](../results/compact_whole_d20_native_five_adaptive_endpoint_repair_20261021.md),
[new finite-head proof](../results/compact_whole_d20_all_finite_large_heads_paid_core_exclusion_20261021.md)
and [new dotted-head proof](../results/compact_whole_d20_dotted_head_adaptive_paid_reference_exclusion_20261021.md)
exclude all three branches. The two new proofs retain all possible left
helpers, select the actual strong endpoints that exist, and pay every
nonzero head/helper cross term. The earlier unsupported individual tail
floors remain explicitly superseded; their historical files are preserved.
Combined with the existing higher-dimensional bound, this gives nineteen.

The previous arithmetic dimension-eighteen frontier is also closed.
The actual A17+I, B17+I and D17+I exclusions leave whole D18.
Its finite long heads were already excluded; the remaining dotted heads
have degree at most four. New complete [quadratic](../results/arithmetic18_whole_d18_quadratic_dotted_head_exclusion_20261021.md),
[cubic](../results/arithmetic18_whole_d18_cubic_dotted_head_two_cycle_exclusion_20261021.md)
and [quartic](../results/arithmetic18_whole_d18_dotted_quartic_exclusion_by_two_cycle_20261021.md)
proofs exclude every remaining field degree. The rational case uses the
explicit external anisotropy criterion. The key improvement is the actual
two-cycle divisibility $UV/4\in\mathcal O_K$. The quartic contradiction
is a short norm inequality; the cubic proof ends with seven hand-derived
polynomial rows; the quadratic proof covers all 49 derived head polynomials.
There is no unexecuted enumeration assumed as a result.

The current frontier dimensions are nineteen and seventeen. There are now
fresh complete necessary lists for [general compact nineteen](../results/compact_dimension19_complete_forest_frontier_candidate_20261021.md)
and [arithmetic seventeen](../results/arithmetic_dimension17_complete_forest_frontier_candidate_20261021.md).
Both start again from the actual farthest-vertex forest and enumerate the
small spectators; neither decrements the closed D20/D18 lists. The general
original list has nine grouped rows, and the conservative arithmetic list fourteen.
I denotes the actual isolated short A1, and the two typed rank-two defects
retain their different old/new meanings. These are necessary lists, not
classifications of realized polytopes.

A new rank-free isolated-I lemma controls every other actual edge, giving
joint I/A1 energy greater than one fifth. A further
[paid-endpoint proof](../results/arithmetic17_a15_isolated_i_and_rank_one_exclusion_20261021.md)
excludes two alternatives in the conservative arithmetic list: A15+I+A1
and A12+I+T3+A1, where T3 is an **irreducible** ordinary rank-three
component. The difficult A15 case uses tightened rational central floors
and one additional endpoint pressure on the same vector. Its exact norm
lower bound is
`1+340306992/1800000000000>1`. I is counted once.

The new [general A16 proof](../results/compact19_a_family_continuation_candidate_20261021.md)
gives component energy greater than .796 with arbitrary spectators and
without arithmeticity. Its complete smaller-rank weak-cap analysis
includes a native-five/native-seven algebraic argument. This excludes
all four A16 spectator alternatives at dimension nineteen.

The new [lower-A arithmetic continuation](../results/arithmetic17_lower_a11_a12_continuation_candidate_20261021.md)
uses the previously reviewed A12 energy greater than .6725 to remove
every A12 alternative at seventeen. The smallest exact margin is 3/2800.
Its fresh N12 cubic-field argument also removes A11+either typed7+T4;
the conjugate interval changes to (-.85,.5), and every isolated foot is
paid simultaneously in one reference. Finally, a new
[A15 half-loss proof](../results/arithmetic17_a15_typed_seven_half_loss_field_strip_exclusion_20261021.md)
uses three explicit integral-basis exclusions in the actual cubic-seven
field to remove both A15+typed7 forests. The held all-field A15 theorem
also removes two ordinary A1 spectators. A15 at that checkpoint retained only
an ordinary irreducible T2. None of these exclusions assumes a universal
ordinary-dihedral energy greater than one sixth.

The latest [endpoint proofs](../results/review_endpoint_continuations_and_prime79_runtime_root_20261021.md)
exclude that ordinary-T2 arithmetic A15 case, both ordinary and short
A1 beside arithmetic A16, and every general A17 rank-two companion.
The general A18 component bound greater than .91 removes its ordinary A1
companion. A subsequent [two-end paid-foot proof](../results/compact19_a18_isolated_i_two_end_paid_foot_exclusion_candidate_20261021.md)
also excludes its actual-I companion, closing the entire A18 row.
The arithmetic A16 extension first proves that its small-energy rank-one
opposite is isolated, then extends the triangle norm argument to its
larger height interval. It covers every actual field, including the rational
case, without using Meyer. These supporting proofs belong in the existing
two compact paper families and have separate full source-bound reviews.

The fresh [arithmetic A14 proof](../results/arithmetic17_a14_all_spectator_energy_and_three_family_exclusions_20261021.md)
proves component energy greater than .76 with arbitrary finite spectators
over every actual field, and a stronger bound in the actual-I case.
It excludes all three remaining A14 spectator alternatives. Its same-vector
couplings do not assume that every ordinary T2 has energy greater than 1/6.

The new [B18+actual-I proof](../results/compact19_b18_isolated_i_full_endpoint_foot_exclusion_candidate_20261021.md)
closes its entire general row through a fresh short-tail classification
and paid head/I pressures. The new
[A11+actual-I proof](../results/arithmetic17_a11_actual_i_quadratic_time_modulus_exclusion_20261021.md)
keeps all three near4 classes that survive at this smaller rank. Actual
quadratic-field two-cycles then force an integral norm both below eight
and divisible by 144. Together with its simultaneous foot estimate this
excludes both remaining A11 companions over every actual field. Both
proofs have full root and outside-author reviews in the latest receipt.

The complete reduced necessary lists are now:

| General compact nineteen: dominant | Entire remaining spectator alternatives |
| --- | --- |
| A19 | none |
| D19 | none |

| Arithmetic Coxeter seventeen: dominant | Entire remaining spectator alternatives |
| --- | --- |
| A17 | none |
| B17 | none |
| D15 | either typed7 |
| D16 | ordinary A1 or actual I |
| D17 | none |

Here T2 and T3 are irreducible ordinary components, and I/typed7 retain
their actual-wall definitions. There are two and five grouped rows,
respectively. These necessary lists make no realization claim. The
[latest checkpoint](../results/whole_b19_exclusion_acceptance_20261022.json)
binds the full reviews and subtraction of exactly whole B19. The previous
batches removed arithmetic B15+either typed7 and B16+I. Earlier batches removed general
D18+I and D17+either typed7. All other necessary groups and the six
printed editions are preserved.
The bounds still leave gaps of eleven and nine above the known compact
example dimension eight. No new polytope example has been found.

The [complete D18+I proof](../results/compact19_d18_isolated_i_complete_centre_and_head_exclusion_candidate_20261021.md)
extends the earlier local tail proof through every center and head case.
A continuous all-label inequality and a shifted mixed-I payment close
the tight center rows. Every point uses actual D18 normals, including
both forks; no virtual short normal is introduced. The row is now closed.

The [complete D17+typed7 proof](../results/compact19_d17_typed_seven_fresh_centre_exclusion_candidate_20261021.md)
recomputes every tail and centre case at rank seventeen. Its two actual
rank-two defects have different foreign-contact lists, but both pay
energy greater than .09. Pure component pressures and one final paid
reference exclude both alternatives. Whole A19 and D19 remain after
the later B19 exclusion below.

For the now-excluded arithmetic B16+I row, the
[short-opposite reduction](../results/arithmetic17_b16_actual_i_short_opposite_isolation_and_degree_restriction_20261021.md)
excludes all six nonisolated white roles. A surviving short cap must be
isolated and dotted. The new
[cubic polynomial argument](../results/arithmetic17_b16_isolated_cubic_cap_unique_polynomial_candidate_20261021.md)
forces a unique cubic field and integral basis. The complete
[actual-I strip and energy proof](../results/arithmetic17_b16_isolated_cubic_cap_actual_i_strip_exclusion_candidate_20261021.md)
then excludes that entire cubic branch. The complete
[quadratic proof](../results/arithmetic17_b16_quadratic_short_cap_actual_i_and_head_exclusion_candidate_20261022.md)
derives all fourteen short-cap and sixteen actual-I possibilities,
excludes positive foreign contacts and treats every actual long head.
Its minimum final norm estimate is1.00026. Full root and outside-author
reviews accept the whole-row exclusion. No relative I-root length class
or unexecuted field census is assumed.

The [B15+typed7 proof](../results/arithmetic17_b15_typed_seven_long_head_exclusion_candidate_20261022.md)
instead avoids the short opposite. The actual defect forces the cubic-seven
field and energy greater than.09; a fresh prefix reference and complete
pure/mixed long-head analysis give total norm at least1.00069125. Both
typed defects are excluded, with their different actual-wall definitions
retained. This has full root and other-author reviews.

The [complete B19 two-wall reduction](../results/compact_whole_b19_two_terminal_wall_reduction_foundations_20261022.md)
and its [short4 exclusion](../results/compact_whole_b19_short_four_terminal_fresh_index17_exclusion_foundations_20261022.md)
first leave only last-long own3/short3, with an isolated native5 long head.
The latter proof audits every cap at one newly improved index17 floor,
pools the terminal slopes to preserve concavity, and uses the same cap's
actual short19 endpoint. Its single-vector total is1.00017135845.
The [final short3 proof](../results/compact_whole_b19_short_three_terminal_fresh_index17_exclusion_foundations_20261022.md)
audits every cap at a new single-index loss.50. It retains the possible
orthogonal dotted cap, uses actual mutual-product gaps for every other
role, and pools the last slopes before lifting the endpoint. One original
reference gives total1.00111946. Full root and other-author reviews
accept the corrected proof and its conjunction with the two-wall reduction.
Both alternatives are excluded, so the whole B19 group is now removed.
This leaves whole A19 and D19, without lowering the global ceiling19.

For finite volume, the [complete Context662 composition](../results/finite57_complete_context662_source_composition_20261021.md)
finances this one selected producer over all whole-parent domains.
Its inequality is against held actual ideal charge after retained
priorities and assigned debt. The [current A17 provenance audit](../results/finite57_current_a17_endpoint_producer_scope_audit_20261021.md)
finds only the two already-financed endpoint producers at the named
complete written interfaces. It is not a fresh all-current target replay,
whole-coupon deletion or a global fifty-six proof. The formal forest
2D17+D23 still has charge 372, independently of A17; a new global budget
argument would be necessary. No new polytope example has been found.
The new [homogeneous large-D quadrilateral obstruction](../results/finite57_homogeneous_large_d_quadrilateral_schur_obstruction_20261021.md)
excludes its precise full-white, unchanged-spectator configuration by a
Schur/rank argument. It supplies no uniform charge saving or lower ceiling.
A new [electrical proof for long affine A cycles](../results/finite57_long_affine_a_distance_five_zero_transfer_candidate_20261021.md)
reproves the earlier held order-at-least-26 distance-four statement.
It uses a cycle-Laplacian gap estimate and the native height gap; its
order25 local Lambert quadrilateral explains the obstruction to simply
extending the deletion. This is an alternative proof of an existing
internal theorem, not a new bound or an accepted change to the frozen
payment tables. It leaves the other affine families and the costly
finite forests unresolved.

The new [full D16 source-orbit proof](../results/finite57_d16_full_orbit_priority_with_retained_context662_credits_20261021.md)
pays all 136 pairs in the varying D16 factor at its specified whole cusp
types, while retaining the previous selected finite-recipient cap.
Its complete reserved surpluses are 74 and 280. This is a restricted
replacement allocation, with no reduction of the formal maximum 372
or the dimension ceiling 57.

A [new actual-quadrilateral restriction](../results/finite57_heavy_rank_transfer_quadrilateral_local_discharge_20261021.md)
shows that the specified long-leaf transfer between two charge372
forests forces a third corner of charge at most357.125. The
[stronger six-orientation proof](../results/finite57_rank_transfer_all_active_attachments_gap_exclusion_candidate_20261022.md)
lowers that charge to353.640625 under exactly the same actual-face
hypotheses, with complete root and other-author reviews. A global
routing and multiplicity bound is still missing; the maximum372 remains.

The [literal norm79 host](../results/compact_h9_nonunit_coupling_prime79_full_base_and_first_transverse_20261021.md)
is now proved nonreflective at its accepted actual-germ and admissible-arithmetic
boundaries. The earlier 38,784-trial interval found the asymmetric A9-7
corner. The subsequent [END0 execution](../results/prime79_actual_a9_end0_interval_execution_20261021.json)
completed all 810,436 scalar trials and six first-simple profiles: six
integral mirrors, with unique minimum D9-right id4. Its
[domain](../results/compact_h9_prime79_actual_a9_end_edge_upper_and_complete_contract_candidate_20261021.md)
and [algorithm contract](../results/compact_h9_prime79_actual_a9_end_interval_checker_spec_20261021.md)
have full source reviews.

The [two-D9 proof](../results/compact_h9_prime79_two_actual_d9_germ_period_candidate_20261021.md)
matches full actual germs and proves original-lattice integrality of an
infinite-order chamber symmetry. A separate bounded
[two-map exact check](../results/prime79_two_d9_germ_period_execution_20261021.json)
passed all matrix/inverse, Gram, glue and orbit identities; it is not another
root census. The [acceptance receipt](../results/a14_a18_and_prime79_period_acceptance_20261021.json)
binds both actual runs and full written reviews. END8 was refused before a
child and is unnecessary for this conclusion. No new example, field-wide,
form-class, other-lattice or infinite-family exclusion follows.

The [good-prime odd golden-glue extension](../results/golden_odd_glue_full_base_and_discriminant_reduction_20261021.md)
proves that a broader literal lattice family has the entire actual
Bugaenko base. An all-divisor collar and the explicitly inherited analytic
upper bound give a finite necessary host reduction, which was not run.
The distinct norm239 lattice has this starting geometry. Its
[first binary edge](../results/compact_h9_prime239_full_base_and_first_binary_edge_candidate_20261021.md)
is proved by complete six/five-row hand enumeration, with unique norm2
mirror and full finite E8+A1 germ. The
[next transverse edge](../results/compact_h9_prime239_first_e7_golden_transverse_edge_candidate_20261021.md)
now has a complete five-profile first-simple-wall proof and full finite
E7+I2(5) endpoint. The [third-edge proof](../results/compact_h9_prime239_second_h2_edge_first_wall_and_e7_a2_germ_candidate_20261022.md)
identifies a unique original norm2 mirror and full E7+A2 endpoint.
Its complete finite algorithm and actual325,636 scalar trials are
documented as mathematical statements. All twenty units and the separate
combine have actual bounded PASS records in the
[runtime manifest](../results/prime239_folded_second_edge_runtime_manifest_20261022.json).
The largest unit took6.579s and used57,950,208B. The full new firstness
and geometric composition have root and outside-author reviews.
Original dyadic glue is retained throughout.
Further chamber completion and finite volume remain open; this is not
a new example.

## Ideal and compact right-angled four-volume

| Statement | Result and limit |
| --- | --- |
| Ideal minimum | Published: 4π²/3, uniquely the regular ideal 24-cell. |
| Second ideal volume | Internal reviewed proof: 8π²/3, uniquely its facet reflection double, with face vector (42,168,164,38). |
| Restricted structure | Octahedral/doubled-octahedral facets force a chain of 24-cells, with opposite cuts in intermediate cells. This is not a classification of arbitrary entirely ideal polytopes or a universal volume-divisibility theorem. |
| Compact constraints | f≥2V+5, neighboring-excess bounds, local facet cutoffs 46/46/96, zero-excess recognition, W=0 or W≥17, and degree two implies W≥20. |
| Compact minimum/second volume | Still open here. The constraints do not prove that the 120-cell minimizes volume or that its double has the second volume. |

These occupy the last two existing papers. Their
[ideal proof map](ideal_ra4_volume_and_structure_proof_map.md) and
[compact reproduction map](compact_ra4_volume_constraints_reproduce_20261019.md)
identify published inputs and actual script/formal evidence.

## Resource discipline and preserved history

Only root starts proof jobs. All scientific, TeX and Lean batches share
one exclusive lock and run serially. Admission requires 40 GiB free disk
and 6 GiB available RAM; abort thresholds stay 8/3 GiB. Scientific/TeX
jobs keep their 59-second, 1-GiB, one-thread caps; Lean uses one module,
one thread and 2,048 MiB. Owned temporary products are removed after each
batch. Below 20 GiB free disk, start cleanup before more jobs.

The former 1,000-line cumulative status page and the predecessor
result map are byte-preserved in the
[documentation archive manifest](../results/source_archive/six_paper_completed_audits_status_20261021/manifest.json).
Historical execution counts and intermediate bounds belong to those
records, not the current status. Lean remains a final verification step;
resource gates do not prevent manual mathematical review.
