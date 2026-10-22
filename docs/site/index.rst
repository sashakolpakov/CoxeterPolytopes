Coxeter polytopes: papers and verification
==========================================

This project studies dimension bounds, volumes, and structures of hyperbolic
Coxeter polytopes. Start with the :doc:`six papers <papers>`, especially the
general finite-volume bound of 57 and the sharp right-angled bounds of 8 and 4.
Each paper has a matching compiled PDF and a proof map identifying its inputs.

The papers are complete at their declared published-result and finite-data
boundaries. Internal mathematical reviews, exact script checks, PDF builds,
and conditional Lean proofs provide different kinds of evidence. None is a
claim of external refereeing or an end-to-end formalization of the geometry.
The :doc:`verification status <verification>` keeps those distinctions explicit.

.. list-table:: Dimension results and their proof homes
   :header-rows: 1
   :widths: 42 26 32

   * - Class
     - Bound in the compiled paper
     - Later internally reviewed written result
   * - General finite volume
     - 57
     - 57
   * - General compact
     - 23
     - 18; supporting proof, not in the 23 PDF
   * - Compact arithmetic Coxeter
     - 20
     - 16 over every actual admissible field; not in the 20 PDF
   * - Right-angled finite volume
     - 8, sharp
     - 8
   * - Entirely ideal right-angled
     - 4, sharp
     - 4

The finite-volume search also has a new
:download:`isolated parallel short-fork theorem <../../results/review_finite57_isolated_parallel_short_fork_financing_root_acceptance_20261022.md>`:
seven actual receiving cusp families finance the full selected transfers,
with minimum remaining curvature3,485. Global coverage for56 remains open.

The other two papers concern ideal and compact right-angled polytopes in
dimension four. The compact volume minimum and second-volume problems remain
open. No new higher-dimensional example is claimed here.

The latest supporting written proof excludes the last whole D19 family
through all six long-head classes, completing the dimension-nineteen frontier
and proving global compact dimension at most18. Its `acceptance and full review bindings
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/compact_dimension18_global_written_acceptance_20261022.json>`_
record the complete manual argument. The arithmetic dimension-seventeen
frontier is now also closed, including whole A17. The entire
D16+actual-I forest is excluded by
the actual-degree reduction and complete quadratic, cubic and quartic
contradictions. Its `full review and root acceptance
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/d16_actual_i_exclusion_and_finite57_hexagon_root_acceptance_20261022.md>`_
also record complete selected-domain financing for finite57 all-heavy
hexagons at actual nonadjacent cusps. Neither result lowers a global
ceiling by itself. The later `whole-D17 exclusion
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/review_arithmetic17_whole_d17_complete_exclusion_root_acceptance_20261022.md>`_
uses complete field reductions and a36-field/120-pair exact census. Its
two survivors both violate the same legally funded head/fork norm budget.
All four actual bounded units passed. A scalar Lean kernel compiled and
passed a separate four-declaration audit; geometry, complete enumeration
and paid-reference applicability remain explicit written inputs. The
arithmetic ceiling at that checkpoint stayed17. The later `complete B17 exclusion
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/review_arithmetic17_whole_b17_complete_exclusion_root_acceptance_20261022.md>`_
has full internal review, including all finite-head branches and the
derived B17 applicability of the unchanged census. Both surviving pairs
contradict the same funded norm inequality. The `A17 endpoint theorem
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/review_arithmetic17_whole_a17_endpoint_reduction_root_acceptance_20261022.md>`_
isolates both ends, excludes the golden field and bounds both edge cosh
lengths below2.1. The later geometric reduction forces one finite
endpoint and derives four possible actual fields. The `complete A17
quadratic exclusion
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/review_arithmetic17_whole_a17_quadratic_exclusion_root_acceptance_20261022.md>`_
has full source/runtime/composition review: three batches and a separate
complete replay checked25 fields and159 ordered pairs, with zero survivors.
This left the two actual cubic fields K7 and K9. Their `complete actual census
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/a17_finite_cubic_endpoint_runtime_manifest_20261022.json>`_
now has two field units and a separate all-record replay: all2,432 coefficient
triples were checked, with zero survivors. The `complete global16 proof
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/arithmetic17_whole_a17_complete_endpoint_census_and_global16_composition_20261022.md>`_
and its `full review and root acceptance
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/review_arithmetic17_whole_a17_global16_complete_composition_root_acceptance_20261022.md>`_
close all fourteen original forest groups. The prior global17 theorem
supplies the higher-dimensional tail, giving the global compact arithmetic
ceiling16. Its conditional final composition now has a
:download:`bounded Lean compile and separate three-declaration audit PASS <../../results/lean_arithmetic16_composition_verification_20261022/manifest.json>`,
with27 explicit external inputs and no geometric instance. The printed
arithmetic20 paper and all six editions are unchanged.

At exact compact dimension18, the
:download:`complete A16 closure <../../results/review_compact18_whole_a16_row_root_acceptance_20261022.md>`
closes both actual typed-seven alternatives without arithmeticity. The new
:download:`universal A12 energy theorem <../../results/review_compact_a12_point65_and_general18_eight_groups_root_acceptance_20261022.md>`
proves energy above13/20 with arbitrary finite spectators and closes its
whole row too. The later
:download:`complete D16 typed-seven exclusion <../../results/review_compact18_whole_d16_typed_seven_root_acceptance_and_seven_group_frontier_20261022.md>`
closes both actual N7/O7 forests without arithmeticity, leaving **seven
necessary groups** after the earlier A15/A16 exclusions. This new geometric
proof has written review credit, not Lean credit.
The :download:`ordinary dihedral floor above3/19 <../../results/review_compact_ordinary_dihedral_three_nineteenths_root_acceptance_20261022.md>`
closes arithmetic dimension16's ordinary D14 spectator. The
:download:`complete D14 acceptance <../../results/review_arithmetic16_whole_d14_root_acceptance_and_thirteen_group_frontier_20261022.md>`
also closes both actual typed-seven defects. No global ceiling or example changes.

At exact dimension16, a fresh complete nineteen-group list and the
`accepted refinement
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/review_arithmetic16_whole_a16_and_universal_a3_root_acceptance_20261022.md>`_
now leave **thirteen necessary groups**. The closures include the selected
A14+actual-I+ordinary-A1 subcase, D13 with rank-three spectators, ordinary
rank-five A11, B14 with either typed-seven defect, whole A16, A11+O7+A3 and D14.
The common advances are arithmetic A11 energy above13/20, nonarithmetic A3
energy above13/50, and one-reference actual endpoint/finite-anchor arguments.
The :doc:`verification guide <verification>` preserves the exact boundaries.
This does not prove dimension15 or realize any remaining forest. A separate
:download:`conditional frontier kernel <../arithmetic_sixteen_frontier_refinement_logic_20261022.md>`
now has a
:download:`bounded compile and separate six-declaration axiom audit PASS <../../results/lean_arithmetic16_frontier_thirteen_group_verification_20261022/manifest.json>`
under nine explicit written interfaces, with no ``sorryAx``. Failed/refused
predecessors remain preserved. Global ceilings remain57/18/16.

The distinct
:download:`glued golden E7/A2 construction <../../results/review_glued_e7_a2_second_ambient_edge_and_quadrilateral_root_acceptance_20261022.md>`
has accepted actual edges, full finite germs and complete compact triangle
and quadrilateral faces. Finite nine-dimensional chamber volume remains
unresolved; no new polytope example is claimed.

The literal norm239 and norm379 construction lattices are proved
nonreflective by original integral infinite-order chamber periods.
The norm379 fourth edge reaches a full E8 corner with an isolated long
wall. Its `complete proof and acceptance
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/review_prime379_fourth_germ_original_period_root_acceptance_20261022.md>`_
bind five actual profile runs and a separate combine, all15 selectors
and294 scalar rows. The earlier RAM refusal remains preserved with zero
execution credit. These conclusions concern the two original lattices only.
No new polytope example is claimed.
The distinct composite551 literal host is now also `nonreflective
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/review_composite551_binary_release_and_literal_nonreflectivity_root_acceptance_20261022.md>`_.
Seven actual transverse profile runs and a separate combine give an
E7+two-A1 corner. The subsequent complete manual next edge and original
integral chamber period close that module; this is not a new worker run,
example, or exclusion of its whole commensurability class.
The compact18 conditional composition now has bounded compilation and a
separate five-declaration actual axiom audit under23 named written premises.
The six compiled paper editions are unchanged.

.. toctree::
   :maxdepth: 2
   :caption: Read and reproduce

   papers
   verification
   reproducibility
   methods
   scope

The source repository is `sashakolpakov/CoxeterPolytopes
<https://github.com/sashakolpakov/CoxeterPolytopes>`_. Repository publication,
Pages deployment, and release-data availability are separate from the local
proof records described here; this site does not certify a completed upload.

