Coxeter polytopes: results and verification
===========================================

This project studies dimension bounds, explicit reflection chambers and
right-angled volume/structure problems. There are **six paper families**,
with matching TeX/PDF editions and separately scoped mathematical reviews,
exact computations and Lean kernels. Updated23 October2026.

Current dimension bounds
------------------------

The new ceilings are internally reviewed written proofs. They are distinct
from published theorems, external refereeing and full geometric formalization.

.. list-table:: Dimension comparison
   :header-rows: 1
   :widths: 43 16 18 23

   * - Class
     - Largest example
     - Published ceiling
     - Repository ceiling
   * - General Coxeter, finite volume
     - 21
     - 995
     - **57**
   * - General compact Coxeter
     - 8
     - 29
     - **18**
   * - Compact arithmetic Coxeter
     - 8
     - 29
     - **16**
   * - Noncompact arithmetic Coxeter
     - 21
     - 21
     - **21**, published and sharp
   * - Right-angled, finite volume, all vertices
     - 8
     - 12
     - **8**, internally sharp
   * - Right-angled, entirely ideal
     - 4
     - 6
     - **4**, internally sharp
   * - Right-angled, compact
     - 4
     - 4
     - **4**, published and sharp

The compact18 and arithmetic16 arguments are supporting work in the existing
compact paper families; their frozen PDFs still state23 and20. Arithmetic
reflection lattices automatically admit the quadratic-form model: “first type”
is not an additional restriction. See the
:download:`current guide <../current_bounds_and_verification.md>` and
:download:`published-input register <../manuscript_assumptions.md>`.

At exact compact dimension18, only **wholeB18 or wholeD18** remains in the
accepted necessary frontier. The complete
:download:`B17+actual-I exclusion <../../results/review_compact18_whole_b17_actual_i_and_four_group_frontier_root_acceptance_20261023.md>`
and :download:`D17+actual-I exclusion <../../results/review_compact18_whole_d17_actual_i_and_three_group_frontier_root_acceptance_20261023.md>`
remove two groups, and the new
:download:`whole-A18 exclusion <../../results/review_compact18_whole_a18_and_two_classical_group_frontier_root_acceptance_20261023.md>`
removes the third. It uses complete actual endpoint classification, the native
angle gap and a positive H4+A1 projection in the last case. This geometric
result is written and internally reviewed, not Lean-verified. The remaining
groups are necessary configurations, not examples; the ceiling stays18.
Arithmetic16 now has :download:`nine original groups <../../results/review_arithmetic16_a13_actual_i_ordinary_dihedral_group_nine_root_acceptance_20261023.md>`.
Its entire A13 row is excluded by universal energy, a direct cubic7-field
integral-basis argument and a simultaneous actual-I endpoint/weak-foot square.
No enumeration was needed for the cubic step. The ceiling stays16.
The :download:`new arithmetic A14 theorem <../../results/review_arithmetic16_a14_point597_four_fifths_root_acceptance_20261023.md>`
gives energy greater than4/5 and removes its two-ordinary-A1 alternative.
OrdinaryT2 and actual typed-seven alternatives remain; the nine IDs are unchanged.

The six papers
--------------

.. list-table:: Current frozen editions
   :header-rows: 1
   :widths: 53 30 17

   * - Paper
     - Result in PDF
     - Read
   * - General finite volume
     - Dimension at most57
     - :download:`48 pages <../../manuscripts/build/finite_volume57.pdf>`
   * - Sharp right-angled dimensions
     - Finite volume8; entirely ideal4
     - :download:`5 pages <../../manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf>`
   * - General compact
     - Dimension at most23
     - :download:`36 pages <../../manuscripts/build/compact_farthest_vertex_bounds.pdf>`
   * - Compact arithmetic
     - Dimension at most20
     - :download:`21 pages <../../manuscripts/build/compact_arithmetic_dimension_bound.pdf>`
   * - Ideal right-angled four-volume and structure
     - Minimum, second volume, restricted chains
     - :download:`17 pages <../../manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf>`
   * - Compact right-angled four-volume constraints
     - Capacity, excess and recognition
     - :download:`14 pages <../../manuscripts/build/compact_ra4_volume_constraints.pdf>`

All six are **written, compiled and internally reviewed**. The57 paper has
12 pages of main exposition and36 pages of essential appendices. The
:download:`reader index <../../manuscripts/README.md>` maps each result to its
technique, source, proof map and verification ID. Earlier notes and refinements
are not additional proposed papers.

Recent verified work
--------------------

The :download:`new B18 terminal theorem <../../results/review_compact18_b18_uniform_terminal_lengths_complete_seventeen_inventory_root_acceptance_20261023.md>`
gives unconditional SAME-vector sums S17<3.760 and S18<3.830, actual edge
length bounds and a complete seventeenth-white inventory. The entire
:download:`near3/own6 role <../../results/review_compact18_b18_near_three_own_six_complete_short_exclusion_root_acceptance_20261023.md>`
and :download:`all finite-own far-helper roles <../../results/review_compact18_b18_all_finite_far_helper17_raw_prefix_root_acceptance_20261023.md>`
are excluded. The :download:`complete dotted far-helper theorem <../../results/review_compact18_b18_dotted_far_helper17_complete_root_acceptance_20261023.md>`
now leaves ONLY the short-wall contact family for the seventeenth opposite.
Its actual neighbouring point is never assumed farthest; a freshly weakened
comparison closes every zero/one/two weak-end pattern. A supporting anchor
price correction is explicit and byte-preserved; all norm margins remain.
The :download:`complete D18 centre inventory <../../results/review_compact18_d18_complete_centre16_inventory_two_holdouts_root_acceptance_20261023.md>`
and :download:`strong15/finite-fork reduction <../../results/review_compact18_d18_nonpure_strong15_finite_fork3_head_root_acceptance_20261023.md>`
leave two stated weak-centre profiles. Neither wholeB18 nor wholeD18 is
excluded, and the global ceiling remains18.

The :download:`fresh arithmetic B16 head theorem <../../results/review_arithmetic16_whole_b16_head_near_short_three_scale_root_acceptance_20261023.md>`
and :download:`complete short-family reduction <../../results/review_arithmetic16_whole_b16_complete_short_white_isolated_dotted_field_root_acceptance_20261023.md>`
now leave both actual opposites isolated, the short one dotted. EVERY
nonisolated short role is excluded over all actual fields. The remaining
pair has an actual quadratic/cubic field, strict conjugate strips and an
additional half-product integrality condition. It is not excluded:
original group15 and the arithmetic ceiling16 remain.

The :download:`zero-exchanged row7/near-row8 financing theorem <../../results/review_finite57_zero_exchanged_all_row7_near_row8_financing_root_acceptance_20261023.md>`
retains every older payment. Full actual restorer products limit the A4
finite receiver to one qualifying arrival across both possible ports.
Far-row8 and positive exchanged products remain open; no56 theorem follows.

The :download:`B18 scalar kernel <../../results/lean_b18_terminal_pressure_verification_20261023/manifest.json>`
now has actual bounded compilation in1.71seconds and a SEPARATE ten-declaration
audit in1.16seconds. The audit reused its object; no successful old module
ran. All ten prints use standard foundational axioms and no ``sorryAx``.
Actual geometry, inventories, original-vector pressure and terminal sums
remain external. This is scalar verification for near3/own6, not wholeB18
or later geometric exclusions. See its :download:`logic map <../b18_terminal_pressure_kernel_logic_20261023.md>`.
Two unused-variable warnings are preserved; temporary proof work is absent.

The earlier :download:`A13 eleven-declaration <../../results/lean_a13_i_endpoint_graft_verification_20261023/manifest.json>`,
:download:`A12 ten-declaration <../../results/lean_a12_adjusted_loss_verification_20261023/manifest.json>`,
:download:`A18 ten-declaration <../../results/lean_a18_two_dotted_endpoint_verification_20261023/manifest.json>`
and :download:`heavy-centre eleven-declaration <../../results/lean_heavy_centre_forest_verification_20261023/manifest.json>`
audits retain their exact conditional scopes. Preparation labels are
historical; actual receipts supply execution status. Earlier failures and
resource deferrals earn zero PASS credit. No end-to-end geometric
formalization is claimed.

Both the :download:`glued H9 host <../../results/review_glued_e7_a2_two_unit_releases_original_period_and_literal_nonreflectivity_root_acceptance_20261023.md>`
and the :download:`squared-unit host <../../results/review_h9_squared_unit_e7_return_original_period_literal_nonreflectivity_root_acceptance_20261023.md>`
are now proved nonreflective by original chamber symmetries of infinite
order. The latter product has trace16phi^6=80+128phi>10; its earlier
order-two symmetry alone did not establish infinite volume. These close
two literal hosts, not their whole commensurability classes.

A different :download:`prime41 host and first actual edge <../../results/review_h9_cubed_unit_prime41_original_host_first_edge_root_acceptance_20261023.md>`
and its :download:`second actual unit edge <../../results/review_h9_prime41_second_actual_unit_edge_original_germ_obstruction_root_acceptance_20261023.md>`
are proved with complete original-lattice firstness arguments. The
:download:`third edge and original order-two symmetry <../../results/review_h9_prime41_third_actual_edge_original_involution_full_hexagon_root_acceptance_20261023.md>`
close the WHOLE compact E7 hexagon, of area4pi/3. Ambient volume remains
open; the symmetry is not an infinite period. There is no new example.

Read and reproduce
------------------

.. toctree::
   :maxdepth: 1

   papers
   methods
   verification
   reproducibility
   scope

The `GitHub repository <https://github.com/sashakolpakov/CoxeterPolytopes>`_
contains the source arguments and evidence. Before interpreting a PASS, read
its declared scope: source review, arithmetic replay, conditional kernel logic
and full geometric formalization are different checks. The optional936-file
bulk payload has not been uploaded and is excluded from Git.
