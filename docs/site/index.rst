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

The :download:`heavy-centre forest kernel <../../results/lean_heavy_centre_forest_verification_20261023/manifest.json>`
has bounded compilation and a separate actual audit of eleven declarations.
It derives ordered-ring cone bounds, finite-chain no-return from two local
geometric bridges and incidence algebra. Actual geometric and graph instances
remain external. The initial RAM refusal is preserved; no older module ran.

The :download:`A18 continuous scalar kernel receipt <../../results/lean_a18_two_dotted_endpoint_verification_20261023/manifest.json>`
records bounded compilation and a SEPARATE actual audit of ten declarations.
It proves continuous ordered-ring inequalities and paid-square logic for the
two-dotted endpoint branch. Actual geometry and legal-price premises remain
explicit. All printed dependencies are standard Lean foundational axioms;
there is no ``sorryAx``. This is not a whole-A18 formalization or a formal
verification of the later two-family frontier.

The finite-volume57 lane now has a
:download:`complete finite-centre classification, safe union and heavy escape theorem <../../results/review_finite57_finite_centre_safe_union_and_heavy_escape_root_acceptance_20261023.md>`.
The new :download:`metric forest theorem <../../results/review_finite57_heavy_centre_forest_and_linear_exit_expansion_root_acceptance_20261023.md>`
forces exactly H+2 actual exits from each heavy component of H vertices.
The :download:`mixed row1/row2 receiver theorem <../../results/review_finite57_row_one_two_variable_cone_component_financing_root_acceptance_20261023.md>`
now finances every full enlarged component whose nonheavy centre neighbours
have either receiving type, retaining all prior credits. The subsequent
:download:`zero-exchanged costly-exit financing theorem <../../results/review_finite57_zero_exchanged_bad_centre_affine_financing_root_acceptance_20261023.md>`
pays four more full-table rows through complete actual cusp-parent recovery,
retaining every shared debt. Zero-exchanged rows7/8 and positive exchanged
products remain; no56 theorem is claimed.
The :download:`latest compact reductions <../../results/review_compact18_b18_strong17_and_d18_point58_root_acceptance_20261023.md>`
bound both D18 fork cosh lengths by12, or10 under a strong5 hypothesis.
A selected :download:`strong17 B18 pure-helper/isolated-short branch <../../results/review_b18_strong_helper_actual_a17_neighbor_root_acceptance_20261023.md>`
is now excluded, including its formerly open orthogonal locus. A selected
:download:`weak pure dotted D18 centre16 branch <../../results/review_compact18_d18_weak_pure_dotted_centre_branch_complete_root_acceptance_20261023.md>`
is also excluded by complete actual fork cases.
The new :download:`unconditional isolated17 theorem <../../results/review_compact18_b18_unconditional_isolated_seventeen_root_acceptance_20261023.md>`
excludes that entire role without strong17 or a short-opposite premise. The
:download:`finite-own pure centre16 theorem <../../results/review_compact18_d18_finite_own_weak_pure_centre_root_acceptance_20261023.md>`
also excludes the finite-own version of the weak pure dotted branch.
Other centre16 roles remain; neither whole family is excluded.

The :download:`A13 scalar source/map <../a13_i_endpoint_graft_scalar_kernel_logic_20261023.md>`
have a :download:`reviewed two-line repair <../../results/review_a13_endpoint_graft_two_identity_rewrites_root_preclear_20261023.md>`
after an :download:`actual failed compile <../../results/lean_a13_i_endpoint_graft_verification_20261023/attempt_2_failed_compile.json>`.
A later admitted retry has an :download:`actual compile PASS in1.33seconds <../../results/lean_a13_i_endpoint_graft_verification_20261023/pre_audit_manifest.json>`.
The separate audit was initially refused before Lean at5.3GiB. A later
:download:`actual eleven-declaration audit PASS <../../results/lean_a13_i_endpoint_graft_verification_20261023/manifest.json>`
took1.28seconds at6.34GiB RAM. It reran no successful compilation and
uses no ``sorryAx``. Earlier failed prints earn zero credit; geometric
premises stay external.

A different :download:`squared-unit host <../../results/review_h9_squared_unit_glued_host_first_two_actual_edges_root_acceptance_20261023.md>`
has a :download:`whole compact E7 quadrilateral face <../../results/review_h9_squared_unit_third_edge_full_e7_quadrilateral_root_acceptance_20261023.md>`
and an :download:`actual ambient D6 edge leaving it <../../results/review_h9_squared_unit_first_d6_ambient_edge_root_acceptance_20261023.md>`.
Two :download:`further actual edges reach a full D6+B3 corner <../../results/review_h9_squared_unit_d6_b2_to_b3_actual_edges_root_acceptance_20261023.md>`.
The subsequent :download:`long D6-leaf edge and original chamber involution <../../results/review_h9_squared_unit_d6_long_leaf_edge_original_involution_root_acceptance_20261023.md>`
are proved. The symmetry has order exactly two and trace6. Ambient finite
volume and compactness remain open; no new example is claimed.

The :download:`glued nine-dimensional host is proved nonreflective <../../results/review_glued_e7_a2_two_unit_releases_original_period_and_literal_nonreflectivity_root_acceptance_20261023.md>`.
Two actual edge proofs and a universal original-module parity identity give
an integral chamber stabilizer of trace5+4phi>10. Its infinite order proves
infinite chamber volume. This closes one literal host and yields no example
or commensurability-class exclusion.

The new :download:`A12 scalar kernel <../a12_adjusted_loss_kernel_logic_20261023.md>`
has :download:`bounded compilation and a separate ten-declaration audit PASS <../../results/lean_a12_adjusted_loss_verification_20261023/manifest.json>`.
Compilation took2.37seconds and the separate audit1.51seconds; no ``sorryAx``
occurs. Earlier refused/failed attempts are preserved. Actual geometry and
field inventories remain external; the later universal A12/A13 theorems are
outside this scalar module. The execution receipt supplies the later status
of the immutable preparation source and logic map.

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
