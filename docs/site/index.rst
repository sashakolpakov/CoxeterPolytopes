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
Arithmetic16 still has thirteen groups.

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
Costly and shared-capacity exits remain, so no56 theorem is claimed.
The :download:`latest glued-host continuation <../../results/review_glued_e7_a2_large_facet_release_full_scope_root_acceptance_20261023.md>`
proves an actual edge leaving the large facet, using all seventeen possible
first-wall profiles. Finite chamber volume is unresolved; this is not a new
nine-dimensional example.

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
