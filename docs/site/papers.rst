Six papers, their techniques and complete reading copies
========================================================

The project has six paper families. The three illustrated dimension editions
below add explanations in figures to the reviewed consolidated sources. Removing
the marked additions recovers those base sources byte for byte. No theorem,
case exclusion or computation is strengthened by adding a figure.

.. list-table:: Current complete reading copies
   :header-rows: 1
   :widths: 29 11 28 32

   * - Paper and result
     - Pages
     - Technique
     - Complete source and proof map
   * - :download:`General finite volume: 57 <../../manuscripts/build/finite_volume57_illustrated.pdf>`
     - 50
     - Corner payments, spherical forests and cusp curvature.
     - :download:`Illustrated TeX <../../manuscripts/finite_volume57_illustrated.tex>`; :download:`proof map <../finite_volume57_lemma_map.md>`.
   * - :download:`General compact: 17 <../../manuscripts/build/compact_farthest_vertex_bounds_illustrated.pdf>`
     - 286
     - One farthest vector, funded component energies and actual-rank tails.
     - :download:`Illustrated TeX <../../manuscripts/compact_farthest_vertex_bounds_illustrated.tex>`; :download:`complete case map <../compact17_consolidated_proof_map_20261025.md>`.
   * - :download:`Compact arithmetic: 15 <../../manuscripts/build/compact_arithmetic_dimension_bound_illustrated.pdf>`
     - 217
     - The general compact theorem, field restrictions, four computed outputs and manual dyadic exclusion.
     - :download:`Illustrated TeX <../../manuscripts/compact_arithmetic_dimension_bound_illustrated.tex>`; :download:`dependency map <../compact_arithmetic_consolidation_proof_map_20261025.md>`.
   * - :download:`Sharp right-angled dimensions: 8 and 4 <../../manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf>`
     - 5
     - Actual projection, heights, a cubic tensor identity and metric equality.
     - :download:`TeX <../../manuscripts/right_angled_sharp_dimension_bounds_publication.tex>`; :download:`proof map <../right_angled_sharp_dimensions_publication_proof_map_20261021.md>`.
   * - :download:`Ideal RA4 volume and structure <../../manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf>`
     - 17
     - Defect count, equality spectrum, actual cutting and restricted recognition.
     - :download:`TeX <../../manuscripts/ideal_ra4_volume_and_structure_publication.tex>`; :download:`proof map <../ideal_ra4_volume_and_structure_proof_map.md>`.
   * - :download:`Compact RA4 volume constraints <../../manuscripts/build/compact_ra4_volume_constraints.pdf>`
     - 14
     - Facet capacity, neighboring excess, marked polygon geometry and graph arguments.
     - :download:`TeX <../../manuscripts/compact_ra4_volume_constraints.tex>`; :download:`proof map <../compact_ra4_volume_constraints_reproduce_20261019.md>`.

How much of a paper is essential?
---------------------------------

These are complete article lengths, not minimum necessary-page claims. The
compact 17 and arithmetic 15 articles have about 19 and 38 pages of common
exposition before their case arguments, respectively. Those front sections
are **not standalone proofs**. Actual-rank tails, exceptional endpoint/fork
cases and the four finite-output contracts remain necessary. Arithmetic 15
also depends on the independent general compact 17 theorem; its 217 pages
must not be presented as eliminating that companion dependency.

The finite-volume appendices likewise supply proofs and finite-table
interfaces used by the main argument. The five-page sharp right-angled article
is already short. The 17-page ideal paper deliberately combines the minimum,
second-volume equality and restricted chain argument in one sequence.
The compact RA4 article proves constraints rather than the requested absolute
minimum or second-volume classification. No smaller self-contained extraction
has been certified for any family by this illustration update.

Read the proof before the logs
------------------------------

* :doc:`finite_volume` follows the single corner rule through finite and ideal
  links and explains why exact dimension 58 and the higher tail are separate.
* :doc:`compact_bounds` starts with the common geometric budget, then identifies
  the classical branches and the additional arithmetic restrictions.
* :doc:`right_angled` writes out the projection and cubic code argument, including
  the strict endpoint, and the different ideal metric-equality branch.
* :doc:`volume_problems` separates the independently proved equality spectrum
  from the restricted Oct/D cutting theorem and the open compact questions.

The :doc:`algorithms` chapter describes finite domains, induction arguments and
computed propositions. :doc:`formal_logic` exposes the remaining hypotheses of
the selected Lean kernels. All six families use external mathematical inputs;
:doc:`scope` explains which ones remain assumed and where authentication is
incomplete.

Preserved earlier editions
---------------------------

Earlier PDFs remain downloadable with their original sources and claims. They
are editions of the same papers, not additional papers.

.. list-table:: Dimension-paper predecessors
   :header-rows: 1
   :widths: 27 24 24 25

   * - Family
     - Frozen edition
     - Consolidated base
     - New illustration content
   * - Finite volume
     - :download:`57, 48 pages <../../manuscripts/build/finite_volume57.pdf>`; :download:`TeX <../../manuscripts/finite_volume57.tex>`.
     - :download:`57, 50 pages <../../manuscripts/build/finite_volume57_consolidated.pdf>`; :download:`TeX <../../manuscripts/finite_volume57_consolidated.tex>`.
     - Two new diagrams; the existing figure remains.
   * - General compact
     - :download:`23, 36 pages <../../manuscripts/build/compact_farthest_vertex_bounds.pdf>`; :download:`TeX <../../manuscripts/compact_farthest_vertex_bounds.tex>`.
     - :download:`17, 285 pages <../../manuscripts/build/compact_farthest_vertex_bounds_consolidated.pdf>`; :download:`TeX <../../manuscripts/compact_farthest_vertex_bounds_consolidated.tex>`.
     - Two new diagrams sharing the original vector and distinguishing the higher tail.
   * - Compact arithmetic
     - :download:`20, 21 pages <../../manuscripts/build/compact_arithmetic_dimension_bound.pdf>`; :download:`TeX <../../manuscripts/compact_arithmetic_dimension_bound.tex>`.
     - :download:`15, 216 pages <../../manuscripts/build/compact_arithmetic_dimension_bound_consolidated.pdf>`; :download:`TeX <../../manuscripts/compact_arithmetic_dimension_bound_consolidated.tex>`.
     - Two new diagrams separating field filters from the geometric norm budget.

The three right-angled reading copies in the first table are unchanged from
the accepted editions and already contain ten figures in total. The six frozen
edition hashes are recorded in the
:download:`original delivery manifest <../../results/six_paper_final_tex_delivery_20261021.json>`.
The compact/arithmetic base sources have
:download:`separate full source acceptance <../../results/compact17_arithmetic15_final_working_source_root_acceptance_20261025.json>`
and :download:`editorial and matching-PDF acceptance <../../results/compact17_arithmetic15_editorial_completion_20261025.json>`.

A further :download:`compact RA4 publication source <../../manuscripts/compact_ra4_volume_constraints_publication.tex>`
consolidates later reviewed local constraints in the same family. Its
:download:`source map <../compact_ra4_publication_proof_and_verification_map_20261024.md>`
states its distinct scope. No replacement matching PDF has been built, so it
is not silently substituted for the 14-page reading copy above. Parked
field-specific refinements and exploratory stronger-bound work are preserved
outside the complete reading editions.

Commands and their scope
------------------------

The historical catalogue identifiers ``compact20`` and
``compact_farthest_bounds`` select the frozen 20 and 23 editions, respectively.
They do not replay every new 15/17 appendix. ``finite57``, ``right_angled_sharp``,
``ideal_ra4_volume_and_structure`` and ``compact_ra4`` also
retain their declared selection scopes. Use :doc:`reproducibility` to inspect
a plan before deciding which actual new run is needed.
