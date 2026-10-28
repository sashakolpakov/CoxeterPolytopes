Six papers, their techniques and complete reading copies
========================================================

The project has six paper families. The three illustrated dimension editions
below add explanations in figures to the reviewed consolidated sources. Removing
the marked additions recovers those base sources byte for byte. No theorem,
case exclusion or computation is strengthened by adding a figure. The two
28 October right-angled submission editions add historical context, proof
roadmaps and precise reproduction references to the existing arguments.

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
   * - :download:`Sharp right-angled dimensions: 8 and 4 <../../manuscripts/build/arxiv_sharp_dimension_bounds_8_and_4_20261028_layout_final.pdf>`
     - 7
     - Actual projection, heights, a cubic tensor identity and metric equality.
     - :download:`Submission TeX <../../manuscripts/revisions/right_angled_sharp_dimension_bounds_arxiv_20261028.tex>`; :download:`proof map <../right_angled_sharp_dimensions_publication_proof_map_20261021.md>`.
   * - :download:`Ideal RA4 volume and structure <../../manuscripts/build/arxiv_ideal_four_dimensional_volume_and_structure_20261028_historical_context.pdf>`
     - 19
     - Defect count, equality spectrum, actual cutting and restricted recognition.
     - :download:`Submission TeX <../../manuscripts/revisions/ideal_ra4_volume_and_structure_arxiv_20261028.tex>`; :download:`proof map <../ideal_ra4_volume_and_structure_proof_map.md>`.
   * - :download:`Compact RA4 volume constraints <../../manuscripts/build/compact_ra4_volume_constraints.pdf>`
     - 14
     - Facet capacity, neighboring excess, marked polygon geometry and graph arguments.
     - :download:`TeX <../../manuscripts/compact_ra4_volume_constraints.tex>`; :download:`proof map <../compact_ra4_volume_constraints_reproduce_20261019.md>`.

.. _arxiv-submission-editions:

The two arXiv submission editions
---------------------------------

Each 28 October source is a single self-contained TeX file: its figures are
inline TikZ and its bibliography is inline. The sharp-dimension paper contains
two figures; the ideal-volume paper contains five. Both name Alexander Kolpakov
and include the requested AI-use acknowledgment. The
:download:`edition guide <../arxiv_submission_editions_20261028.md>` and
:download:`source/PDF release manifest <../../results/arxiv_publication_release_manifest_20261028.json>`
record the exact editions and matching isolated, zero-warning PDF builds.

The sharp paper adds an ideal-branch diagram and a short introduction explaining
how its common projection closes the published dimension gaps. The ideal paper
expands the equality reduction into full prose and puts the second ideal volume
in its line of antecedents: Potyagailo--Vinberg's compact facet question,
Guo--Zheng's compact facet and volume minima, and Kolpakov's earlier ideal
facet and volume minima. The compact result is cited as a related preprint;
it is not an input to the ideal proof.

Their reproduction sections link to the exact companion suites
``right_angled_sharp`` and ``ideal_ra4_volume_and_structure`` at public snapshot
``04c9bab7``. The two and four script checks, respectively, use small tracked
inputs rather than the unuploaded bulk collection. Existing mathematical
outputs and conditional Lean records retain their original scope. The
submission PDF builds are new typesetting evidence, not new proof executions;
:ref:`repro-arxiv-editions` explains the distinction and the commands.

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
interfaces used by the main argument. The sharp right-angled argument has a
preserved five-page edition and a seven-page submission edition with a second
figure and reproduction section. The nineteen-page ideal submission paper
deliberately combines the known minimum, second-volume equality and restricted
chain argument in one sequence.
The compact RA4 article proves constraints rather than the requested absolute
minimum or second-volume classification. No smaller self-contained extraction
has been certified for any family by these presentation revisions.

Read the proof before the logs
------------------------------

* :doc:`finite_volume` follows the single corner rule through finite and ideal
  links and explains why exact dimension 58 and the higher tail are separate.
* :doc:`compact_bounds` starts with the common geometric budget, then identifies
  the classical branches and the additional arithmetic restrictions.
* :doc:`right_angled` writes out the projection and cubic code argument, including
  the strict endpoint, and the different ideal metric-equality branch.
* :doc:`volume_problems` separates the independently proved equality spectrum
  from the restricted Oct/D cutting theorem, the related compact minimum
  preprint and the remaining compact second-volume questions.

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

.. list-table:: Preserved right-angled submission predecessors
   :header-rows: 1
   :widths: 30 35 35

   * - Family
     - Earlier accepted PDF
     - Its unchanged source
   * - Sharp dimensions 8 and 4
     - :download:`Five-page edition <../../manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf>`
     - :download:`Earlier TeX <../../manuscripts/right_angled_sharp_dimension_bounds_publication.tex>`.
   * - Ideal volume and structure
     - :download:`Seventeen-page edition <../../manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf>`
     - :download:`Earlier TeX <../../manuscripts/ideal_ra4_volume_and_structure_publication.tex>`.

The fourteen-page compact RA4 reading copy remains the accepted edition with
its four figures. The six frozen edition hashes are recorded in the
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
a plan before deciding which actual new run is needed. In particular, a
catalogue ``--pdf`` selection for either right-angled suite still typesets
its preserved reference edition, not the later submission copy above.
