The three manuscripts and their proofs
======================================

The project has three current arXiv submission folders.
One article joins the finite-volume, compact and compact arithmetic dimension
proofs. The sharp right-angled dimensions and second ideal four-volume remain
separate articles. Compact RA4 constraints remain a supporting paper outside
these submissions. The mathematical conclusions are unchanged.

.. list-table:: Current complete reading copies
   :header-rows: 1
   :widths: 29 11 28 32

   * - Paper and result
     - Pages
     - Technique
     - Complete source and proof map
   * - :download:`Dimension bounds: 57 / 17 / 15 <../../manuscripts/build/arxiv_hyperbolic_coxeter_dimension_bounds_20261028_current_20261030.pdf>`
     - 551
     - Corner payments and cusp curvature; farthest-vertex component energies; arithmetic field restrictions.
     - :download:`Submission TeX <../../arxiv/hyperbolic_coxeter_dimension_bounds/main.tex>`; :download:`part and proof guide <../../arxiv/DIMENSION_BOUNDS.md>`.
   * - :download:`Sharp right-angled dimensions: 8 and 4 <../../manuscripts/build/arxiv_sharp_dimension_bounds_8_and_4_20261028_current_20261030.pdf>`
     - 6
     - Actual projection, heights, a cubic tensor identity and metric equality.
     - :download:`Submission TeX <../../arxiv/sharp_dimension_bounds_8_and_4/main.tex>`; :download:`proof map <../right_angled_sharp_dimensions_publication_proof_map_20261021.md>`.
   * - :download:`Second ideal RA4 volume and structure <../../manuscripts/build/arxiv_ideal_four_dimensional_volume_and_structure_20261028_current_20261030.pdf>`
     - 18
     - Defect count, equality spectrum, actual cutting and restricted recognition.
     - :download:`Submission TeX <../../arxiv/ideal_four_dimensional_volume_and_structure/main.tex>`; :download:`proof map <../ideal_ra4_volume_and_structure_proof_map.md>`.

The :download:`compact RA4 source <../../manuscripts/compact_ra4_volume_constraints.tex>`
and :download:`proof map <../compact_ra4_volume_constraints_reproduce_20261019.md>`
supply supporting constraints on facet capacity, neighboring excess and metric
recognition. They do not settle the compact second-volume problem.

.. _arxiv-submission-editions:

The three arXiv submission editions
-----------------------------------

Each source is a single self-contained TeX file: its figures are
inline TikZ and its bibliography is inline. The sharp-dimension paper contains
two figures; the ideal-volume paper contains five; the combined dimension
paper contains nine. All name Alexander Kolpakov
and include the requested AI-use acknowledgment. The
:download:`edition guide <../arxiv_submission_editions_20261028.md>` and
:download:`current source/PDF manifest <../../results/arxiv_submission_current_manifest_20261030.json>`
record the exact editions and matching isolated, zero-warning PDF builds.

The sharp paper adds a diagram for the ideal case and an introduction explaining
how its common projection closes the published dimension gaps. The ideal paper
expands the equality reduction into full prose and puts the second ideal volume
in its line of antecedents: Potyagailo--Vinberg's compact facet question,
Guo--Zheng's compact facet and volume minima, and Kolpakov's earlier ideal
facet and volume minima. The compact result is cited as a related preprint;
it is not an input to the ideal proof.

The repository guide identifies the companion suites
``right_angled_sharp`` and ``ideal_ra4_volume_and_structure``. The two and four script checks, respectively, use small tracked
inputs rather than the unuploaded bulk collection. Existing mathematical
outputs and conditional Lean records retain their original scope. The
submission PDF builds are new typesetting evidence, not new proof executions;
:ref:`repro-arxiv-editions` explains the distinction and the commands.

The combined article brings the three general/arithmetic proof families into
one document, with the compact theorem now an internal input to the arithmetic
argument. Its main exposition occupies pages 1--71. Complete case appendices
begin on page 72; the finite-volume, compact and arithmetic appendix groups
begin on pages 72, 107 and 372, respectively. The
:download:`dimension-paper guide <../../arxiv/DIMENSION_BOUNDS.md>` identifies
the original proof maps and the scope of the source-preservation review.

How much of a paper is essential?
---------------------------------

These are complete article lengths, not minimum necessary-page claims. The
combined dimension paper has 71 pages of main exposition before its case
appendices. Those front sections are **not a standalone proof**. Actual-rank
tails, exceptional endpoint/fork cases and the four arithmetic finite-output
contracts remain necessary. The independent compact theorem is included
before its use in the arithmetic argument.

The finite-volume appendices likewise supply proofs and finite-table
interfaces used by the main argument. The six-page sharp right-angled paper
contains both complete proofs and two figures. The eighteen-page ideal paper
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
