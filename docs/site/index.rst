Coxeter polytopes: proofs and verification
==========================================

This repository studies dimension bounds, reflection chambers and
right-angled volume and structure problems. The documentation explains the
arguments, their finite algorithms and the limits of the verification.
Updated 28 October 2026.

There are **six paper families**. The latest reading copies include illustrated
revisions of the three general/arithmetic dimension papers and two prepared
arXiv editions: sharp right-angled dimensions and ideal four-dimensional volume
and structure. Earlier editions remain available. These presentation revisions
preserve the mathematical conclusions; no new polytope example is claimed.

Current dimension comparison
----------------------------

The repository ceilings are internally reviewed written arguments at their
listed external and computed premises. They are not externally refereed results
or complete geometric Lean formalizations. Read :doc:`scope` alongside the table.

.. list-table:: Examples, published ceilings and repository arguments
   :header-rows: 1
   :widths: 40 14 16 30

   * - Class
     - Largest example
     - Published ceiling
     - Repository ceiling and boundary
   * - General Coxeter, finite volume
     - 21
     - 995
     - **57**; current-rule applications checked, inherited tail inputs explicit.
   * - General compact Coxeter
     - 8
     - 29
     - **17**; complete written cases, selected conditional Lean.
   * - Compact arithmetic Coxeter
     - 8
     - 29
     - **15**; computed premises and the published golden-lattice theorem.
   * - Noncompact arithmetic Coxeter
     - 21
     - 21
     - **21**, published and sharp.
   * - Right-angled, finite volume, all vertex types
     - 8
     - 12
     - **8**, internally sharp at named published inputs.
   * - Right-angled, entirely ideal
     - 4
     - 6
     - **4**, internally sharp at named published inputs.
   * - Right-angled, compact
     - 4
     - 4
     - **4**, published and sharp.

The :download:`primary-source check of Bugaenko's golden-lattice theorem <../bugaenko_golden_lattice_cutoff_authentication_20261026.md>`
closes the earlier citation-authentication warning on the arithmetic 15 row.
The theorem is still an external mathematical input.

The general arguments assume no arithmeticity. For arithmetic reflection
lattices the quadratic-form model is automatic; “first type” adds no restriction
here. Compact and noncompact arithmetic cases have different proofs. Exact
references and assumptions are in the
:download:`comparison register <../current_bounds_and_verification.md>` and
:download:`assumption register <../manuscript_assumptions.md>`.

Choose a reading route
----------------------

Start with :doc:`methods` for the common approach and where each proof branches.
:doc:`papers` maps results to complete manuscripts and gives current page counts.
The mathematical chapters develop the arguments without requiring readers to
navigate execution logs first.

.. toctree::
   :caption: Mathematics and manuscripts
   :maxdepth: 1

   papers
   methods
   finite_volume
   compact_bounds
   right_angled
   volume_problems

The verification chapters state each algorithm's domain, termination,
completeness and exact output before discussing commands and receipts. A
successful script is credited only for its declared proposition. Conditional
Lean checks the internal deductions from explicitly supplied premises.

.. toctree::
   :caption: Evidence and reproduction
   :maxdepth: 1

   verification
   algorithms
   formal_logic
   reproducibility
   scope

.. toctree::
   :caption: Examples and sources
   :maxdepth: 1

   examples
   sources

What this release contains
--------------------------

The :ref:`two 28 October submission editions <arxiv-submission-editions>`
have complete, self-contained TeX sources and matching reading PDFs:

* :download:`Sharp dimensions 8 and 4, seven pages <../../manuscripts/build/arxiv_sharp_dimension_bounds_8_and_4_20261028_layout_final.pdf>`:
  two inline figures, a proof roadmap and exact reproduction references.
* :download:`Ideal volume and structure, nineteen pages <../../manuscripts/build/arxiv_ideal_four_dimensional_volume_and_structure_20261028_historical_context.pdf>`:
  five inline figures, the second-volume theorem, restricted chains and the
  history of the compact and ideal extremal questions.

Both sources include the author and AI-use acknowledgment. Their isolated PDF
builds have zero warnings. The
:download:`submission-edition guide <../arxiv_submission_editions_20261028.md>`
and :download:`release manifest <../../results/arxiv_publication_release_manifest_20261028.json>`
identify the exact source/PDF pairs and retain the existing mathematical and
conditional Lean verification scopes. Updating a reading copy does not
constitute a new execution of its supporting checks.

The :download:`26 October release note <../documentation_illustration_release_20261026.md>`
records the source-preservation and illustration audit boundaries.

Five original vector diagrams illustrate the proof dependencies. Six additional
TikZ figures appear in the three illustrated general/arithmetic dimension
manuscripts. Their marked illustration blocks can be removed to recover each
reviewed base source byte for byte. The current right-angled reading copies
contain two, five and four inline figures, respectively. Compilation and
rendered-page reviews are recorded separately from mathematical source review.

Finite57 has a successful **259-unit current-rule conjunction**, including all
64 selected application shards. This does not mean a fresh combined replay of
the historical 146-job catalogue or regeneration of every ancestor enumeration.
Successful unchanged Lean modules were not recompiled for this documentation
revision. The earlier separate compilation and actual axiom audits retain their
conditional scopes.

The `GitHub source repository <https://github.com/sashakolpakov/CoxeterPolytopes>`_
contains the manuscripts, code and selected evidence. The manifest describes
936 scientific files in 33 additional archives, but **no bulk asset has been
uploaded**. Third-party originals are also excluded. A source checkout therefore
is not a promise of complete transitive reproduction; see :doc:`sources`.
