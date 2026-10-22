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
     - 19; supporting proof, not in the 23 PDF
   * - Compact arithmetic Coxeter
     - 20
     - 17 over every actual admissible field; not in the 20 PDF
   * - Right-angled finite volume
     - 8, sharp
     - 8
   * - Entirely ideal right-angled
     - 4, sharp
     - 4

The other two papers concern ideal and compact right-angled polytopes in
dimension four. The compact volume minimum and second-volume problems remain
open. No new higher-dimensional example is claimed here.

The latest supporting written batch excludes whole A19 through all three
central-end cases. Only whole D19 remains in general dimension nineteen;
four arithmetic seventeen groups remain. The `acceptance and full review bindings
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/whole_a19_exclusion_and_construction_progress_acceptance_20261022.json>`_
keep these necessary forests separate from global dimension bounds.
The literal norm239 construction lattice is now proved nonreflective
by an original integral infinite-order chamber period. A distinct norm379
host has a reviewed first binary edge and a complete prospective transverse
domain; no transverse worker or new example is claimed.
A new six-declaration D15 conditional Lean kernel passed compilation and a
separate actual axiom audit. The six compiled paper editions are unchanged.

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

