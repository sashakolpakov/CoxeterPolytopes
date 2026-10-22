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

The latest supporting written proof excludes the last whole D19 family
through all six long-head classes, completing the dimension-nineteen frontier
and proving global compact dimension at most18. Its `acceptance and full review bindings
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/compact_dimension18_global_written_acceptance_20261022.json>`_
record the complete manual argument. Four arithmetic seventeen groups remain;
the D16+actual-I branch now requires a quadratic actual field after complete
quartic and cubic exclusions. The entire quadratic forest remains open.
The literal norm239 construction lattice is now proved nonreflective
by an original integral infinite-order chamber period. A distinct norm379
host now has reviewed E7+I2(5) and subsequent E7+2A1 full finite germs,
each with five actual profile runs and a separate combine. The earlier
RAM refusal is preserved with zero execution credit.
No new polytope example is claimed.
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

