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
The literal norm239 and norm379 construction lattices are proved
nonreflective by original integral infinite-order chamber periods.
The norm379 fourth edge reaches a full E8 corner with an isolated long
wall. Its `complete proof and acceptance
<https://github.com/sashakolpakov/CoxeterPolytopes/blob/main/results/review_prime379_fourth_germ_original_period_root_acceptance_20261022.md>`_
bind five actual profile runs and a separate combine, all15 selectors
and294 scalar rows. The earlier RAM refusal remains preserved with zero
execution credit. These conclusions concern the two original lattices only.
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

