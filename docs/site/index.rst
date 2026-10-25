Coxeter polytopes: results and verification
===========================================

This project studies dimension bounds, explicit reflection chambers and
right-angled volume and structure problems. There are **six paper families**,
with frozen TeX/PDF editions, supporting proofs and separately scoped reviews,
exact computations and Lean kernels. Updated25 October 2026.

The current task is :download:`auditing and consolidating the same six papers <../audit_and_consolidation_status_20261024.md>`.
Fresh audits, same-family revisions and bounded execution have separate
statuses. New scoped check receipts are listed separately from unexecuted
plans and resource deferrals. The six frozen editions remain unchanged.

Current dimension bounds
------------------------

The repository improvements are internally reviewed written proofs, distinct
from publication, external refereeing and full geometric formalization.

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
     - **17**
   * - Compact arithmetic Coxeter
     - 8
     - 29
     - **15**
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

Arithmetic reflection lattices automatically admit the quadratic-form model;
“first type” adds no restriction here. General results assume no arithmeticity.
The :download:`current comparison guide <../current_bounds_and_verification.md>`
and :download:`assumption register <../manuscript_assumptions.md>` give exact
premises and references. **No new polytope example has been established.**

Current audit and publication work
----------------------------------

Finite57 current-rule verification, both right-angled packages and the
compact17/arithmetic15 article revisions have current completion records.
Selected Lean, matching PDFs, strict documentation, actual Git-index inspection
and source delivery retain separate evidence. Further dimension and construction
searches are stopped.

The arithmetic working source now contains its reviewed local proofs and
complete forest dispatch. Its
:download:`working proof map <../compact_arithmetic_consolidation_proof_map_20261025.md>`
identifies the general companion, literal literature assumption and four
finite-output premises. Its matching216-page revision has separate source
and sampled visual acceptance, with no overfull boxes or undefined references.
Eight underfull spacing warnings remain; the frozen PDF still states20.

The compact working source now has complete mathematical and editorial
source acceptance: every exact 18 row, the actual-rank19–29 tail, all three
whole rank19 consumers and the independent final conjunction. Its
:download:`working proof map <../compact17_consolidated_proof_map_20261025.md>`
identifies the shared method and the case-specific arguments. Its matching
285-page revision has separate source and sampled visual acceptance, with no
overfull boxes or undefined references. Nine underfull spacing warnings remain.
The frozen PDF still states23. The common exposition occupies about19 pages,
followed by complete case arguments; no arithmeticity is used.

Finite57 now has :download:`all 64 actual selected application PASS receipts <../../results/finite57_longer_tail_all64_actual_application_completion_20261025.json>`
and a :download:`successful 259-unit current-rule reader <../../results/finite57_final_current_rule_conjunction_execution_20261025.json>`.
They cover28,028 maps and 83 prefixes under one frozen source/premise snapshot.
The :download:`reproduction guide <../finite_volume57_reproduce.md>` records
what was checked and which geometric and finite-data premises remain explicit.
This does not regenerate native enumeration, provide a combined historical
146-job replay or formalize the complete geometry. The current 50-page working
article includes the reviewed terminal-tail lemma; the frozen48-page edition
is unchanged. Bounded59–163 tail inputs retain their inherited numerical scope.

The :download:`right-angled release records <../right_angled_release_packages_20261025.md>`
have :download:`current root acceptance <../../results/right_angled_current_release_root_acceptance_20261025.json>`
for source, computation and actual conditional-Lean receipts for both papers.
They distinguish predecessor-edition checks from current source reviews. The two new compact17 Lean modules now have
:download:`actual compilation and separate axiom audits for ten declarations <../../results/compact17_targeted_lean_completion_20261025/manifest.json>`.
Their geometric and complete-data parameters remain external; older successful
modules retain their scopes.

The :download:`audit guide <../audit_and_consolidation_status_20261024.md>`
tracks each remaining gate. Earlier research toward stronger bounds and
constructions is preserved in the
:download:`detailed comparison guide <../current_bounds_and_verification.md>`;
those partial results do not alter the table above.

The six papers
--------------

All six editions are written, compiled and internally reviewed. Supporting
improvements remain in these families; earlier notes and atlases are not
additional proposed papers. The
:download:`paper index <../../manuscripts/README.md>` and
:download:`result correspondence <../result_manuscript_correspondence.md>`
identify the common techniques and where their cases separate.

.. list-table:: Current frozen editions
   :header-rows: 1
   :widths: 53 30 17

   * - Paper
     - Result in PDF
     - PDF
   * - General finite volume
     - At most57
     - :download:`48 pages <../../manuscripts/build/finite_volume57.pdf>`
   * - Sharp right-angled dimensions
     - Finite volume8; entirely ideal4
     - :download:`5 pages <../../manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf>`
   * - General compact
     - At most23; supporting proof17
     - :download:`36 pages <../../manuscripts/build/compact_farthest_vertex_bounds.pdf>`
   * - Compact arithmetic
     - At most20; supporting proof15
     - :download:`21 pages <../../manuscripts/build/compact_arithmetic_dimension_bound.pdf>`
   * - Ideal right-angled four-volume and structure
     - Minimum; unique second-volume double; restricted chains
     - :download:`17 pages <../../manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf>`
   * - Compact right-angled four-volume constraints
     - Capacity, excess and recognition
     - :download:`14 pages <../../manuscripts/build/compact_ra4_volume_constraints.pdf>`

The ideal four-dimensional paper does not classify all ideal polytopes as
24-cell chains. The compact four-dimensional paper establishes neither an
absolute120-cell minimum nor a unique second-volume double.

Read and reproduce
------------------

Finite57 has 25 actual bounded replacement PASS units and catalogue integration.
Its full136-script/nine-Lean-root/one-PDF selection has no combined146-job
runtime PASS or complete transitive regeneration;111 other command matches
retain inherited execution scope. Core11 is unchanged. See the
:download:`reproduction guide <../finite_volume57_reproduce.md>`.

Selected Lean kernels have actual compilation and SEPARATE actual axiom audits.
The new dyadic and B16 scalar modules have14 and 12 audited declarations;
their field, census, geometry and real-instance premises remain external.
The :download:`formal guide <../../formal/README.md>` preserves logic maps,
exact receipts, failures and resource deferrals. No selected conditional kernel
is called a full geometric formalization.

One root-owned scientific/Lean/TeX process runs at a time under the shared lock.
The user's25October policy admits at 2% free disk and 15% available RAM,
measured against total filesystem capacity and physical RAM. Stop active
jobs below 1% disk or6% RAM; start cleanup below 2% disk. Exact ratios and
policy identities are recorded. Script/TeX batches use59 seconds/1 GiB/one
thread; Lean2048 MiB/one thread followed by a separate audit. Owned unused
products are removed after each batch.
Publication I/O and strict Sphinx use their separate monitored watchdog and
must not launch mathematical checks. Generated sites, compiler objects and
caches remain untracked. The optional936-file bulk payload is **unuploaded**.

.. toctree::
   :maxdepth: 1

   papers
   methods
   verification
   reproducibility
   scope

The `GitHub repository <https://github.com/sashakolpakov/CoxeterPolytopes>`_
contains source arguments and evidence. Read each check's declared scope before
interpreting PASS: written review, arithmetic replay and conditional Lean
logic supply different evidence.
