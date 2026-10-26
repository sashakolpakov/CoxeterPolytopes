Conditional Lean verification
=============================

The formal project checks selected implications inside the written arguments.
Its sources use Lean 4.32.2 and the standard library, pinned by
:download:`lean-toolchain <../../formal/lean-toolchain>` and
:download:`the Lake configuration <../../formal/lakefile.toml>`.
The project does not claim a formal construction of every geometric object,
a verified library of number fields, or a complete formalization of any long
manuscript. The boundary is part of each theorem's meaning.

.. figure:: _static/evidence-layers.svg
   :alt: Written proof, exact computation and conditional Lean are distinct evidence layers; publication checks bind their sources without proving their mathematical premises.
   :width: 100%

   Source bindings connect the evidence layers. A PDF build or matching hash
   does not supply a missing geometric, classification or enumeration argument.

What a conditional theorem says
-------------------------------

A typical kernel has a theorem of the form

.. math::

   H_{\mathrm{geometry}}\;\land\;H_{\mathrm{cover}}\;\land\;
   H_{\mathrm{output}}\quad\Longrightarrow\quad C.

Lean verifies the implication. To apply it to an actual polytope, the written
argument must construct the objects and establish all three hypotheses with
the exact meanings expected by the theorem. A predicate called ``ActualWall``
is not, by its name, an instantiated hyperbolic wall. A complete list supplied
as an argument remains a completeness premise until its producer and geometric
coverage are justified.

The distinction is particularly transparent in the actual audited theorem
``QuadraticDyadicTripleKernel.global_dimension_fifteen``:

.. code-block:: lean

   theorem global_dimension_fifteen (n : Nat)
       (prior_global_ceiling : n ≤ 16)
       (dimension_sixteen_excluded : n ≠ 16) : n ≤ 15 := by
     omega

This proves a useful final numerical implication. It does not prove the
prior ceiling or the exclusion of dimension sixteen. The complete nineteen-group
dispatch, actual quadratic field and geometric interpretation remain written
inputs, identified in the
:download:`dyadic kernel map <../quadratic_dyadic_triple_kernel_logic_20261023.md>`.

Compilation, axiom audit and interpretation
--------------------------------------------

The records distinguish three separate acts.

#. **Compilation** checks that the named source elaborates to a valid Lean
   object with the recorded source, imports, configuration and toolchain.
#. **An actual axiom audit** requires matching compiled-target, source and
   dependency fingerprints, then separately prints the named declarations'
   axiom dependencies. It has its own command, output and receipt.
#. **Interpretation** relates the theorem's ordinary arguments to the actual
   geometry, number field or finite computation used in the paper.

The first two have machine execution records. The third is not obtained
merely by printing an axiom list. In particular, a theorem can have no custom
axiom and still assume a geometric proposition as an ordinary argument.
An implication from an uninstantiated or inconsistent collection of hypotheses
would not establish the desired existence or dimension statement about polytopes.
The written interface review is therefore indispensable.

The preserved audited outputs use Lean's standard foundational axioms such as
``propext``, ``Quot.sound`` and, where required, ``Classical.choice``. They are
not an inventory of the literature the manuscript cites. Their absence would
not authenticate an external classification either. Successful records report
no ``sorryAx`` for the selected declarations.

Current compact and arithmetic kernels
---------------------------------------

.. list-table:: Selected current internal deductions
   :header-rows: 1
   :widths: 29 9 31 31

   * - Module or package
     - Audited declarations
     - Internal result
     - What remains outside
   * - D18UniformDottedHeadScalarKernel
     - 7
     - Scalar sealing, individual anchor affordability, far-next exclusion and near/far contradictions.
     - Actual head points, continuous estimates, original coefficients and geometric instance.
   * - CompactDimensionSeventeenComposition
     - 3
     - Same-vertex discrete exclusions and the numerical compact bound of 17.
     - Complete forest cover, every geometric consumer and higher-dimensional tail.
   * - QuadraticDyadicTripleKernel
     - 14
     - Valuation/parity, bounded integer norms and the final numerical implication yielding 15.
     - Actual maximal orders, splitting interpretation, field squares and exclusion of exact dimension 16.
   * - B16IsolatedPairScalarKernel
     - 12
     - Cubic and paid-prefix scalar contradictions.
     - Actual cap/field domains, census interpretation and original pressures.
   * - B18TerminalPressureKernel
     - 10
     - Continuous scalar terminal-pressure implications.
     - Actual terminal geometry, complete roles and legally available coefficients.
   * - A13IEndpointGraftScalarKernel
     - 11
     - Continuous graft gains and paid-square implications.
     - Actual endpoint/foot configuration and full domain coverage.
   * - A12AdjustedLossKernel
     - 10
     - Adjusted loss and finite-list square deductions.
     - The actual list and its geometric hypotheses.
   * - A18TwoDottedEndpointScalarKernel
     - 10
     - Continuous endpoint inequalities.
     - Actual endpoint cores, original-vector pressures and native domains.
   * - HeavyCentreForestScalarKernel
     - 11
     - Cone, conditional chain and incidence algebra.
     - Complete forest interpretation and geometric inequalities.

Each row has a separate execution package, rather than a new aggregate formal
proof of the paper. The direct records are:

* :download:`two compact 17 modules and their separate audits <../../results/compact17_targeted_lean_completion_20261025/manifest.json>`;
* :download:`quadratic dyadic triples <../../results/lean_quadratic_dyadic_triple_verification_20261023/manifest.json>`;
* :download:`B16 pair scalars <../../results/lean_b16_isolated_pair_scalar_verification_20261023/manifest.json>`;
* :download:`B18 terminal pressure <../../results/lean_b18_terminal_pressure_verification_20261023/manifest.json>`;
* :download:`A13 endpoint graft <../../results/lean_a13_i_endpoint_graft_verification_20261023/manifest.json>`;
* :download:`A12 adjusted loss <../../results/lean_a12_adjusted_loss_verification_20261023/manifest.json>`;
* :download:`A18 dotted endpoints <../../results/lean_a18_two_dotted_endpoint_verification_20261023/manifest.json>`;
* :download:`heavy-centre forest <../../results/lean_heavy_centre_forest_verification_20261023/manifest.json>`.

Earlier arithmetic 16 and compact 18 kernels retain their own results. Their
names do not relabel them as geometric formalizations of the later 15/17 papers.
The first D18 scalar compilation failed on a reserved binder. The reviewed
three-token rename, original failure and successful retry are all preserved in
the compact 17 package. Earlier resource deferrals are not successful runs.

Right-angled and finite-volume interfaces
-----------------------------------------

The sharp dimension paper's release package records three selected kernels
with thirty audited declarations, plus a separately recorded regular-24-cell
supplement with five. The ideal RA4 package selects eight declarations in
``ManuscriptIdeal4``. Its geometric cutting, recognition and gluing premises
remain explicit; an induction on a supplied one-step cutting theorem does not
prove that theorem or extend its facet-type hypothesis to arbitrary polytopes.
See :doc:`right_angled` and :doc:`volume_problems`.

For finite volume 57, forest induction and table arithmetic, contextual covers,
concentration and shape-curvature logic have distinct kernels. For example,
a Bellman inequality proof bounds every forest built from its allowed items;
it does not establish that those items and weights cover every actual vertex.
The mathematical bridge belongs to :doc:`finite_volume` and the
:download:`assumption register <../manuscript_assumptions.md>`.
The historical selection of nine Lean roots and 146 total jobs is not a fresh
combined runtime PASS. Current rule applications have the separate 259-unit record described
in :doc:`algorithms`.

The complete predecessor module inventory, including generated tables,
legacy dependencies, parked field-specific kernels and their stated limits,
is retained as
:download:`detailed formal-verification records <../formal_verification_records_20261025.md>`.
Those dated records are not silently updated into new execution credit.

Replaying one target safely
----------------------------

Use a named target only when changed inputs or a specific reproduction task
make a replay necessary. Choose new ledger names and run compilation before
a separate actual audit:

.. code-block:: sh

   python3 -B scripts/run_proof_batch_20261024.py lean-compile \
     CoxeterProofs.CompactDimensionSeventeenComposition \
     --ledger results/local_compact17_compile.json
   python3 -B scripts/run_proof_batch_20261024.py lean-audit \
     CoxeterProofs.CompactDimensionSeventeenComposition \
     --ledger results/local_compact17_audit.json

The public controller monitors one owned process tree and uses the shared
exclusive lock. The Lean compiler cap is 2048 MiB with one thread; the outer
watchdog retains its 59-second time limit and 1 GiB RSS limit. Admission uses
exact byte ratios against total filesystem capacity and physical RAM: at least
2% free disk and 15% available RAM. Stop an active job below 1% free disk or
6% available RAM. See :doc:`reproducibility`.
Do not substitute an uncapped whole-project ``lake build``.

Caches and cleanup are part of provenance
-----------------------------------------

A source-bound object needs its guarded compilation receipt and matching source,
imports and toolchain. Legacy objects without that link are labelled
``LEGACY_CACHED_UNCAPPED``; reuse earns no new capped compilation credit.
Objects required for a pending separate audit remain until that audit ends.
After the work item, rebuildable products can be removed while sources,
actual receipts, final PDFs and necessary legacy bytes are retained.

The cleanup archive records original paths, hashes and nanosecond timestamps.
It is an ignored local preservation archive, not a Git dependency download.
A fresh checkout can give sources newer timestamps, and the runner must refuse
stale reuse. Never change timestamps to force cache acceptance. Restored object
bytes do not prove a missing source-to-legacy-object relationship.
The :download:`cleanup record <../../results/tex_lean_repro_only_cleanup_20261025.json>`
and :doc:`sources` state what is available. A complete clean-checkout execution
of every historical formal dependency is not claimed.
