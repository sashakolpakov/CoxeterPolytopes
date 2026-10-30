.. _reproducibility-chapter:

Reproducing a stated result
===========================

Reproduction starts by choosing a proposition and its exact edition. A
saved-certificate check, regeneration of a finite search, a conditional Lean
proof and a PDF build establish different things. Their receipts retain
those differences even when they contribute to the same paper.

The combined submission brings the finite-volume 57, compact 17 and compact
arithmetic 15 proof families together in 551 pages, with 71 pages of
main exposition followed by complete proof appendices. Its compact theorem is
an internal input to its arithmetic argument.

The two right-angled submission editions remain separate current
reading copies: six pages with two figures for sharp dimensions 8/4,
and eighteen pages with five figures for ideal volume and structure. Their
isolated TeX builds have zero warnings. These presentation builds retain the
recorded mathematical and conditional Lean evidence; they do not rerun it.
See :ref:`repro-arxiv-editions` for the exact companion selections.

This page explains the available execution routes. The mathematical
contracts are developed in :doc:`algorithms`; :doc:`scope` records the
geometric, literature and finite-data assumptions that those routes leave
explicit. Successful unchanged jobs need no repeat for a documentation or
illustration update.

.. _repro-illustration-identity:

Source assembly and current editions
--------------------------------------

The combined dimension article includes the complete arguments and appendices
of its three proof families. Its source-preservation review checks their
assembly, label changes and internal references. The
:download:`current manifest <../../results/arxiv_submission_current_manifest_20261030.json>`
identifies the exact submission sources, matching PDFs and isolated build
receipts. The :download:`dimension guide <../../arxiv/DIMENSION_BOUNDS.md>`
links the retained proof inputs and mathematical maps.

Source preservation and PDF compilation do not establish a fresh mathematical
replay. The catalogue and evidence records below retain their declared source
bindings and scopes.

.. _repro-edition-selection:

Choose the edition before the command
-------------------------------------

Run commands from the repository root. The
:download:`catalogue <../../manuscripts/catalogue.json>` binds an identifier
to a particular TeX source, script selection and conditional Lean scope.
An identifier does not automatically follow a later working article.

.. list-table:: Primary catalogue selections
   :header-rows: 1
   :widths: 30 31 39

   * - Identifier
     - Recorded edition
     - Important boundary
   * - ``finite57``
     - Frozen finite-volume 57 article
     - The current application conjunction and illustrated source have separate records.
   * - ``compact_farthest_bounds``
     - Frozen compact 23 article
     - The compact-17 appendices have separate source and verification records.
   * - ``compact20``
     - Frozen compact arithmetic 20 article
     - It does not replay every premise and case of the arithmetic-15 article.
   * - ``right_angled_sharp``
     - Sharp right-angled 8/ideal 4 article
     - Selected exact and conditional logic checks retain their listed inputs; ``--pdf`` typesets the preserved five-page edition.
   * - ``ideal_ra4_volume_and_structure``
     - Ideal RA4 volume and structure article
     - Restricted chain recognition remains restricted; ``--pdf`` typesets the preserved seventeen-page edition.
   * - ``compact_ra4``
     - Compact RA4 constraint article
     - Its checks do not settle unrestricted volume minimality.

Inspect the available identifiers and a proposed selection without
launching mathematical children:

.. code-block:: sh

   python3 -B scripts/run_proof_batch_20261024.py manuscript --list
   python3 -B scripts/run_proof_batch_20261024.py manuscript finite57 \
     --mode full --pdf --plan
   python3 -B scripts/run_proof_batch_20261024.py manuscript compact_farthest_bounds \
     --mode core --pdf --plan
   python3 -B scripts/run_proof_batch_20261024.py manuscript compact20 \
     --mode core --pdf --plan

The modes select scripts, conditional logic, their core combination or the
full declared combination; ``--pdf`` adds the recorded typesetting stage.
Read the returned source paths and coverage fields. A plan authenticates
the proposed inventory, not its execution. The guarded public entry also
rejects Python optimization: assertions must remain enabled, with neither
``-O``/``-OO`` nor an optimizing ``PYTHONOPTIMIZE`` setting.

.. _repro-arxiv-editions:

Reproducing the three submission editions
------------------------------------------

The :download:`28 October edition guide <../arxiv_submission_editions_20261028.md>`
and :download:`current source/PDF manifest <../../results/arxiv_submission_current_manifest_20261030.json>`
bind the current public scientific sources to their reading PDFs:

.. list-table:: Submission source and supporting suite
   :header-rows: 1
   :widths: 32 18 30 20

   * - Source and reading PDF
     - Pages / figures
     - Companion identifier
     - Script-only jobs
   * - :download:`Sharp-dimension TeX <../../arxiv/sharp_dimension_bounds_8_and_4/main.tex>`; :download:`PDF <../../manuscripts/build/arxiv_sharp_dimension_bounds_8_and_4_20261028_current_20261030.pdf>`
     - 6 / 2
     - ``right_angled_sharp``
     - Two exact model/algebra checks.
   * - :download:`Ideal-volume TeX <../../arxiv/ideal_four_dimensional_volume_and_structure/main.tex>`; :download:`PDF <../../manuscripts/build/arxiv_ideal_four_dimensional_volume_and_structure_20261028_current_20261030.pdf>`
     - 18 / 5
     - ``ideal_ra4_volume_and_structure``
     - Four exact model/incidence/root checks.
   * - :download:`Combined-dimension TeX <../../arxiv/hyperbolic_coxeter_dimension_bounds/main.tex>`; :download:`PDF <../../manuscripts/build/arxiv_hyperbolic_coxeter_dimension_bounds_20261028_current_20261030.pdf>`
     - 551 / 9
     - Three existing proof families
     - :download:`Proof and computation maps <../../arxiv/DIMENSION_BOUNDS.md>`; no combined fresh replay claimed.

Each source is self-contained for typesetting with standard pdfLaTeX packages:
its figures and bibliography are inline. Each upload folder already contains
its ``main.tex``. Package the three folders separately; keep PDFs outside them.
The :doc:`papers` chapter identifies the proof maps for the current sources.

All inputs selected by the two right-angled script-only suites are small,
tracked files. Neither selection requires the separate unuploaded bulk
collection. The workers use Python 3.11 or later and only the standard library.
The guarded launcher at this snapshot uses macOS resource measurements;
printing the plans is portable, while actual guarded execution requires those
probes.

Inspect the jobs without launching mathematical checks:

.. code-block:: sh

   python3 -B scripts/run_proof_batch_20261024.py manuscript \
     right_angled_sharp --mode scripts --plan
   python3 -B scripts/run_proof_batch_20261024.py manuscript \
     ideal_ra4_volume_and_structure --mode scripts --plan

For a necessary new execution, choose ledger paths that do not already exist:

.. code-block:: sh

   python3 -B scripts/run_proof_batch_20261024.py manuscript \
     right_angled_sharp --mode scripts \
     --ledger results/sharp_replay_new.json
   python3 -B scripts/run_proof_batch_20261024.py manuscript \
     ideal_ra4_volume_and_structure --mode scripts \
     --ledger results/ideal4_replay_new.json

These are sequential model or algebra checks, with the scopes explained in
:doc:`right_angled` and :doc:`volume_problems`. The ``core`` selection also
includes conditional Lean stages; compilation precedes a separate actual axiom
audit where declared. Geometry, literature and finite-data interfaces remain
explicit premises. The ideal suite's ``full`` selection additionally offers
its historical double-model producer.

The optional catalogue ``--pdf`` job follows the source already bound to that
catalogue identifier. It compiles the earlier five- or seventeen-page reference
edition, not the current six- or eighteen-page submission. The new isolated
build records establish typesetting of the matching submission sources, not
a fresh combined execution of their Python and Lean companions. No successful
unchanged mathematical check was repeated merely for this presentation update.

For the combined dimension article, the source-preservation review covers the
three main texts and all three appendix groups. Existing script results and
conditional Lean kernels keep their original hypotheses and scopes. The
compact and arithmetic case proofs do not acquire additional replay credit by
being placed in one document. The source is self-contained for typesetting;
the :doc:`sources` chapter separately states computational data availability.

The :download:`edition guide <../arxiv_submission_editions_20261028.md>` gives
isolated build commands for all three sources. The short papers use
``scripts/build_arxiv_submission_20261028.py``; the combined article uses
``scripts/build_dimension_bounds_arxiv_20261028.py``. Both run through the
same serial resource-monitored entry. A PDF build launches no mathematical
checker or Lean module.

.. _repro-finite57-ledgers:

Finite57 has two different inventories
---------------------------------------

The historical manuscript inventory contains 136 scripts, nine Lean roots
and one PDF job: 146 jobs, with an unchanged core inventory of eleven.
The two formerly eager loader selections have complete bounded
replacements: all 25 distinct actual PASS units were integrated into the
primary catalogue. The other 111 script-command matches retain inherited
execution records. These are separate records, not a combined fresh
146-job runtime PASS.

The assertion-safe successor explicitly separates compilation from actual
axiom audits. Its recorded full plan has 151 stages: 136 script stages,
nine compilation stages, five separate audit stages and one PDF stage.
Its core plan has sixteen stages. Only five of the nine selected roots
contain the relevant print commands. Splitting those phases changes the
plan accounting; it supplies no new proof execution.

The current-rule conjunction is another, completed, scope. It joins 195
prerequisite, control and index units to 64 actual selected application
units, for 259 actual scoped units. Their coverage is 28,028 maps: 15,510
direct units and 12,518 nonunit applications, with 83 prefixes and 24
terminal-tail applications. All 64 selected shards and the final reader
have actual PASS records at the explicitly registered ordinary and
supplemental premises.

The evidence homes are the
:download:`bounded-loader integration <../../results/finite57_complete_loader_catalogue_integration_20261022.json>`,
the :download:`current-rule completion <../../results/finite57_current_rule_release_completion_20261025.json>`,
the :download:`all-64 application record <../../results/finite57_longer_tail_all64_actual_application_completion_20261025.json>`
and the :download:`final reader execution <../../results/finite57_final_current_rule_conjunction_execution_20261025.json>`.
The :download:`finite57 reproduction guide <../finite_volume57_reproduce.md>`
states the current algorithm contract, source pins, failure history and
exact commands within one guide.

.. _repro-individual-obligations:

What an individual algorithm must establish
--------------------------------------------

A command list is not a completeness argument. Before using an output,
identify its input domain, the inclusion of every required geometric case,
the meaning of each accepted or rejected branch, termination and the exact
conclusion consumed by the proof. The main instances require the following
different arguments.

.. list-table:: Mathematical obligations behind the replay routes
   :header-rows: 1
   :widths: 23 45 32

   * - Algorithm or certificate
     - Completeness and termination obligation
     - Meaning of the checked output
   * - Additive forest tables
     - Every component has positive rank. Every reachable last-component transition is dominated, and every finite positive-rank entry has a tight predecessor. Context interactions must be fixed or retained in the state before additivity is used.
     - An exact table maximum at the specified alphabet and rank. A single maximizing witness proves only attainability.
   * - Ideal-link terminal profiles
     - Sorted positive-part recursion lowers the remaining sum. The concentration proof first reduces every required profile to the terminal domain. Every product-polynomial coefficient is reconstructed by exact convolution.
     - The stated curvature inequalities throughout that terminal domain, with the original fiber multiplicities retained when lifted to actual links.
   * - Native labelled maps and selected applications
     - Backtracking tries every unused target vertex and preserves all labelled edges and nonedges. Classification and order-cap premises justify the target domain. Shards and prefix checks must cover their exact physical keys.
     - The registered map and application coverage. Reading their receipts does not regenerate the native-map census.
   * - Arithmetic finite outputs
     - Full maximal-order strips or proved cubic coefficient boxes precede enumeration. All embeddings, required signs, cycles and field-square identities are retained. Finite integer boxes and exact decision predicates justify termination.
     - The four named output premises of arithmetic 15. The B/D17 conjunction rechecks survivors and counters at source-bound batch coverage; it does not independently re-enumerate every rejected trial.
   * - Compact forest covers
     - Every required rank partition enters the necessary cover; removal uses a proved energy inequality for the same original vector. Proper forests and whole classical components need their own consumers.
     - Necessary alternatives for geometric proofs, rather than realized examples or a census of polytopes.
   * - Right-angled exact models and identities
     - The checked model, combinatorial domain or algebraic identity must match the named proposition. Its written geometric inclusion and cited equality or recognition theorem remain separate.
     - The precise model or identity check, rather than unrestricted geometric classification.

See :ref:`alg-five-obligations` for the common definition and the subsequent
sections for the induction and domain arguments. The manual whole-D16
dyadic exclusion in arithmetic 15 is a written proof, not an executed D16
census. A timeout or a memory refusal never certifies an empty search.

.. _repro-input-availability:

Check input availability before replay
---------------------------------------

The source repository contains manuscripts, code, final PDFs and selected
certificates and receipts. The separate
:download:`bulk manifest <../../repro/public-data-manifest.json>` describes
936 scientific files in 33 prepared archives. **No bulk asset has been
uploaded.** The restore interface therefore does not yet provide a
demonstrated public restoration of those inputs.

Read the :download:`bulk-data guide <../../repro/PUBLIC_DATA.md>` and inspect
the local plan:

.. code-block:: sh

   python3 -B scripts/fetch_public_data.py --plan

This planning mode makes no network requests and writes no restored files.
The existence of a manifest, an upstream release URL or a local research
archive does not make every replay input available from a source checkout.
Do not substitute a newly generated expected answer for a missing pinned
certificate during verification.

When the exact authorized archives are available, restoration can select
one or more literal manifest paths. For already downloaded archives, the
publication-only watchdog route is:

.. code-block:: sh

   python3 -B tools/run_publication_io.py -- \
     python3 -B scripts/fetch_public_data.py --restore \
     --from-directory release-assets \
     --path results/fresh_ambient57_payment_frontier_20261018.json

Repeat ``--path`` for additional exact paths; it does not accept shell
globs. Omitting ``--from-directory`` selects the configured download route,
which requires the assets actually to exist. Omitting ``--path`` selects
all missing manifest entries. The helper verifies archive and restored-file
hashes, refuses conflicting existing files and restores sequentially.
Its own implementation also retains an additional 20-GiB free-space
refusal and checks the space needed for the selected restoration. That
tool-specific check does not replace the current percentage policy below.
Downloaded archives are normally removed after use; ``--keep-cache``
retains them.

The :download:`external-input guide <../../repro/EXTERNAL_INPUTS.md>` and
:download:`external inventory <../../repro/external-inputs-manifest.json>`
record 37 excluded third-party files. The bulk helper does not restore
them. Their upstream references and local hashes support deliberate
retrieval, not a promise that the same bytes remain obtainable or a
redistribution licence. The six primary papers' selected verification
plans do not require those excluded downloads; historical research routes
may do so.

Successful restoration proves byte identity. Mathematical completeness,
geometric interpretation and any conditional literature inputs still
need their own arguments. Neither these instructions nor a local cleanup
receipt establishes complete clean-checkout replay.

.. _repro-necessary-new-run:

Run only the new scope that is needed
-------------------------------------

Use fresh output paths, retain failed attempts and resource deferrals, and
keep the actual command, source identities, input snapshot, stdout and
monitoring ledger together. A changed source needs its own review; an old
PASS must keep its original source binding.

For a substantively needed new finite57 final-reader replay, after checking
its required local inputs, the guarded entry is:

.. code-block:: sh

   python3 -B scripts/run_proof_batch_20261024.py check \
     --ledger results/finite57_final_current_rule_new.json -- \
     scripts/check_finite57_current_rule_conjunction_20261025.py \
     --batch final --parent-lock \
     > results/finite57_final_current_rule_new.stdout

This reads the registered conjunction and checks its application evidence.
It does not rerun all native indexes, all 64 application jobs or all
inherited inputs. Its exact implementation and helper are pinned by the
current contract. Representation safeguards include at most 320 actual
units, depth 48, 8 MiB per expanded payload, 128 MiB aggregate expanded
payload, and separate 256-MiB source-hashing and parser-read budgets.
Exceeding a cap refuses the input rather than truncating it and reporting
complete coverage. A successful execution establishes the exact admitted
snapshot; it does not guarantee successful completion for arbitrary
modified inputs.

For a necessary catalogue-core replay, an example is:

.. code-block:: sh

   python3 -B scripts/run_proof_batch_20261024.py manuscript right_angled_sharp \
     --mode core --pdf --ledger results/right_angled_core_new.json

Review its plan first. A core success covers that selected core, including
its declared conditional interfaces. It does not add unselected scripts,
new working-paper appendices or external geometry to its scope.

.. _repro-resource-policy:

Resource gates and serial ownership
------------------------------------

Scientific, Lean and TeX jobs share one active process and
``formal/.proof-batches/exclusive.lock``. The current policy is the user's
25 October percentage override. Let :math:`D` be total filesystem
capacity, :math:`d` free bytes, :math:`M` physical RAM and :math:`m`
available RAM. Admission requires

.. math::

   100d\geq2D,\qquad100m\geq15M.

Stop an active job if :math:`100d<D` or :math:`100m<6M`. Below 2% free
disk, start cleanup before further jobs. Decisions use exact byte ratios,
not rounded displays. Receipts preserve total and available bytes,
percentages and the resource-policy source hash. The former scientific
40/6-GiB admission rules remain historical in their dated records.

Scientific and TeX work uses the existing bounded research runner with
59 seconds, 1 GiB and one thread. Lean uses one named module with
2048 MiB and one thread through
``scripts/run_lean_batches_signal_safe_20261018.py``, followed by a
separate actual axiom audit. Keep process-tree and RSS monitoring enabled.
If monitoring access is denied, obtain the required execution permission
before launching; do not bypass it.

Publication-only I/O uses ``tools/run_publication_io.py`` and its
allowlist. It samples the owned process group's RSS against 256 MiB and
uses the same percentage resource gates. This is sampled monitoring,
rather than an enforced heap-allocation cap, and is separate from the
proof runner. It cannot serve as an alternative mathematical launcher.

After a batch remove only owned, unused build products. Preserve sources,
certificates, final PDFs, receipts, failed attempts and history. Check an
exact byte round trip when compressing preserved material. Never kill
unrelated processes or launch a parallel proof farm.

.. _repro-lean-two-stages:

Compile Lean, then audit the actual target
------------------------------------------

For a necessary new run, the public stages for one compact scalar module
are separate:

.. code-block:: sh

   python3 -B scripts/run_proof_batch_20261024.py lean-compile \
     CoxeterProofs.D18UniformDottedHeadScalarKernel \
     --ledger results/d18_needed_compile_new.json
   python3 -B scripts/run_proof_batch_20261024.py lean-audit \
     CoxeterProofs.D18UniformDottedHeadScalarKernel \
     --ledger results/d18_needed_axioms_new.json

The audit requires matching compiled-target, source and import
fingerprints, then separately invokes Lean to print the declarations. It
does not silently rebuild a missing target. Cleanup of rebuildable objects
does not erase the preserved successful audit, but a new actual audit
still needs the correct target.
Avoid a raw whole-project Lean build and do not rerun successful unchanged
modules merely to refresh documentation.

The :doc:`formal_logic` chapter shows what these kernels assume. Ordinary
theorem parameters can encode geometry, a complete case cover or an exact
computed output without appearing in the printed axiom list. A successful
audit of standard foundational axioms therefore remains a conditional
internal proof, not a formalization of every geometric input.

.. _building-only-this-site:

Building only this site
-----------------------

The documentation environment is separate from the proof environment.
The :download:`documentation requirements <../../repro/requirements-docs.txt>`
pin Sphinx 9.1.0; the site uses its built-in Alabaster theme and MathJax.
Proof dependencies are instead pinned in
:download:`proof requirements <../../repro/requirements-proof-20261018.txt>`;
Lean's version and package configuration are described in
:doc:`formal_logic`.

In a suitable documentation environment, the local strict build is:

.. code-block:: sh

   python3 -m pip install -r repro/requirements-docs.txt
   python3 -B tools/run_publication_io.py -- \
     python3 -B tools/audit_sphinx_inputs.py --require-indexed
   python3 -B tools/run_publication_io.py -- \
     python3 -m sphinx -W --keep-going -b html docs/site docs/site/_build/html

The input audit checks publication boundaries, including files exposed by
download directives. A locally present bulk file must not be republished
through Sphinx while its payload remains unuploaded and unauthorized.
The Sphinx configuration launches no mathematical checks. A strict site
build proves neither a paper's mathematics nor its Lean interfaces.

Before every push, audit the actual Git index through the publication
watchdog:

.. code-block:: sh

   python3 -B tools/run_publication_io.py -- \
     python3 -B tools/audit_public_tree.py

Keep generated sites, Lean objects and dependencies, TeX intermediates,
Python caches and packaging output untracked. Intentional final PDFs and
verification stdout are retained evidence. The Pages workflow builds
documentation after an authorized push; it performs no proof checking,
TeX compilation or bulk-data restoration.

