Reproducibility without unbounded runs
======================================

Run commands from the repository root. Inspect a plan first; planning launches
no proof child. Use the paper IDs on :doc:`papers`.

.. code-block:: sh

   python3 -B scripts/verify_manuscript.py right_angled_sharp --mode core --pdf --plan
   python3 -B scripts/verify_manuscript.py finite57 --mode full --pdf --plan

The second command is an inventory only. **Keep** ``--plan`` on finite57's full
selection until the two eager loaders have complete bounded replacements and
the catalogue has been updated. The current replacement has 3 of 25 actual
PASS units; a partial ledger is not a successful full replay.

Data and environment
--------------------

The public source tree and the bulk proof-data archive have separate roles.
Before replaying a check, restore the exact inputs required by its plan using
the manifest and :download:`public-data guide <../../repro/PUBLIC_DATA.md>`.
The :download:`bulk manifest <../../repro/public-data-manifest.json>` binds
the release assets and original bytes. Historical third-party downloads have
a separate :download:`retrieval guide <../../repro/EXTERNAL_INPUTS.md>` and
:download:`inventory <../../repro/external-inputs-manifest.json>`; the bulk
helper does not restore them.
Inspect the restoration plan first:

.. code-block:: sh

   python3 -B scripts/fetch_public_data.py --plan

The planned release assets must not be assumed available merely because these
instructions exist. A restore plan is not a mathematical replay.

Once the release assets are available, restore all missing manifest-listed
inputs sequentially, or select an exact path:

.. code-block:: sh

   python3 -B scripts/fetch_public_data.py --restore
   python3 -B scripts/fetch_public_data.py --restore \
     --path results/fresh_ambient57_payment_frontier_20261018.json

Repeat ``--path`` to select several inputs; it does not accept shell globs.
For archives already downloaded into a local ``release-assets`` directory:

.. code-block:: sh

   python3 -B scripts/fetch_public_data.py --restore --from-directory release-assets

By default, downloaded archives are removed after use; ``--keep-cache`` retains
them. Restoration verifies byte identity and refuses conflicting existing
files. It does not verify the mathematical claims encoded by those bytes.

Proof-environment requirements are separate from the small Sphinx environment.
Use the :download:`proof-environment instructions <../../repro/README.md>` and
:download:`pinned proof requirements <../../repro/requirements-proof-20261018.txt>`.
Do not replace a missing certificate by regenerating an expected result during
verification.

Admitted sequential checks
--------------------------

After reviewing the plan, restoring its required data, and satisfying the
resource gates, the sharp right-angled core can be run serially with a new
ledger:

.. code-block:: sh

   python3 -B scripts/verify_manuscript.py right_angled_sharp --mode core --pdf \
     --ledger results/right_angled_core_new.json

The smaller finite57 endpoint/logic/typesetting core is also distinct from
the blocked full script selection:

.. code-block:: sh

   python3 -B scripts/verify_manuscript.py finite57 --mode core --pdf \
     --ledger results/finite57_core_new.json

Choose unused ledger filenames. Preserve failures and resource deferrals;
do not overwrite old receipts or silently credit a success under a different
source pin. The finite57 core does not replay the 113 script obligations.
Its :download:`reproduction guide <../finite_volume57_reproduce.md>` supplies
the completed bounded parts and their exact conjunctions. The shared
:download:`verification guide <../manuscript_verification.md>` explains
``scripts``, ``logic``, ``core``, and ``full``.

Resource contract
-----------------

Scientific, Lean, and TeX work shares **one process at a time** and the lock
``formal/.proof-batches/exclusive.lock``.

* Admission requires at least **40 GiB free disk and 6 GiB available RAM**.
* Abort thresholds are **8 GiB disk and 3 GiB available RAM**.
* Scientific and TeX children use the existing 59-second, 1-GiB, one-thread
  bounded wrappers. Lean uses one module, 2,048 MiB, and one thread, followed
  by a separate axiom audit.
* Keep process-tree/RSS monitoring enabled. If it is unavailable, defer the
  run; do not bypass the watchdog or lower a resource floor.
* A failed or deferred child stops the coordinated run. Remove only owned,
  unused temporary build products after a batch; retain sources, certificates,
  PDFs, and evidence. Below 20 GiB free disk, cleanup precedes further jobs.

No parallel proof farm, raw project-wide Lean build, or repeated successful
batch is needed to read the papers or build this documentation.

Building only this site
-----------------------

The documentation uses Sphinx 9.1.0 and its built-in Alabaster theme; it needs
neither MyST nor a separate theme package. In a suitable documentation
environment:

.. code-block:: sh

   python3 -m pip install -r repro/requirements-docs.txt
   python3 -m sphinx -W --keep-going -b html docs/site docs/site/_build/html

The Sphinx configuration launches no proof checks. A documentation build does
not compile a paper, replay mathematics, or verify a Lean module.

