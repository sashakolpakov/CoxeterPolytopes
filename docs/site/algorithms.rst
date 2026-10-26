.. _algorithms-chapter:

Finite algorithms, certificates and replay
==========================================

A computation enters these proofs through a mathematical proposition:
a specified domain, a complete algorithm and an exact output. The
implementation, certificate and execution receipt then identify which
instance of that proposition was checked. These layers answer different
questions.

The :download:`computational proof contract <../computational_proof_contract.md>`
gives the common requirements. This chapter develops the principal examples
from :doc:`finite_volume` and :doc:`compact_bounds`, then explains the
recorded execution boundaries and guarded commands.

.. _alg-five-obligations:

Five obligations before a result is used
----------------------------------------

Domain
   Specify the objects, exact arithmetic, bounds and equivalence relation.
   Prove that every geometric object needed in the argument reaches this
   domain. If a label or rank starts unbounded, retain a symbolic family
   or prove its reduction; a convenient experimental cutoff is insufficient.

Soundness
   Every accepted row, transition or survivor has the stated meaning.
   Every rejected branch must fail a necessary condition, or belong to a
   branch already handled by a uniform proof.

Completeness
   Every admissible object is visited or represented by a justified
   quotient. A symmetry reduction must preserve the marked roles and all
   predicates used later. An example of a successful completion does not
   show that all completions were considered.

Termination
   Give a finite state space, bounded loops, a decreasing integer or a
   convergent exact decision procedure with a proof that it eventually
   decides. A resource timeout is a refusal to finish that execution, not
   a proof that its intended mathematical search is complete.

Computed output
   State the exact maximum, list of exceptions or family of inequalities
   needed by the theorem. Preserve the input, implementation, complete
   certificate and actual receipt. An unexecuted plan supplies none of
   these output assertions.

The geometric inclusion into a finite domain is particularly important.
A correct exhaustive search of the wrong domain proves the wrong
proposition. Conversely, an excellent geometric reduction with an
unexecuted finite output remains conditional on that output.

.. _alg-forest-certificate:

Exact additive forest maximization
----------------------------------

Let :math:`\mathcal A` be a finite alphabet of component types. Each
:math:`T` has positive integral rank :math:`r(T)` and exact weight
:math:`w(T)`. For a fixed total rank bound :math:`N`, consider all
multisets of these types. Repetition is permitted.

The recurrence is

.. math::

   V(0)=0,\qquad
   V(s)=\max_T\{V(s-r(T))+w(T)\},\quad1\leq s\leq N,

where a candidate is permitted only when :math:`r(T)\leq s` and its
predecessor is reachable. An empty maximum is an explicit unreachable
state, not an ordinary negative number.

Why the recurrence is exact
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Evaluate in increasing rank. Every candidate is attainable by appending
its last component to an attaining predecessor. Conversely, remove any
component :math:`T` from a nonempty multiset of rank :math:`s`. By induction
its remaining cost is at most :math:`V(s-r(T))`, so its full cost is at
most a candidate examined by the recurrence. These two inequalities prove
the exact maximum and the reachability statement together.

There are finitely many ranks and finitely many candidates at each rank.
Every predecessor has smaller rank because :math:`r(T)>0`. Thus the
algorithm terminates, and following a stored predecessor eventually
reaches the rank-zero base.

What a saved table must certify
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A certificate claiming an upper bound must dominate *every* admissible
last-item transition, with reachability reconstructed from the same
alphabet. A certificate claiming the exact maximum must additionally
provide a tight predecessor for every positive-rank finite entry.
Rank zero is reachable with value zero.

Only checking the chain of one maximizing witness would prove
attainability, not an upper bound for omitted multisets. Only checking
upper inequalities could permit a table with unnecessarily high finite
values. The transition inequalities and tight witnesses have different
jobs and together prove equality.

.. _alg-contexts:

Context-dependent costs need more state
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The actual vertex rules can include pair interactions, multiplicity
denominators, presence bonuses and guarded replacements. Such a cost is
not additive merely because its rows are stored in a table. The proof
first fixes all interacting counts, or retains enough state to determine
every incremental cost.

In finite57, multiplying by 64 produces integer tables. Context-sensitive
counts are fixed before the additive residual problem is used. Candidate
classes capable of exceeding the target are expanded into their actual
type fibers. Counts of small components and the remaining presence guard
are retained. A presence bonus is paid once; reserving a first copy does
not authorize that bonus for every later copy.

Compact forest covers have a related but different output: they enumerate
necessary rank partitions and discard a forest only by a proved
same-vector energy inequality. Their retained entries are hypotheses for
later geometric consumers. They are not constructed examples or a census
of realizable polytopes.

.. _alg-partitions:

Sorted partitions and ideal-link polynomials
--------------------------------------------

Let :math:`\mathcal P(n,\ell)` denote nondecreasing partitions of
:math:`n` with parts at least :math:`\ell`. The recursive definition is

.. math::

   \mathcal P(0,\ell)=\{()\},\qquad
   \mathcal P(n,\ell)=
      \bigcup_{q=\ell}^{n}\{(q,p):p\in\mathcal P(n-q,q)\}.

Every call lowers the remaining sum by a positive integer. The unique
first part of a sorted partition, followed by its unique tail, proves
both completeness and absence of duplicates by induction.

For finite57's terminal cusp calculation, choose an optional exceptional
part :math:`\ell\in\{0,\ldots,7\}`, omit zero, and generate the rest
from :math:`\mathcal P(57-\ell,8)`. The exceptional part determines its
outer branch uniquely. This gives the exact terminal domain justified by
the concentration argument in :ref:`fv-positive-concentration`.

For each profile, start with coefficient list :math:`(1)`. Adjoining a
rank :math:`r` replaces the list by convolution with
:math:`(1,\ldots,1)` of length :math:`r+1`:

.. math::

   h^{(k+1)}_j=\sum_{a=0}^{r}h^{(k)}_{j-a},

where coefficients outside the previous range are zero. Induction proves
that the list is the full product polynomial, not merely a degree or
symmetry invariant. Integer coefficients and rational curvature formulas
then give exact margins.

The recorded 2,552-profile proposition is checked by reconstructing the
domain and every convolution. A check of a few attaining profiles,
polynomial symmetry or aggregate count would have a narrower scope.

.. _alg-induced-maps:

Complete labelled induced maps
-------------------------------

The current finite57 route domain fixes a physical high-core row, its
ordered base and selected marks, and its accepted whole-component cap.
For the long-component branch, the frozen route algorithm considers
finite :math:`A,B,D` and affine :math:`A,B,C,D` targets in every allowed
order, at most seventeen. The classification and order-cap premises
justify why this target list covers the required geometric alternatives.

Backtracking assigns each pattern vertex to an unused target vertex.
Before descending, it compares the new assignment with every earlier
one, checking both the required labelled edges and the required nonedges.
A degree filter is used only as a necessary condition for an induced
injection. Assigning the vertices in a different fixed order changes
efficiency but not the domain.

At depth :math:`d`, every compatible assignment of the first :math:`d`
vertices is retained. Trying all permitted images of the next vertex
proves the corresponding statement at depth :math:`d+1`. Deleting the
last assignment proves that every complete injection occurs. Strictly
increasing depth bounded by the pattern size proves termination.

The output key contains all nine physical/native fields described in
:ref:`fv-actual-coverage`. The current application reader reuses the
actual route-index receipts; it does not rerun this enumeration. Direct
units, nonunit applications and physical prefixes retain separate
coverage checks.

.. _alg-arithmetic-domains:

Full integer-ring domains and exact field tests
------------------------------------------------

Quadratic strips bound every coefficient
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For squarefree :math:`d>1`, use the maximal-order basis

.. math::

   \omega=\begin{cases}
     (1+\sqrt d)/2,&d\equiv1\pmod4,\\
     \sqrt d,&\text{otherwise}.
   \end{cases}

If :math:`Y=a+b\omega` generates the quadratic field, has
:math:`Y>Y'>0` and :math:`Y<B`, then
:math:`0<b\sqrt\Delta=Y-Y'<B`. Thus :math:`\Delta<B^2`, and the
strict strips give finite integer intervals for :math:`a,b`. This proves
coverage of the full order, including half-integral elements. Squarefree
trial division and exact integer square roots supply finite generation
and floor operations.

For an element :math:`R=a+b\sqrt d`, a field-square test uses both signed
norm branches. When :math:`b\ne0`, a square must have rational candidates

.. math::

   x^2=(a+n)/2\quad\hbox{or}\quad x^2=(a-n)/2,
   \qquad n^2=\mathrm N(R).

For a nonzero rational square candidate, set :math:`y=b/(2x)` and verify
the full identity :math:`(x+y\sqrt d)^2=R`. The rational case and zero
are handled separately. The final identity check and physical positive
sign are part of the predicate; a square rational norm alone is not a
certificate of a square in the field.

Cubic boxes precede the search
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two finite-anchor :math:`A_{17}` cubic domains use their proved full
integer rings. Rational isolating intervals locate each embedding.
Lagrange interpolation of the bounded embedding values gives finite
integer boxes for all three coefficients of
:math:`V=a+bU+cU^2`. Negative coefficients are retained. The box proof
comes before the trial loop, so it does not infer its limits from the
models that happened to survive.

Arithmetic reduces modulo the cubic polynomial. A nonzero element has a
representative of degree at most two, so it cannot vanish at a cubic
generator. Interval refinement of a separate root bracket therefore
eventually determines its sign; zero is first recognized by its three
zero coefficients. Field-square decisions retain the bounded trace and
both signed norm possibilities, then verify the square identity and the
physical branch.

The article states the four outputs in :ref:`arithmetic-four-outputs`.
Their replay scopes are deliberately unequal: some conjunctions revisit
every trial, while the B/D17 conjunction rechecks survivors and counters
at complete source-bound batch coverage. Neither same-code replay nor a
receipt parser is an independent second implementation.

.. _alg-evidence-levels:

Reading a certificate and an execution ledger
---------------------------------------------

The following distinctions prevent an identity check from acquiring a
larger mathematical meaning than it has.

.. list-table:: Different evidence answers different questions
   :header-rows: 1
   :widths: 27 38 35

   * - Evidence
     - What it can establish
     - What still needs its own argument
   * - Mathematical domain theorem
     - Every required actual case maps into a specified finite or symbolic domain.
     - The exact output and execution of its implementation.
   * - Complete saved-certificate check
     - The named certificate satisfies the checked transition, identity or application predicates.
     - Unchecked ancestors and geometric interpretations.
   * - Regeneration
     - The named producer recomputed its stated domain and output at its inputs.
     - Completeness of a domain not proved by that producer and inputs outside its scope.
   * - Source hash
     - Identity of the bytes that a record names.
     - Truth of the mathematical proposition in those bytes.
   * - Parser/control fixture
     - A particular schema, association or refusal behavior on its fixtures.
     - Actual scientific applications or an exhaustive geometric theorem.
   * - PDF build
     - The stated source produced a typeset document with the recorded warnings.
     - Mathematical correctness or inspection of every rendered page.

For finite57, an actual unit binds its exact command, worker source,
bounded PASS ledger, public stdout and input snapshot. The public output
also records the owned-process monitor and enabled assertions. The final
reader requires supported exact inner statuses; a similarly named status
or a source-only clearance cannot substitute for execution.

The selected and final snapshots are separate. Registering the selected
receipts in the final snapshot must not rewrite the input consumed by the
earlier selected jobs. Selected ordinary premises remain canonically
identical in the final record, and shared expanded source contracts must
agree. Earlier failed or interrupted attempts remain separate evidence.

.. _alg-fail-closed:

Complete output or explicit refusal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The current reader bounds its actual-unit inventory by 320. It has a
depth-48 JSON boundary, an 8 MiB complete expanded-payload boundary,
a 128 MiB aggregate expanded-payload boundary, a 256 MiB source-hashing
budget and a separate 256 MiB parser-read budget. The frozen parser also
has its semantic-visit cap, and graph procedures retain their fuel cap.

These are representation and execution safeguards, not a new mathematical
cutoff on the intended geometric domain. Reaching a cap refuses the
input; it does not truncate a list and report complete coverage. Large
complete results use a pinned lossless encoding with length, digest and
end-of-stream checks. Byte and object round trips establish that the
encoded result is complete, not that its propositions are true.

The actual recorded PASS receipts demonstrate that their exact snapshots
completed within these bounds. They do not guarantee that arbitrary
modified inputs will complete under the same caps.

.. _alg-necessary-replay:

Inspecting and running a necessary replay
-----------------------------------------

Run commands from the repository root. Successful unchanged computations
and Lean modules need no repeat to update a date or documentation. Before
a necessary new run, identify its exact scope, confirm required inputs and
choose output paths that do not already exist.

Inspect the historical manuscript selection without executing proof jobs:

.. code-block:: sh

   python3 -B scripts/run_proof_batch_20261024.py manuscript finite57 \
     --mode full --pdf --plan

The full selection has 136 scripts, nine conditional Lean roots and one
PDF job. Its historical 146-job catalogue is separate from the current
259-unit receipt/application conjunction. A core selection is smaller
and does not replay all selected script obligations.

To inspect the smaller plan, use:

.. code-block:: sh

   python3 -B scripts/run_proof_batch_20261024.py manuscript finite57 \
     --mode core --pdf --plan

When a new final-reader replay is substantively required and its inputs
are available, the documented public entry is:

.. code-block:: sh

   python3 -B scripts/run_proof_batch_20261024.py check \
     --ledger results/finite57_final_current_rule_local_replay.json -- \
     scripts/check_finite57_current_rule_conjunction_20261025.py \
     --batch final --parent-lock \
     > results/finite57_final_current_rule_local_replay.stdout

Those filenames are examples of fresh destinations. Replaying this reader
checks the registered conjunction; it does not rerun the 64 native route
indexes, all 64 application jobs or every ancestor calculation. The
:download:`finite57 reproduction guide <../finite_volume57_reproduce.md>`
specifies the source identities and current availability boundary.

.. _alg-resource-policy:

Serial resource gates
~~~~~~~~~~~~~~~~~~~~~

One scientific, Lean or TeX process runs at a time under
``formal/.proof-batches/exclusive.lock``. Scientific and TeX work uses the
bounded research runner with 59 seconds, 1 GiB and one thread. A new
batch requires at least 2% free filesystem capacity and 15% available
physical RAM. An active job stops below 1% free disk or 6% available RAM.
Cleanup starts below 2% free disk and removes only owned unused products.

Admission compares exact integer byte ratios, not rounded percentages or
the historical fixed GiB thresholds. Receipts record the total capacities,
available bytes, displayed percentages and policy-source identity.
Process-tree/RSS monitoring remains required; inaccessible monitoring is
an execution-permission issue, not a reason to bypass the watchdog.

Sources, certificates, final PDFs and verification history are preserved.
Compression of preserved material requires a byte-exact round trip.
Removing rebuildable objects does not remove their source, but a later
axiom audit still needs a matching compiled target.

.. _alg-lean-stages:

Lean compilation and an actual axiom audit are separate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lean uses one named module, 2048 MiB and one thread through the signal-safe
bounded runner. After successful compilation, a separate actual axiom audit
reads the matching target. For a necessary new run of the compact scalar
module, the two public stages are:

.. code-block:: sh

   python3 -B scripts/run_proof_batch_20261024.py lean-compile \
     CoxeterProofs.D18UniformDottedHeadScalarKernel \
     --ledger results/d18_needed_compile_new.json
   python3 -B scripts/run_proof_batch_20261024.py lean-audit \
     CoxeterProofs.D18UniformDottedHeadScalarKernel \
     --ledger results/d18_needed_axioms_new.json

The compile stage ultimately uses
``scripts/run_lean_batches_signal_safe_20261018.py``. The audit stage
does not silently rebuild a missing target. These commands are instructions
for a necessary future replay; the successful unchanged module already has
its actual compilation and separate audit records.

Printed foundational axioms and ordinary theorem parameters answer
different questions. A kernel may use only the expected foundational
axioms while still taking geometry, exact inequalities, complete
inventories or a computed output as parameters. Its logic map, not the
printed axiom list alone, describes that conditional scope.

.. _alg-availability:

Byte restoration and public availability
-----------------------------------------

The source repository includes code, proof texts, final PDFs and selected
evidence. The bulk-data manifest describes additional scientific inputs,
but no bulk asset has been uploaded. A clean source checkout therefore
does not supply every input required by every replay route.

The :download:`public-data guide <../../repro/PUBLIC_DATA.md>` explains the
manifest and restoration interface. The
:download:`external-input register <../../repro/EXTERNAL_INPUTS.md>` identifies
excluded third-party downloads and their retrieval boundaries. Restoring
matching bytes establishes identity; it does not prove the completeness or
geometric interpretation of a dataset. No redistribution licence is implied.

The right-angled package guide and formal records also distinguish removed
rebuildable targets from preserved local legacy dependency bytes. A local
archive is not supplied automatically by a source checkout, and a restored
old object does not acquire a new source-to-object proof. A complete
clean-checkout replay has not been demonstrated by these documentation
chapters.
