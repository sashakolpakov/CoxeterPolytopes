# Reproducing the finite-volume57 paper

The current [manuscript](../manuscripts/finite_volume57.tex) proves dimension
at most57 for every finite-sided, finite-volume real hyperbolic Coxeter
polytope. Its current [48-page PDF](../manuscripts/build/finite_volume57.pdf)
matches source `c300db69…ca82a0` and PDF `1f3c020c…` in the
[build receipt](../results/six_paper_pdf_completion_20261020.json).
The [lemma map](finite_volume57_lemma_map.md) separates the written
geometry, finite certificates and conditional Lean results. No
arithmeticity, compactness, simple-cusp or facet-count assumption is used.

## Coordinated reproduction

From the repository root, inspect the selected checks before running them:

```sh
python3 -B scripts/verify_manuscript.py finite57 --mode full --pdf --plan
```

The current full selection has **123 jobs:113 scripts, nine conditional
Lean roots and one temporary PDF build**. The core selection has eleven
jobs: the endpoint checker, the same nine Lean roots and the PDF.
`--plan` runs no proof check. The saved
[shape-counting plan](../results/publication_shape_counting_verification_plan_20261020.json)
has the same ordered jobs and commands as the current catalogue.
It records an earlier catalogue hash; a plan is not runtime evidence.

The full selection is currently an inventory, not a safe complete replay:
its B13 and actual finite-frontier commands still call the eager loaders.
Their reviewed memory-bounded replacement must complete execution and
be integrated before removing `--plan` from the full command. The recorded bounded parts are
described below.

A [25-unit replacement](../scripts/check_finite57_b13_frontier_batches_20261021.py)
now has complete source-equivalence reviews, including correction of a
cached nested-input check. Its tiny parser unit has an actual bounded
PASS (0.139s,60,080,128B) in the
[runtime record](../results/finite57_loader_replay_runtime_20261021.json).
The [post-finishing successor ledger](../results/finite57_loader_replay_post_finishing_20261021.json)
reuses that preserved parser receipt and adds the actual `b13:pins-0`
PASS:14.154s and150,093,824B. The next unit was refused at4.6GiB available
RAM before a child started. The [next successor ledger](../results/finite57_loader_replay_quadratic_b16_20261022.json)
preserves all those records and adds the actual `b13:pins-1` PASS:
20.269s and130,203,648B, after admission at6.24GiB available RAM.
Three of25 units have runtime PASS records;22 remain without one.
The next unit is `b13:pins-2`. This prepared replacement is not yet the
catalogue selection, and neither producer has a new complete runtime PASS.
Its exact conjunction and assumptions are in the
[replacement contract](finite_volume57_lemma_map.md#prepared-bounded-b13-and-actual-frontier-replacement).

For the smaller endpoint, conditional-logic and typesetting selection,
run the core serially and save a new ledger:

```sh
python3 -B scripts/verify_manuscript.py finite57 --mode core --pdf \
  --ledger results/finite57_core_replay_new.json
```

Choose a new ledger filename if that path already exists. A core PASS
does not replay the 113 selected script obligations.
Each scientific child uses the existing shared lock, one thread and a
59-second/1GiB process-tree watchdog. Admission requires40GiB free disk
and6GiB available RAM; low resources defer work rather than lower the
floor. Lean uses its separate bounded runner. Process-monitor access must
be available; never bypass the watchdog. A failed or deferred child stops
the coordinated run. Preserve its ledger and resolve that specific issue
before resuming; a partial run is not complete verification.

The PDF check builds in a temporary directory, removes its intermediates
and preserves the delivered PDF and build ledger. No raw `pdflatex`
invocation is needed. The orchestration runs no optimizer and has no
parallel proof farm.

## Exact selected scope and existing evidence

These are preserved execution records, not new runs for this guide.

| Mathematical or technical obligation | Selected replay and actual evidence |
| --- | --- |
| Final finite377 closure; all2552 terminal ideal profiles, exact polynomial products and margins; all28 stop-loss transfers; constants377/36/4205 | The endpoint actually passed as the first job of the [original run](../results/finite57_publication_full_replay_20261020.json). It checks inherited interfaces and the eight old/six increasing forests; it does not rebuild the historical finite DP or native census. |
| Complete new D13 support, local interfaces and finite repair | All19 bounded parts and the current parser fixture pass in the [D13 ledger](../results/finite57_d13_batched_replay_runtime_20261020.json). The [equivalence contract](finite_volume57_lemma_map.md#memory-bounded-d13-replay-contract) covers every original comparison and pin obligation. |
| Coarse3298-choice forest screen and the759-row first increasing family | Both [ancestor phases](../results/finite57_ancestor_forest_replay_execution_summary_20261020.json) pass:8114 tables and10,195,241 transition inequalities. The [contract](../results/finite57_ancestor_forest_replay_contract_20261020.md) excludes exactly four separate actual-refinement fields from the coarse comparison. |
| Actual finite fibers, second increasing family, B13 repair and low8/higher-tail interfaces | Separate full-tier scripts are selected; these are certificate replays with the inherited boundaries given in the [interface audit](../results/review_six_manuscript_verification_interfaces_resource_20261020.md). A source hash alone does not replay an ancestor calculation. |
| Thirty original higher-tail commands | Actual PASS records are in the [tail ledger](../results/finite57_higher_tail_replay_batches_20261020.json) and its [remaining-check ledger](../results/finite57_higher_tail_remaining_replay_20261020.json). Their historical top statuses also preserve earlier interruptions. |
| The remaining original108-range obligation | The equivalent [distance9 replay](../results/finite57_distance9_partitioned_replay_runtime_20261020.json) passes all2196 required units in56 complete intervals, including overlap generation/classification and isolation. |
| Selected internal Lean logic | The [logic map](publication_logic_formalization_20261020.md) distinguishes actual axiom-audited kernels from legacy cached targets and lists explicit premises. Compilation does not instantiate the full geometric theorem. |

Follow-up: root's [separate low8 bounded replay](../results/finite57_audit_edition_low8_bounded_replay_20261021.json)
now passes in4.903s with462,585,856B peak RSS, checking2552 profiles and
the87/10240 normalized/36 absolute minima. With this addition,111 of the113
selected script commands have identified current-source bounded receipts.
B13 and the actual finite-frontier producer retain historical evidence;
their eager357MB serialized-input loads need reviewed partitioning before
safe current replay. They were not rerun. This is an execution-evidence
boundary, not a new mathematical gap.

**No combined123-job runtime PASS is claimed.** The original combined run
passed the endpoint and then exceeded the eager D13 loader's1GiB limit.
Its preserved failure is not a current D13 failure: the reviewed19-part
replacement now has complete actual coverage. The original eager D13
command is historical and is not a recommended reproduction entrypoint.

The full selection is also **not whole transitive regeneration**. Some
native local censuses and ambient59–65 forest/local-table computations
remain explicit inherited finite-data propositions, bound to their
sources and saved certificates. Actual face coverage, whole-factor and
fiber interpretation, convex-order geometry and the geometric higher
tail require the manuscript's written proofs. A successful hash, plan,
parser fixture, metadata coverage check or PDF build proves none of these
by itself. The [full source audit](../results/review_finite57_audit_edition_completion_20261021.md)
checks the current exposition at that declared boundary.

## Smallest endpoint and bounded parts

For only the final arithmetic composition:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py \
  --parent-lock -- scripts/fresh_d13_context_dimension57_20261018.py --verify
```

Its expected child headline is:

```text
PASS_COMPLETE_D13_CONTEXT_DIMENSION57_COMPOSITION dimension 57 finite 377 surplus 36 constant 4205
```

For one D13 part use the bounded replacement, not the original producer:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py \
  --parent-lock -- scripts/check_d13_transfer_batches_20261020.py \
  --verify --batch support
```

The [D13 contract](finite_volume57_lemma_map.md#memory-bounded-d13-replay-contract)
lists all19 required names. One part proves only its own comparison.

The [distance9 coordinator](../scripts/resume_distance9_replay_20261020.py)
can inspect the saved complete coverage or resume one admitted interval:

```sh
python3 -B scripts/resume_distance9_replay_20261020.py \
  --ledger results/finite57_distance9_partitioned_replay_runtime_20261020.json --plan
python3 -B scripts/resume_distance9_replay_20261020.py \
  --ledger results/finite57_distance9_partitioned_replay_runtime_20261020.json --max-jobs 1
```

The existing ledger is already complete, so same-source successes need no
repeat. A job-limit exit is not complete coverage; the coordinator checks
receipts rather than proving mathematics. Its
[contract](finite_volume57_lemma_map.md#resuming-the108-distance-nine-replay)
specifies all eight domains and unchanged gates. Generate no replacement
expected certificate to conceal a mismatch; use replay-only modes and
preserve frozen sources, data and reports.

## Mathematical certificate boundary

The finite forest algorithm fixes every interacting/context count, then
uses positive-rank additive Bellman transitions. Every stored finite
positive-rank entry needs both all-transition upper bounds and a tight
predecessor. Rank zero is the reachable base. Smaller predecessor rank
proves termination and exactness; context exhaustion and geometric
domination remain separate obligations.

The ideal checker independently enumerates sorted rank57 terminal
profiles, reconstructs every polynomial by convolution and evaluates exact
rational margins. Unique first parts prove domain completeness. The
written low8 descent terminates because each transfer strictly increases
the squared-rank sum, bounded by57². Its convex-order comparison is applied
only after positive terminal margins are established. The computed minima
are87/10240 normalized and36 absolute, and need not have the same witness.
The [algorithm contract](computational_proof_contract.md) and manuscript
give the full soundness, termination and completeness arguments.

The final global subtraction is

```text
D5 ≤ 377 vf + sum C_actual
D5 ≥ 377 vf + sum kappa + 4205
kappa − C_actual ≥ 36
therefore 0 ≥ 4205 + 36 vi > 0.
```

The two-extrema correction4205 is counted once and also excludes vi=0.
The completed higher tail handles every dimension at least59; no
unproved inheritance of a desired lower-dimensional Coxeter face is used.
Selected Lean kernels check internal implications with explicit
geometric/script inputs. No end-to-end geometric Lean verification,
external refereeing, dimension56 improvement or new example is asserted.
