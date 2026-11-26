# Changing and solving the global constraints

The scientific input contains the geometric and incidence conditions. Presolve
and conditioning produce an equivalent numerical representation, with exact
reconstruction maps and replayable proofs. Keep these artifacts distinct.

## Adding or enhancing a family

1. State the complete parameter domain and uniform map from every hypothetical
   polytope. Include unbounded ranks and labels, source hypotheses and the
   treatment of rank equalities. Selected diagrams do not establish coverage.
2. Use one moment registry per occupation measure. Only pointwise identities
   may be multiplied by that measure's state coordinates; aggregate incidence
   equations between measures are not pointwise identities.
3. Preserve all 26 feature coefficients, including contextual and cross-component
   terms. At fixed dimension scope, retain every parent equation and cone, or its
   exactly checked principal-submatrix image. Record current source-bound local
   lift regressions. Finite fixtures check implementation, not uniform coverage.
4. Write a family record containing `family_id`, `parameter_domain`,
   `uniform_lift_argument`, `source_references`, `moment_registry`,
   `implementation_modules`, `local_lift_receipt`, `rank_conditions_treatment`
   and `independence_status`. Run
   `scripts/register_global_constraint_change_20261125.py` on the parent and new
   inputs. Pass its result as `--registration` when preparing a changed model.
   This checks retention and test provenance, not the truth of a geometric proof.
5. Distinguish validity, independence and numerical usefulness. An old feasible
   witness may prove separation when it has no extension to the new moments.
   Lack of such separation is not proof of redundancy. Dimension/parity changes
   are separate uniform instantiations, not fixed-scope constraint extensions.

Independent source review and root acceptance remain required before changing a
geometric bound. Research exports can explicitly retain pending review status.

## Preparing solver input

Every scientific command uses the serial watchdog, resource limits and shared
lock in `AGENTS.md`. Use a fresh output directory and monitor for each batch.

1. Regenerate the scientific input. For the current model,
   `scripts/regenerate_global_constraints_20261125.py` compares the complete
   level-7 dimension-57 model with the frozen input. It does not rerun unchanged
   historical geometric verification jobs.
2. Run `scripts/prepare_conic_model_20261125.py --original INPUT --output PASS1`.
   An `EXACT_REDUCTION_CHECKPOINT_NOT_SOLVER_READY` result requires another batch
   with `--previous PASS1` and a fresh output. Continue until
   `presolve_fixed_point` is true. Reductions include exact positive scalar
   circuits, zero PSD rows, symbolic column kernels, constants, duplicate
   positive cone scalings and unused gauges. Face implications are propagated
   on the original cone rows until stable before rewriting the full matrices;
   each round is proved modulo only the preceding verified equations. Older
   one-round proof records remain replayable. Never erase coefficients by tolerance.
3. Run `scripts/verify_conic_reduction_chain_20261125.py --prepared FINAL_PASS
   --output VERIFIED`. It replays every source-cone proof from the frozen input
   and reconstructs every original coordinate, without trusting intermediate models.
   A batch replays two stages by default. If it returns
   `EXACT_REDUCTION_CHAIN_REPLAY_CHECKPOINT`, use the same final `--prepared`,
   `--previous REPLAY_CHECKPOINT` and a fresh output until the complete-chain
   status is recorded. Partial replays never authorize export.
4. Run `scripts/export_ready_global_sdp_20261125.py --prepared FINAL_PASS
   --verified VERIFIED --reference POINT.json.gz --output READY`. The reference
   selects units only and adds no feasible-point restrictions. Exact interval
   coordinates and positive column/block/diagonal units have checked inverse
   identities. Interval coordinates use sparse exact elimination and prefer
   response coordinates with equal-and-opposite coefficients, preserving shared
   populations when that option exists. Check homogeneous and overlapping
   intervals, including dependent selected pivots, with
   `scripts/check_global_conditioning_20261125.py`. Every affine entry is checked
   under the inverse transformation. Progress and an equivalent representation
   checkpoint are saved before the final audit; that checkpoint is not solver-ready.
   A new fixed-point audit runs after conditioning. Further exposed
   reductions block export rather than authorizing a partial model.
5. Only a matching `SOLVER_READY` receipt authorizes production SDPA. Use
   `scripts/run_sdpa_checkpoint_20261125.py --native-input READY/input.dat-s
   --initial READY/initial.dat --preflight READY/solver_ready.json --output RUN`.
   Continue later batches with `--resume RUN`, the same input and the same
   preflight. An iteration limit is a checkpoint boundary, not a decision.

The primal fixed-point audit does not establish strict feasibility of the
certificate cone. For an explicit certificate search, continue with these
separate serial jobs:

6. Encode the unchanged native input using
   `scripts/compare_conic_solvers_20261126.py encode`, with matching `--preflight`.
   SCS uses lower-column PSD vectors and Clarabel uses upper-column vectors;
   both require the square-root-of-two off-diagonal factor. The backend fixtures
   check these conventions against nonuniform three-by-three matrices.
7. Run `scripts/prepare_dual_face_model_20261126.py prepare --native-input INPUT
   --preflight PREFLIGHT --encoding ENCODING --exact-model EXACT --output FACE`.
   It discovers coordinate directions numerically but verifies every facial
   implication over the rational source. On each already proved face, a diagonal
   PSD coefficient direction forces the corresponding nonnegative dual weights
   to vanish. A zero PSD diagonal forces its entire row to vanish. Iterate to the
   recognized coordinate fixed point. Never delete a coefficient by tolerance.
8. Run the same script with `verify --prepared FACE` and a fresh output. This
   replays every exact implication and checks all reduced coefficients, cone
   embeddings and deleted zero equations. Preparation alone cannot authorize a
   certificate solve. This reduction preserves certificate feasibility; no
   primal feasibility equivalence or minimal-face claim follows from it.
9. Use `scripts/compare_conic_solvers_20261126.py solve --mode certificate` with
   both `--dual-face-model FACE` and `--dual-face-replay REPLAY/result.json`.
   The gate checks the separate replay and all its bound inputs and sources.
   The `--unreduced-certificate-diagnostic` option is an explicit benchmark
   exception, not the production path. Source or model changes require new
   matching preparation/replay records.

10. Reduce the certificate equalities over Q with
    `scripts/reduce_certificate_rational_20261126.py prepare`, supplying the
    exact prepared source, numerical encoding and separate dual-face replay.
    Use actual matrix entries and off-diagonal trace multiplicity two. Clear
    each equality's denominators and divide its integer content. The default
    eliminates only scalar columns occurring in one original equation, so
    reconstruction cannot fill other equations. Preserve any coupled core.
    Record numerator and denominator bit lengths, total map terms and surviving
    equations. A growth guard stops preparation; it does not authorize a solver.
11. Run the same script with `verify --prepared Q_MODEL` in a separate batch.
    Replay the triangular witnesses, check every source equation over Q, compare
    the full rational reconstruction, and check every exported numerical entry.
    This is equivalence of certificate feasibility. It neither decides that
    problem nor proves feasibility of the original necessary-data model.
12. For the current homogeneous model with one remaining normalization, run
    `scripts/balance_rational_certificate_20261126.py export --prepared Q_MODEL
    --verified Q_REPLAY/result.json`. This verifies the exact homogeneous
    reconstruction and uses positive units to balance the normalization. Keep
    its inverse when interpreting solver variables and infeasibility rays.
    Do not silently reuse the old native SDPA layout or its common-margin role
    for this explicit certificate problem.

13. To work on the normalized dual of that reduced certificate search, use
    `scripts/prepare_balanced_dual_20261126.py prepare` and a separate `verify`.
    The parametrization and exact trace calculation are in
    `notes/lp/global_formulation.md`, under “Dual certificate after exact equality
    reduction”. Every source coefficient, pivot inverse, cone partition and
    numerical entry is replayed. Its 61,513 coordinates enforce the dual linear
    identities before solving; only cone membership remains to be established.
14. Audit scalar row scales and PSD matrix-row scales, including constants.
    Bounded column norms do not ensure that small constraints are numerically
    resolved. `scripts/rescale_balanced_dual_20261126.py export` uses positive
    dyadic scalar units and repeated PSD diagonal congruences based on complete
    matrix-row envelopes. It checks every inverse and drops no coefficient.
    Readiness requires nonzero scalar and PSD matrix-row envelopes in `[0.5,2]`.
    This is a coefficient-scale condition, not strict feasibility or an
    eigenvalue bound. Diagonal-only scaling can enlarge off-diagonal entries;
    the retained first attempt illustrates why full matrix rows are used.
15. `scripts/solve_balanced_certificate_20261126.py` solves the balanced
    certificate search. `scripts/search_balanced_dual_field_20261126.py solve`
    solves its smaller explicit dual. The cone-rescaled variant uses
    `scripts/rescale_balanced_dual_20261126.py solve`. Record chordal-decomposition
    settings, actual arithmetic, all tolerances and returned arrays. The latter
    saves native scaled arrays and separately maps primal slacks and dual
    variables back to the original explicit-dual units.
16. Field reconstruction is a separate job. The explicit-dual point is checked
    by `scripts/search_balanced_dual_field_20261126.py attempt`, followed by
    `verify` if it produces a candidate. The opposite certificate side can be
    checked by `scripts/recover_balanced_certificate_field_20261126.py`, with its
    source role specified. A cone-rescaled seed must use the checked inverse
    cone/free-coordinate units. Kernel proposals restrict candidate recovery
    only. LLL height bounds and the tested quadratic fields are search choices,
    not exhaustive field classifications or new geometric constraints.
    Free coefficients use a shared denominator grid to avoid unrelated
    denominators across a large matrix. Exact affine projection has fill and
    coefficient-height guards; a stopped proposal is not a mathematical
    infeasibility result. Quadratic PSD checks use exact LDL and rational norm
    sign comparisons at `sqrt(d)>0`. A rational point cannot contain a nonzero
    radical component.

17. For SDPA-GMP, use `scripts/run_balanced_dual_gmp_20261126.py prepare`, a
    separate `verify`, then `solve`. Supply the same checked exact dual and cone
    units, balanced-certificate gate and original preparation preflight. The
    exporter reads rational actual matrix entries; it does not reuse binary64
    matrix values. All 100-significant-digit decimal coefficients and the exact
    trace objective are replayed before a native launch. The objective preserves
    the full feasible set and has the explicit strictly positive dual point
    `h I`; see the algebraic derivation in `notes/lp/global_formulation.md`.
    Native dual iterates now have adjoint target `c`, not zero: account for this
    objective before any Farkas recovery. Record requested and reported precision
    bits separately from tolerances and printed decimal places. Preserve native
    decimal arrays without a binary64 conversion. A new tiny input-shift fixture
    checks that the installed solver resolves a `1e-40` perturbation at 256 bits;
    do not rerun unchanged installation tests. Memory/time stops are resource
    outcomes, not solver infeasibility conclusions.

18. Inspect the backend's structural cost as well as coefficient scales.
    `scripts/diagnose_gmp_schur_20261126.py` counts the verified full model's
    block-support pair slots and source-bound integer-cache payload. It does not
    estimate actual total RSS or compute Cholesky fill. The current SDPA-GMP
    representation has 627,025,775 pair slots and 12.54 GB of integer-cache
    payload. `scripts/sdpa_gmp_sparse_graph_20261126.py build` prepares an isolated
    solver with faster assembly of the same symbolic adjacency; compiled
    old/new/independent tests and unchanged numerical-source hashes are required.
    It adds no constraint and changes no matrix coefficient. Preserve the stock
    attempt and identify the private executable and patch in every later run.
19. Native iteration limits write checkpoints; an external wall-time termination
    can leave no arrays. The actual GMP regression is
    `scripts/check_gmp_iteration_checkpoint_20261126.py`. A restart must contain
    exactly the three unlabeled brace arrays extracted from the report, with its
    dates, parameter summaries and iteration table removed. The current full
    attempt did not finish one Newton step within its native batch limit.
    `scripts/prepare_gmp_long_checkpoint_20261126.py plan` binds a proposed
    one-iteration job to the replayed input and private binary. Its longer `run`
    requires an explicitly authorized wall-time policy and matching monitor
    environment; the current 900-second watchdog does not authorize that run.
    The user's subsequent 26 November approval supplies that override for this
    SDP: use `scripts/run_bounded_gmp_hour_10gb_20261126.py`, restricted to the
    approved native-input hash, with one thread, 10 GB process-tree RSS and the
    unchanged percentage gates. The original watchdogs and old receipts remain
    unchanged. The first call requests one native iteration to write a full
    checkpoint; its 3,300-second native time limit expired during the first
    factorization, so it returned no arrays. The new continuation adapter's
    small `check` action passed under the serial watchdog; see
    `results/global_gmp_continuation_checks01/result.json`.
    `continue_gmp_checkpoint_20261126.py run`
    verifies the bound full input and a completed native checkpoint, checks its
    cone-state positivity, and preserves all three native arrays without a
    binary64 conversion. Cone-state positivity alone does not assert affine
    feasibility of the actual model. No full-model continuation was launched
    without a completed checkpoint. Do not repeat an unchanged cold job under
    the same budget when it cannot finish its first factorization.
20. When a full GMP vector is available,
    `scripts/check_gmp_literal_dual_20261126.py` is the prepared direct check of
    every exact G(u) cone, interpreting all printed digits over Q. This adapter
    has not yet been exercised on a full-model output. A failed check rejects
    that vector; it does not decide feasibility. A successful check writes a
    candidate for the separate existing `search_balanced_dual_field_20261126.py
    verify` replay, which must check the full exact adjoint identities as well.

For the smaller explicit dual, final verification checks every cone and all
300,139 adjoint identities directly in exact matrix coordinates. For a recovered
point on the opposite certificate side, verification reinserts the removed dual
rows as zero and checks every prepared-source Farkas coefficient and its negative
constant. Both require the retained reduction proofs. A successful exact check
would still not authorize a geometric bound update before the separate source
review required above.

A time or memory limit never means convergence. Incomplete preparation cannot
authorize a run. Starting arrays need the correct shape and numerically positive
definite cones; native output must actually report all four full-precision formats.
Missing or stale registrations, models, proofs, reconstruction maps or matrix
hashes block launch. Do not invoke the native binary or retired wrappers manually
to bypass the gate.

For a restart, inspect the native residual-ratio thresholds as well as absolute
residuals. In the audited SDPA release, initialization uses `epsilonDash`, while
updates divide by the initial residual above `epsilonDash/100`. A residual
between those thresholds can activate a centering guard and produce tiny steps.
The matched 26 November test changed only `epsilonDash` from `1e-10` to `1e-8`:
the same arrays advanced for 16 iterations instead of stopping after two.
`scripts/run_sdpa_parameter_probe_20261126.py` records the actual parameter diff,
initial-array hash and preflight. Keep this as an explicit diagnosed parameter
choice; do not silently change historical wrappers or relax exact acceptance.

The explicit fixture exception is limited to 16 variables and one MiB of input.
It exists for sign, cone-order and restart regressions, never for research decisions.

## Reading results and retaining evidence

Check candidates in every equality and cone of the frozen original input. Check
certificates by exact coefficient identities and PSD witnesses; reduced-problem
certificates also need the verified reduction chain. Auxiliary feasibility, an
unfinished negative margin, a near-zero dual objective and a small-step error
establish no original decision.

Always identify which problem the solver classified. SDPA updates a coupled
primal/dual pair for the common-margin auxiliary. The explicit normalized
certificate feasibility problem is `A^T y = 0, b^T y = -1, y in K*`.
Its feasibility can certify infeasibility of `Ax + s = b, s in K`; its
infeasibility is not a geometric nonexistence result. After dual facial reduction,
an infeasibility ray may yield only a compressed primal candidate. Check it
exactly, then reconstruct and check the full original model before claiming a
feasible point. `scripts/check_reduced_certificate_result_20261126.py` performs
the compressed check and explicitly leaves the original decision open. Its
optional scalar polishing imposes restrictions on a candidate point only;
failure of those restrictions proves nothing about the necessary-data model.

Distinguish an unsuccessful reconstruction from a refuted solver classification.
`PrimalInfeasible` is a numerical infeasibility conclusion; `NumericalError` is a
separate status. Check the saved binary64 ray in exact source coefficients before
blaming rational rounding. Report scalar violations relative to their row's
absolute term sum as well as global residuals. Increasing the decimal grid for
the same saved binary64 values adds no solver precision. Record actual arithmetic
precision, solver tolerances and positive units in every refinement run.

Keep exact constraint reduction separate from certificate-field selection. A
Q-coefficient parametrization works over any chosen real extension. Recovering
kernel relations over Q or a justified number field can guide a new candidate;
final verification needs exact identities and PSD checks at the specified real
embedding. Such candidate relations are not universal constraints unless proved
independently. No rational or quadratic feasible point is guaranteed merely
because the input coefficients are rational.

Retain scientific inputs, proof chains, current solver input and exact certificates.
Keep compact failure receipts, with watchdog completion separate from mathematical
status. Remove superseded caches using the recovery rules in
`results/repository_search_cleanup_20261125/README.md`.

The gate removes recognized exact degeneracies. It does not prove that the
smallest cone face has been found or that the Schur system is well conditioned.
Later numerical failures require diagnosis while preserving the geometric problem.
