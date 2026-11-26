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
   identities. A new fixed-point audit runs after conditioning. Further exposed
   reductions block export rather than authorizing a partial model.
5. Only a matching `SOLVER_READY` receipt authorizes production SDPA. Use
   `scripts/run_sdpa_checkpoint_20261125.py --native-input READY/input.dat-s
   --initial READY/initial.dat --preflight READY/solver_ready.json --output RUN`.
   Continue later batches with `--resume RUN`, the same input and the same
   preflight. An iteration limit is a checkpoint boundary, not a decision.

A time or memory limit never means convergence. Incomplete preparation cannot
authorize a run. Starting arrays need the correct shape and numerically positive
definite cones; native output must actually report all four full-precision formats.
Missing or stale registrations, models, proofs, reconstruction maps or matrix
hashes block launch. Do not invoke the native binary or retired wrappers manually
to bypass the gate.

The explicit fixture exception is limited to 16 variables and one MiB of input.
It exists for sign, cone-order and restart regressions, never for research decisions.

## Reading results and retaining evidence

Check candidates in every equality and cone of the frozen original input. Check
certificates by exact coefficient identities and PSD witnesses; reduced-problem
certificates also need the verified reduction chain. Auxiliary feasibility, an
unfinished negative margin, a near-zero dual objective and a small-step error
establish no original decision.

Retain scientific inputs, proof chains, current solver input and exact certificates.
Keep compact failure receipts, with watchdog completion separate from mathematical
status. Remove superseded caches using the recovery rules in
`results/repository_search_cleanup_20261125/README.md`.

The gate removes recognized exact degeneracies. It does not prove that the
smallest cone face has been found or that the Schur system is well conditioned.
Later numerical failures require diagnosis while preserving the geometric problem.
