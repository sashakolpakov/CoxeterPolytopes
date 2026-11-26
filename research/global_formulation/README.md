# Global Coxeter LP/SDP research

This is the global necessary-data relaxation for a hypothetical finite-volume
Coxeter polytope. The [mathematical specification](../../notes/lp/global_formulation.md)
contains the actual-polytope lift; the [source notes](../../notes/lp/global_sources.md)
record the literature and hypotheses. The accepted 57/17/15 bounds are unchanged.
Independent geometric review remains pending for the research formulation.

The [execution manifest](experiments.json) distinguishes actual execution,
exact verification and numerical solver status. Historical failures are retained
there as compact records, rather than as current search instructions.

## Current scientific model

The frozen **level-7 edge compatibility model** has 97,515 moment variables,
29,259 equalities, 5,322 PSD blocks and 387,685 scalar inequalities. Its
[primal and normalized Farkas dual](../../results/global_edge_closure_current/result.json)
are preserved in [byte-verified archives](../../results/global_edge_closure_current/archive_roundtrip.json).
All 26 manuscript feature coefficients, including contextual and cross-component
terms, remain free in the dual. There are no selected lists of vertex forests or
marked face diagrams.

The occupation measures retain complete supporting-component populations,
attachment profiles and shared inverse-Gram responses. Every triangle or
quadrilateral uses the same supporting walls at all its corners. Common edge
populations link face sides to their two endpoints. Complete rooted spherical
schemas tie these attachments to edge diagrams, while classified complete cusps
supply the ideal vertex data. Sparse moment and polynomial-ideal closure use a
consistent registry for each measure. Unbounded A/B/D ranks and labels remain
symbolic within the proved envelopes.

The model retains every level-5 constraint. Its
[exact comparison](../../results/global_edge_closure_current/strict_strengthening_audit.json)
finds 465 violated equations involving only the previous witness's moments, so
that witness has no extension. The
[local lift checks](../../results/global_edge_schema_lifts02/result.json) cover
nonzero attachment responses and an ideal quadrilateral. These checks establish
neither global feasibility nor geometric nonexistence.

This export treats dimension 57, using the completed n >= 58 theorem as an
anchor. The general formulation can instead represent a dimension band or an
unbounded parity tail. Changing that scope requires a corresponding uniform
instantiation check; it is not an ordinary fixed-scope constraint extension.

**Level-7 feasibility remains undecided.** The corrected, gated run completed
16 + 16 + 2 iterations before SDPA returned
[`step length is too short`](../../results/global_edge_closure_sdpa10/result.json).
The saved common margin is about `-0.00603938`; the recomputed ordinary dual
adjoint residual is `3.17e-6`. These are auxiliary-problem diagnostics.

The [original-model check](../../results/global_edge_closure_terminal_check04/result.json)
found neither an exact feasible point nor an infeasibility certificate. All
29,259 original equalities hold exactly, but the saved point violates
`uniform_two_unit_triangle_cover0`. The attempted exact Farkas completion
requires 34,403 negative coordinate weights. The
[array diagnostic](../../results/global_edge_closure_iterate_diagnostics04/result.json)
records the precise rejection and distinguishes the updated arrays from cached
native summaries. No cutoff or geometric bound changed.

The [earlier representation audit](../../results/global_edge_closure_representation_audit01/result.json)
motivated the replacement pipeline. Its old large numerical representations have
been removed; compact execution records remain historical evidence.

## Constraint and solver workflow

Follow [CONSTRAINT_WORKFLOW.md](CONSTRAINT_WORKFLOW.md). It specifies the
uniform family argument, parent-constraint retention, source-bound checks,
exact reduction to a fixed point, proof replay, numerical units and mandatory
native launch receipt. Aggregate incidence equations must not be treated as
pointwise identities when multiplying by state coordinates.

The [completed regeneration](../../results/global_edge_closure_regenerated01/result.json)
exactly matches the frozen scientific input. The
[fixed-point preparation](../../results/global_edge_closure_presolve06/result.json)
has 61,787 variables and 289,596 cones. Its
[post-reduction audit](../../results/global_edge_closure_presolve06/fixed_point_audit.json)
reports zero remaining recognized reductions. Each face-propagation round and
every deleted or compressed source cone has an exact reconstruction proof.
The [complete-chain replay](../../results/global_edge_closure_replay02/result.json)
passed. The [scaled export](../../results/global_edge_closure_ready03/result.json)
has a matching [launch receipt](../../results/global_edge_closure_ready03/solver_ready.json)
and a clean audit after conditioning. Scaling checks both coordinate inverses
and every affine matrix entry. It uses sparse elimination and preserves shared
population coordinates when a response-coordinate pivot is available.

The production path consists of these separate serial jobs:

1. `scripts/regenerate_global_constraints_20261125.py`
2. `scripts/prepare_conic_model_20261125.py`, repeated using `--previous` until
   the recorded fixed-point audit passes
3. `scripts/verify_conic_reduction_chain_20261125.py`
4. `scripts/export_ready_global_sdp_20261125.py`
5. `scripts/run_sdpa_checkpoint_20261125.py`, with explicit native input,
   checked initialization and matching `--preflight`
6. `scripts/recover_sdpa_checkpoint_20261125.py`, with the new prepared source,
   native checkpoint and frozen original input explicitly supplied

Every scientific command runs under the shared lock and the monitored runner in
[AGENTS.md](../../AGENTS.md), with fresh output and ledger names. For example:

```sh
python3 -B scripts/run_bounded_scientific_long_10gb_20261123.py --parent-lock \
  --ledger results/research_prepare_pass1_monitor.json -- \
  scripts/prepare_conic_model_20261125.py \
  --original results/global_edge_closure_current/necessary_primal.exact.json \
  --output results/research_prepare_pass1
```

Restore the frozen JSON from its retained `.gz` archive if needed. A completed
pass can still be `NOT_SOLVER_READY`; repeat it until no recognized reductions
remain. Neither a timeout nor a strictly feasible auxiliary starting point
substitutes for this test. Original candidates must pass every original equation
and cone. Certificates require exact identities, PSD witnesses and the verified
reduction chain. The gate does not promise a well-conditioned Schur system or
that all possible cone faces have been found.

The old exploratory solve wrappers are retired. Their mathematical helpers
remain where current reduction or recovery needs them. Sign, cone-order,
full-precision output and restart checks use explicitly marked tiny fixtures.

## Retained exact obstruction to weaker formulations

The full **level-5 SDP is exactly feasible**. Its
[rational point](../../results/global_shared_cusp_scalar_circuits01/original_feasible_point.exact.json)
satisfies all 9,388 original equalities, 3,643 PSD blocks and 134,833 scalar
inequalities. An [independent integer Bareiss check](../../results/global_shared_cusp_independent/result.json)
also verifies the certificate ruling out its normalized Farkas certificate dual,
including all 193,440 adjoint equations. This witness proves that the weaker
relaxation cannot exclude dimension 57; it asserts no geometric realization.
These accepted research inputs and certificates remain intact.

## Cleanup and restoration

The [25 November cleanup record](../../results/repository_search_cleanup_20261125/README.md)
lists the removed obsolete searches, notes and numerical payloads. Superseded
tracked payloads restore from the pushed tag `pre-search-cleanup-20261125` at
`70b23827`. Untracked obsolete generated searches were deliberately deleted;
the tag does not back them up. Compact historical status and source records are
preserved with their original scope. Current scientific inputs, exact witnesses,
certificates and proof-chain evidence take priority over disposable numerical
caches.

The [archive receipt](../../results/global_constraint_pipeline_archive_20261126.json)
binds every reduction proof, the final exact models/maps, numerical input,
initialization and all three new native checkpoints by complete byte round trips.
Large raw files are ignored local working copies. Restore their `.gz` archives
before a replay or continuation in a fresh checkout; the regenerated input's
[restoration note](../../results/global_edge_closure_regenerated01/README.md)
uses the already-retained frozen archive. Intermediate compiled models are
reconstructed from the original and reduction proofs, so their obsolete raw
copies are not required for a fresh replay. Restoring bytes does not rerun any
mathematical verification.

The [completion record](../../results/global_constraint_pipeline_completion_20261126.json)
collects the verified preparation, terminal solver status and cleanup receipts.
