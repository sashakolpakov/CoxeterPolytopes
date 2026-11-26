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

**Level-7 feasibility remains undecided.** The
[26 November diagnosis](../../results/global_solver_diagnosis_completion_20261126.json)
separates a native restart problem from genuine degeneracy of the certificate
cone. Neither has been resolved by adding or dropping geometric constraints.

SDPA updates both primal and dual variables. Its current input maximizes a common
margin `alpha` subject to `F_b(x) - alpha I >= 0` and `alpha <= 1`. This auxiliary
problem has strict primal starting points even when the original `alpha = 0`
problem is infeasible. A native `pFEAS` label therefore does not decide the
original problem.

The [matched restart test](../../results/global_sdpa_restart_tolerance_matched02/result.json)
used exactly the old input and initial arrays, changing only `epsilonDash` from
`1e-10` to `1e-8`. The old run stopped after two iterations; the matched run
completed its 16-iteration budget without the short-step error. The
[source audit](../../results/global_sdpa_restart_mechanism_20261126.json)
identifies inconsistent initialization and update thresholds for a tiny residual
ratio, followed by the centering guard. Exact acceptance tolerances were unchanged.
The new branch comprises 16 + 16 + 16 iterations, with saved margin about
`-0.000824115`. Its [recomputed diagnostics](../../results/global_sdpa_tolerance_diagnostics02/result.json)
give an ordinary dual adjoint residual of `4.03e-8` and an original, unshifted
prepared-cone minimum of `-0.000823461`. It is not a feasible point. This diagnostic
did not repeat the earlier check of all 29,259 original equalities.

The [exact scalar-row minor](../../results/global_conic_rank02/result.json)
has full column rank 61,788, so ordinary linear redundancy in that coefficient
map does not explain the stop. Nevertheless, exact coefficient checks prove
forced zero rows in the certificate cone. Five facial-reduction rounds verify
274 coordinate directions and remove 767 dual coordinates. The
[separate replay](../../results/global_dual_face_replay01/result.json) checks every
facial implication over the rational source model and every reduced numerical
coefficient. It preserves certificate feasibility; it does not assert that the
minimal face has been found or that the reduced cone has a strict feasible point.

Independent backends were run on the same frozen model:

| Backend and problem | Native result | Verified conclusion |
| --- | --- | --- |
| [SCS, common-margin auxiliary](../../results/global_scs_auxiliary01/result.json) | Inaccurate solution at time limit, 1,075 iterations | Affine residual `0.242`; candidate rejected |
| [Clarabel/QDLDL, unreduced certificate search](../../results/global_clarabel_certificate01/result.json) | `MaxTime`, 8 iterations | No decision |
| [Clarabel/Faer, exactly reduced certificate search](../../results/global_clarabel_reduced_faer01/result.json) | `PrimalInfeasible`, 33 iterations | Numerical classification remains uncertified |

The last row concerns the normalized certificate problem
`A^T y = 0, b^T y = -1, y in K*`, for an original problem `Ax + s = b, s in K`.
A feasible certificate would exclude the original problem. Infeasibility of the
certificate search is a different statement. The returned ray gives a candidate
for a compressed primal problem; its [exact rounding attempts](../../results/global_clarabel_certificate_ray_check01/result.json)
and [scalar polishing attempts](../../results/global_clarabel_certificate_ray_polish02/result.json)
did not produce an exact point. This does not disprove the numerical infeasibility
classification. `PrimalInfeasible` and `NumericalError` are separate native statuses.
Even an exact compressed point
would require a separate lift and check before establishing feasibility of the
full original SDP. No exact feasibility or infeasibility decision, cutoff, or
geometric bound follows from these numerical statuses.

The [earlier original-model check](../../results/global_edge_closure_terminal_check04/result.json)
and [array diagnostic](../../results/global_edge_closure_iterate_diagnostics04/result.json)
remain historical evidence: all 29,259 original equalities held for that older
candidate, but its cones and proposed Farkas completion failed.

The [earlier representation audit](../../results/global_edge_closure_representation_audit01/result.json)
motivated the replacement pipeline. Its old large numerical representations have
been removed; compact execution records remain historical evidence.

## Rational certificate reduction and the saved ray

The [rational reduction record](../../results/global_rational_certificate_completion_20261126.json)
adds an exact reduction of the certificate equalities, after the proved dual-face
reduction. It uses actual symmetric-matrix entries and the trace factor two;
there are no svec square roots in the exact arithmetic. Clearing denominators
and removing integer content gives primitive integer equations.

The [separate replay](../../results/global_rational_certificate_replay01/result.json)
checks all source equalities and the complete reconstruction map:

| Quantity | Before | After |
| --- | ---: | ---: |
| Nonzero certificate equalities | 61,514 | 1 |
| Certificate variables | 361,652 | 300,139 |
| Largest denominator, in bits | 205 | 203 in the reconstruction map |

The remaining normalization has nine terms, including its constant. All scalar
and PSD cones are retained through the exact map. The largest primitive integer
coefficient has 261 bits. Selecting scalar columns that occur in only one
original equation prevents substitution into other equations and bounds map
size by the source coefficient count. The earlier cascading attempt hit its
[15-million-term guard](../../results/global_rational_certificate_prepare01/result.json);
it is retained as an unsuccessful preparation record and authorizes no solve.

The [balanced export](../../results/global_rational_certificate_ready01/result.json)
uses the exact positive unit `62,500,000` for the free variables and the entire
homogeneous cone. Dividing the normalization by `125,000,000,000` gives right
side `-1` and maximum absolute coefficient `1`. All homogeneous cone coefficients
remain unchanged. This is a change of units, with an exact inverse, and preserves
certificate feasibility. Its [Clarabel rerun](../../results/global_balanced_clarabel01/result.json)
ended with `NumericalError` at iteration 45; it obtained no exact decision.

The [saved-ray diagnosis](../../results/global_clarabel_exact_ray_diagnosis01/result.json)
evaluates the binary64 values as exact rationals, before decimal-grid rounding.
The derived compressed point has 2,726 negative scalar rows. Its triangle-cover
value is about `-1.5695e-10`; the worst scalar value is about `-7.1896e-9`, or
7.44% of that row's absolute term sum. Thus the saved ray needs refinement or
structured reconstruction. Its exact interpretation already violates necessary
scalar conditions, while the truth of Clarabel's infeasibility classification
remains undecided. Finer decimal grids of the same binary64 output do not add
solver precision.

For reconstruction, [Dostert–de Laat–Moustrou, Section 3](https://optimization-online.org/wp-content/uploads/2020/01/7553.pdf)
combine exact affine projection with kernel recognition, including quadratic
fields. The rational constraint parametrization remains valid over real number
fields. No field has been identified from this saved ray, and no number-field
certificate has yet been computed. Candidate kernel relations may guide rounding;
they require a separate proof before becoming constraints on every feasible point.

## Explicit dual and field reconstruction

The [completed retry record](../../results/global_dual_field_completion_20261126.json)
records the next actual runs and exact reconstruction attempts. Denote the
balanced certificate search by `C`. Its explicit normalized dual has an exact
[61,513-coordinate parametrization](../../results/global_balanced_dual_prepare01/result.json),
whose [separate replay](../../results/global_balanced_dual_replay01/result.json)
checks every source coefficient and cone. The
[algebraic proof](../../notes/lp/global_formulation.md#dual-certificate-after-exact-equality-reduction)
shows that cone membership in this parametrization gives a normalized dual
certificate proving `C` infeasible. It does not assert full original-primal
feasibility. The original certificate side remains available as well.

| Actual Clarabel/Faer run | Native result | Exact decision |
| --- | --- | --- |
| [Balanced `C`](../../results/global_balanced_clarabel01/result.json), 300,139 variables | `NumericalError`, iteration 45 | None |
| [Explicit dual](../../results/global_explicit_dual_clarabel01/result.json), 61,513 variables | `InsufficientProgress`, iteration 41 | None |
| [Cone-rescaled explicit dual](../../results/global_cone_balanced_dual_clarabel01/result.json), full PSD blocks | `NumericalError`, iteration 54 | None |

All three used binary64. The requested feasibility tolerance was `1e-11` and
the requested infeasibility tolerance was `1e-12`; the native residuals and
returned-array diagnostics remain separate. These runs do not confirm the
earlier numerical `PrimalInfeasible` classification.

The [unit audit](../../results/global_primal_dual_units01/result.json) found
moderate column maxima but whole rows as small as `1.4e-15`. The
[current cone-unit export](../../results/global_dual_cone_units02/result.json)
uses positive dyadic scalar units and PSD diagonal congruences. Four rounds
bring every nonzero scalar or PSD matrix-row envelope, measured in the solver's
symmetric-entry units, into `[0.5,2]`, with
maximum coefficient `2`. All numerical entries have checked inverses and no
coefficient was discarded. This is a coefficient-scale test, not an eigenvalue
bound or a proof of strict feasibility. The final native configuration disabled
chordal decomposition and retained full PSD blocks. It made further progress
before stopping; the simultaneous configuration changes do not isolate a single
cause of the numerical failure.

Exact reconstruction was actually attempted on both sides, over `Q`,
`Q(sqrt(2))`, `Q(sqrt(3))` and `Q(sqrt(5))`. It included literal binary64 dyadics,
exact rational and quadratic unpacking, common grids, LLL kernel proposals and
exact affine projection. Shared grids avoid unrelated rational denominators.
The FLINT arithmetic uses the positive square-root embedding and exact LDL;
its [fixtures](../../results/global_dual_field_checks03/result.json) reconstruct
a nonrational rank-one matrix and reject an indefinite one. A
[separate fixture](../../results/global_source_field_checks02/result.json)
checks nonrational Farkas identities, positive normalization and inverse units.
These fixtures are not certificates for the research model.

The [certificate-side attempts](../../results/global_source_field_attempt01/result.json)
reject all 13 candidates at `finite_vertices_nonnegative_guard_a1`, with value
about `-0.0015467`. The
[explicit-dual attempts](../../results/global_dual_field_attempt01/result.json)
also reject all 12 candidates. The literal dyadic point has triangle-cover value
about `-0.0201117`; stricter kernel projection gives about `-0.0127193`.
Candidate algebraic kernels were detected, but none produced a globally feasible
point. Failed reconstruction does not prove infeasibility in the tested fields,
nor exclude other fields. No full-problem field certificate passed, so no full
candidate-certificate replay or new geometric bound is claimed.

The [local cleanup receipt](../../results/global_dual_field_local_cleanup_20261126.json)
removes only the unused numerical cache of the first, superseded congruence
export. Its execution record, exact unit description and source version remain.
Current exact dual inputs, solver arrays and reconstruction evidence are retained.

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

For the explicit certificate search, encode the same native input with
`scripts/compare_conic_solvers_20261126.py encode`. Then run
`scripts/prepare_dual_face_model_20261126.py prepare` on the exact prepared model,
followed by a separate `verify` invocation. Production certificate solves through
`scripts/compare_conic_solvers_20261126.py solve --mode certificate` require both
`--dual-face-model` and `--dual-face-replay`. Missing proofs block launch; the
explicit `--unreduced-certificate-diagnostic` exception is for comparison only.
The [encoding and gate checks](../../results/global_conic_backend_production_gate_checks03/result.json)
cover signs, both PSD triangle orders, off-diagonal trace multiplicity, known
feasible and infeasible fixtures, and rejection of missing facial proofs.

Continue the explicit certificate path with
`scripts/reduce_certificate_rational_20261126.py prepare` and a separate `verify`,
then `scripts/balance_rational_certificate_20261126.py export`. Each export binds
the prior proof replay and checks its source and artifact hashes. These steps
change the representation of the certificate problem; they add no geometric
restriction and do not alter the frozen necessary-data model.

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
The [later diagnosis record](../../results/global_solver_diagnosis_completion_20261126.json)
adds the matched restart, exact dual-face replay and alternative solver results.
Its [archive receipt](../../results/global_solver_diagnosis_archives_20261126.json)
preserves both new native checkpoints by complete byte round trips. The sparse
encoding, dual embedding and returned alternative-solver arrays are retained
inputs for the current replay and candidate checks.
