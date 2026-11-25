# Global Coxeter LP/SDP research

The mathematical specification and actual-polytope lift are in
[global_formulation.md](../../notes/lp/global_formulation.md). The retained
[source notes](../../notes/lp/global_sources.md) identify literature inputs.
This is an initial research relaxation; no improved dimension bound or
separate author's acceptance is claimed.

The [execution manifest](experiments.json) records the actual first runs,
their source and model hashes, failures and remaining work. The current
dimension-57 anchored model has 3,591 moment variables, 444 equalities,
907 PSD blocks and 8,044 scalar inequalities. Its exact
[necessary-data input](../../results/global_tail57_current/necessary_primal.exact.json)
and [certificate dual](../../results/global_tail57_current/certificate_dual.exact.json)
are retained. It uses the established n>=58 exclusion to cover the
complementary tail; the same compiler also exported both unbounded parity models.

**The frozen current SDP is exactly feasible.** A
[rational point](../../results/global_current_exact_recovery/original_feasible_point.exact.json)
satisfies every one of its 444 equalities, 907 PSD blocks and 8,044 scalar
inequalities. The original matrix input has SHA256
`e35594b359e85a02ca706125b5e8259d372c381df3c8ed83122dbcbc3a5b5ebb`.
No geometric constraint was changed. A second
[exact verification](../../results/global_current_independent_feasibility/result.json)
uses integer, symmetric-pivoted Bareiss elimination, independently of the
rational Schur/LDL checker used during point recovery.

The normalized infeasibility-certificate dual is therefore infeasible.
Its [explicit infeasibility certificate](../../results/global_current_independent_feasibility/normalized_dual_infeasibility_certificate.exact.json)
uses potentials (-x,-1) and PSD cone matrices F_b(x), where x is the
verified primal point. All 15,641 adjoint equalities were checked exactly;
its cone matrices are exactly the independently verified primal slacks.
Here the dual means the Farkas certificate problem with constant pairing
-1. The ordinary zero-objective Lagrange dual also has its zero solution.

The initial LP was numerically feasible. SDPA-GMP first failed on
redundant dual coordinates; exact removal of 40 null coordinates got it
through that factorization, but its ten-minute batch produced no certificate.
The stronger model then completed 250 rational PSD-separation
rounds without a verified separator or a verified full-SDP feasible point.
Those runs were inconclusive. Solving the necessary-data orientation with
exact cone-face reductions gave SDPA-GMP phase `pFEAS` in 35 iterations
(334 seconds, 128-bit arithmetic). A further exact four-term nonnegative
identity exposed three remaining redundant coordinates. Rounding the
free coordinates to six decimal places and projecting onto the exact
equalities then gave the rational witness above. Both PSD checkers passed.
The accepted geometric dimension bound remains 57.

The organizing objects are occupation measures for complete vertex
populations, all triangles and quadrilaterals, residual corners, and
global incidence counts. A single shared Lorentz Gram schema covers
every small face. Complete component families have symbolic rank tails.
All26 intrinsic manuscript features remain common marginals, with free
dual coefficients. This code contains no selected list of vertex forests
or marked face diagrams and performs no diagram-family pricing.

The initial sparse moment level is deliberately explicit. It can be
strengthened by raising moment/localizer degrees or adding proved
uniform identities. The LP version replaces PSD blocks by declared
necessary rational-direction tests. Its failure to find a certificate
does not resolve the SDP or the geometric problem.

The present feasible witness rules out a separator at this moment level.
Further progress requires stronger valid constraints or a stronger
moment level; longer solver runs on this unchanged model cannot exclude it.

Verify the retained point and its implication for the normalized dual
without running a numerical solver:

```sh
python3 -B scripts/run_bounded_scientific_long_10gb_20261123.py --parent-lock \
  --ledger results/global_feasibility_recheck_monitor.json -- \
  scripts/run_global_coxeter_formulation.py verify-current \
  --output results/global_feasibility_recheck
```

Use fresh output names. This check needs Python's standard library and
reads the original exact matrices. `solve-current --solver sdpa
--orientation necessary` solves that frozen input; `--orientation
certificate` selects its normalized certificate problem. `--solver recover
--resume-from results/global_current_sdpa_primal` reconstructs a rational
point from the saved numerical checkpoint.

From the repository root, run one monitored batch at a time:

```sh
python3 -B scripts/run_bounded_scientific_long_10gb_20261123.py --parent-lock \
  --ledger results/global_tail57_export_monitor.json -- \
  scripts/run_global_coxeter_formulation.py export \
  --threshold 57 --parity odd --upper-dimension 57 \
  --output results/global_tail57_export
```

Use fresh output names. `--upper-dimension57` is written as
`--upper-dimension 57` on the command line. It uses the completed
n>=58 theorem as an anchor. Thus one global model for the band N<=n<=57
suffices for a proposed cutoff N; separate diagram families or separate
runs for each dimension are not needed. Omitting the option keeps the
whole unbounded parity tail in the new model itself.

After the source/cleanup commit is pushed, `explore-lp` exports both the
necessary-data problem and its certificate dual and searches the latter.
`explore-sdpa` sends the complete conic certificate problem to SDPA after
exact equality elimination. Both retain numerical status separately from
exact rational identities and PSD checks. Every accepted result must be
checked in its original complete matrix; support reduction is only a
candidate-recovery device.

`explore-psd-cuts` adds rational quadratic tests from violated original
PSD blocks. It lifts any recovered LP certificate back to the original
SDP and checks it there. This adds constraints to the same global
matrices, without growing a list of diagram types. The model includes
the uniform old two-unit triangle cover, contextual rank-support rows,
and positive-coordinate support localizers. Old quadrilateral transfers
whose premises require the dimension-58 LP remain outside this first model.

The angle-square observable couples inverse-Gram sums over all finite
vertex components to the face cosines. Its A/B/D formulas and rational
unbounded-tail envelopes are proved in the specification. Its dual
coefficient is free along with the 26 old feature coefficients.

The generic compiler includes affine equalities in its Farkas alternative
and includes the factor two in off-diagonal trace pairings. It does not
fix corner coefficients, impose coefficient boxes or silently replace a
partial vertex/face list by an exhaustive one. Results store source hashes
instead of duplicating the whole source tree at every run.
