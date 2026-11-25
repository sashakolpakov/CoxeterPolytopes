# Global Coxeter LP/SDP research

The mathematical specification and actual-polytope lift are in
[global_formulation.md](../../notes/lp/global_formulation.md). The retained
[source notes](../../notes/lp/global_sources.md) identify literature inputs.
This is an initial research relaxation; no improved dimension bound or
separate author's acceptance is claimed.

The [execution manifest](experiments.json) records actual runs, their
source and model hashes, failures and remaining work.

The current dimension-57 formulation implements the three requested
extensions in order:

1. Common supporting forests and complete attachment families, with
   symbolic A/B/D ranks and marked positions, shared inverse-Gram
   responses, and component count and rank marginals tied separately to
   finite vertices and classified complete cusps. Spectator components
   remain in those populations.
2. Joint triangle and quadrilateral states. All corners reuse the same
   walls and responses; full corner-kind patterns condition their
   component marginals. An edge component statistic links face sides to
   both endpoints. Strict opposite-side incidence is retained.
3. Sparse polynomial-ideal closure in a shared moment registry:
   available products L(q*x)=0, variance equations L(q^2)=0, and explicit
   kernel equations in the associated PSD cliques.

Its [necessary-data input](../../results/global_shared_cusp_current/necessary_primal.exact.json)
has **36,818 variables, 9,388 equalities, 3,643 PSD blocks and 134,833
scalar inequalities**. The [normalized certificate dual](../../results/global_shared_cusp_current/certificate_dual.exact.json)
finds all feature weights and compatibility multipliers together.
Levels 3, 4 and 5 select the successive extensions; level 5 includes all
three. Every addition has an actual-polytope lift valid uniformly in
dimension. This export uses the established n>=58 result as an anchor.

The [exact checks](../../results/global_shared_cusp_checks/result.json)
verify nested retention, family-deletion formulas, a nonzero-response
local Gram example through all three layers, and rejection of a mixture
that satisfies only the mean of a pointwise identity. The previous
second-level point has no extension preserving all its old moments:
seven new equations fail. These checks do not decide global feasibility.
The current SDP has no verified feasible point or infeasibility
certificate. Before the separate cusp marginals, a
[15-round PSD-cut search](../../results/global_shared_star_search/result.json)
retained violated PSD blocks and supplied no certificate. An earlier
version's SDPA-GMP run exhausted 650 seconds without a checkpoint; its
monitor PASS means worker completion.

The finite truncation still relaxes integer ranks and discrete labels,
uses rational envelopes for H3/H4, and omits many higher moments. It
matches full small-face states but only one shared edge statistic;
complete edge-diagram matching and joint three-face states remain open
extensions. This is a global family formulation, without selected
vertex-forest or marked-face lists.

The retained second-level dimension-57 model has **3,598 moment variables, 450
equalities, 1,592 PSD blocks and 8,996 scalar inequalities**. Its
[necessary-data input](../../results/global_level2_current/necessary_primal.exact.json)
and [certificate dual](../../results/global_level2_current/certificate_dual.exact.json)
retain every baseline constraint and every free feature coefficient.
The new families apply to every existing face or component-count schema:

* All mixed order-three spatial Gram minors, the full negative third
  compound of the Lorentz Gram matrix, a spatial tensor-square PSD
  block, and cosine-weighted spatial PSD localizers.
* Integer falling factorials, total-rank support inequalities, and
  adjacent-integer-gap PSD localizers for all fixed and tail counters.

The [source-bound exact comparison](../../results/global_level2_audit/result.json)
checks retention of all 444 baseline equalities and 8,951 baseline cone
constraints. It proves that the old feasible witness has **no extension**
to the stronger model: some new equations and inequalities involving only
old moments already fail (six equations, 96 scalar inequalities and one
PSD block). The PSD failure is in the full negative Lorentz third compound.
Only seven new moments are needed. The
[specification](../../notes/lp/global_formulation.md) gives uniform
pointwise proofs, before taking moments; the additions are not selected
diagram exclusions. Universal polynomial checks also verify the compound
sign and every spatial mixed minor for arbitrary matrix entries.

**The second-level SDP is also exactly feasible.** A different
[rational point](../../results/global_level2_face_recovery/original_feasible_point.exact.json)
satisfies all 450 equalities, 1,592 PSD blocks and 8,996 scalar
inequalities. The original strengthened input has SHA256
`7f62e85735e146dd660594a09111401754ce0ced24c64033be67ae85b4e3df9a`.
The [independent integer Bareiss check](../../results/global_level2_independent/result.json)
passed, and an [explicit certificate](../../results/global_level2_independent/normalized_dual_infeasibility_certificate.exact.json)
rules out the normalized Farkas dual: all 18,910 adjoint equations pass
exactly. These additions therefore do not lower the dimension ceiling.

The first 35-iteration necessary-data run was inconclusive; the direct
certificate-dual run reached its 650-second numerical limit. A native
restart with an identical reduced input reached `pdOPT` after 12 further
iterations, but its relative errors still allowed a failed original
scalar row. An exact four-term positive identity among that inactive-count
row and three existing bound rows exposes the required face. Rounding
free coordinates to six decimal places and projecting onto the exact
equalities gives the point above. Both original-input PSD checkers pass;
the numerical status alone was not accepted.

The retained first-level dimension-57 model has 3,591 moment variables, 444 equalities,
907 PSD blocks and 8,044 scalar inequalities. Its exact
[necessary-data input](../../results/global_tail57_current/necessary_primal.exact.json)
and [certificate dual](../../results/global_tail57_current/certificate_dual.exact.json)
are retained. It uses the established n>=58 exclusion to cover the
complementary tail; the same compiler also exported both unbounded parity models.

**The frozen first-level SDP is exactly feasible.** A
[rational point](../../results/global_current_exact_recovery/original_feasible_point.exact.json)
satisfies every one of its 444 equalities, 907 PSD blocks and 8,044 scalar
inequalities. The original matrix input has SHA256
`e35594b359e85a02ca706125b5e8259d372c381df3c8ed83122dbcbc3a5b5ebb`.
No geometric constraint was changed. A second
[exact verification](../../results/global_current_independent_feasibility/result.json)
uses integer, symmetric-pivoted Bareiss elimination, independently of the
rational Schur/LDL checker used during point recovery.

Its normalized infeasibility-certificate dual is therefore infeasible.
Its [explicit infeasibility certificate](../../results/global_current_independent_feasibility/normalized_dual_infeasibility_certificate.exact.json)
uses potentials (-x,-1) and PSD cone matrices F_b(x), where x is the
verified primal point. All 15,641 adjoint equalities were checked exactly;
its cone matrices are exactly the independently verified primal slacks.
Here the dual means the Farkas certificate problem with constant pairing
-1. The ordinary zero-objective Lagrange dual also has its zero solution.

The initial LP was numerically feasible. SDPA-GMP first failed on
redundant dual coordinates; exact removal of 40 null coordinates got it
through that factorization, but its ten-minute batch produced no certificate.
A later first-level model then completed 250 rational PSD-separation
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

The first-level feasible witness rules out a separator at that moment
level. The second-level point independently gives the same conclusion
for level 2. Neither witness establishes feasibility of levels 3--5.

Run the current complete model and its dual search, using fresh output
names and the shared watchdog:

```sh
python3 -B scripts/run_bounded_scientific_long_10gb_20261123.py --parent-lock \
  --ledger results/global_shared_star_export_monitor.json -- \
  scripts/run_global_coxeter_formulation.py export --level 5 \
  --threshold 57 --parity odd --upper-dimension 57 \
  --output results/global_shared_star_export

python3 -B scripts/run_bounded_scientific_long_10gb_20261123.py --parent-lock \
  --ledger results/global_shared_star_solver_monitor.json -- \
  scripts/run_global_coxeter_formulation.py solve-current \
  --input results/global_shared_cusp_current/necessary_primal.exact.json \
  --orientation certificate --solver sdpa \
  --output results/global_shared_star_solver
```

Use `search-frozen-psd` with the same `--input` and a fresh output to
search LP approximations strengthened by rational PSD directions. Any
recovered separator is lifted and checked in the complete original SDP.
Failure of that search does not decide its feasibility.

Verify the strengthened point and its dual implication with a fresh
output directory:

```sh
python3 -B scripts/run_bounded_scientific_long_10gb_20261123.py --parent-lock \
  --ledger results/global_level2_recheck_monitor.json -- \
  scripts/run_global_coxeter_formulation.py verify-current \
  --input results/global_level2_current/necessary_primal.exact.json \
  --point results/global_level2_face_recovery/original_feasible_point.exact.json \
  --output results/global_level2_recheck
```

Verify the retained **first-level** point and its implication for the normalized dual
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
With `--solver sdpa`, `--resume-from` instead restarts from the checkpoint's
native vector and both cone arrays; it requires byte-identical ordered
reduced matrix inputs. `audit-strengthening --input
results/global_level2_current/necessary_primal.exact.json` repeats the
constraint-retention and old-witness separation audit, using the pinned
completed baseline verification. Both commands use the same watchdog and
a fresh `--output` directory.

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

Use `--level 2` when exporting the stronger model. The default remains
level 1, preserving the baseline formulation. To solve the retained
second-level input, use `solve-current --input
results/global_level2_current/necessary_primal.exact.json --solver sdpa
--orientation necessary` inside the same watchdog. The certificate
orientation selects the normalized Farkas problem for that same input.

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
