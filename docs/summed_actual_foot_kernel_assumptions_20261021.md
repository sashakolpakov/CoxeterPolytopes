# Summed actual-foot energy: formal scope and application inputs

The [Lean source](../formal/CoxeterProofs/SummedActualFootKernel.lean),
SHA256 `41d183ca03e0e405ce3d86ee838bdc3859c03f23ca19bb2cf9cd5239da050c28`,
has actual bounded compilation and a separate ten-declaration axiom audit
PASS. Compilation took0.76s and the audit0.66s; measured process-tree
peaks were774,914,048B and771,833,856B. All ten prints use only
`propext`, `Classical.choice` and `Quot.sound`, with no `sorryAx` or
warnings. The [verification receipt](../results/summed_actual_foot_kernel_verification_20261021.json)
(`f93f631d`) binds the source, object, three local imported objects,
compiler/configuration, immutable ledgers and exact1229-byte
[actual output](../results/SummedActualFootKernel_axiom_stdout_20261021.txt).
That output has the same hash as the successful compile's child output.
The current logic index counts eight kernels and72 distinct printed
declarations. The seven-kernel publication snapshot and catalogue
selections are unchanged.

## Exact theorem and primitive hypotheses

Work with finite coordinate vectors over a linearly ordered commutative
ring; in the geometric application these are Euclidean real vectors.
Let indices range over `0<=i<n`, coordinates over `0<=c<d`, and use the
ordinary finite coordinate dot product. Set

```text
w = sum_i beta_i omega_i,  E = w.w,  H_i = w.omega_i.
```

The application supplies exactly these hypotheses on the summed range:

* `beta_i>=0` and `u_i>=0`.
* `n_i.n_i=1`, and `n_i.n_j<=0` for distinct indices.
* The reference calibration `w.n_i=beta_i`.
* Each actual-foot pressure `T_i<=x.(omega_i-u_i n_i)`.

Then `summed_actual_foot_energy` proves

```text
||x||^2 >= E + sum_i [2 beta_i(T_i-H_i)
                              +(2u_i-u_i^2)beta_i^2].
```

The calibration and individual pressures are supplied identities and
inequalities; the weighted vector, its shifted norm and the resulting
energy bound are derived. There is no assumed energy conclusion,
shift-square bound, nonadjacency, `u_i<=1`, strict pressure or positive
number of indices. Values outside the finite ranges do not matter.

## What Lean derives

Finite-sum bilinearity gives `E=sum beta_i H_i` and the exact single vector

```text
Z = sum_i beta_i u_i n_i,
W = sum_i beta_i(omega_i-u_i n_i) = w-Z,
sum_i beta_i T_i <= x.W.
```

For ANY nonnegative coefficients `a_i`, the unit and obtuse hypotheses
prove `||sum a_i n_i||^2<=sum a_i^2`: every distinct cross contribution
is nonpositive. Applying this to `a_i=beta_i u_i` gives the shift bound.
Calibration gives `w.Z=sum beta_i^2 u_i`; coordinate bilinearity then
proves `||W||^2=E-2sum beta_i^2 u_i+||Z||^2`.
Finally `||x-W||^2>=0` proves `||x||^2>=2x.W-||W||^2`, and the displayed
conclusion follows. All steps use actual finite coordinate sums.
Adjacent feet are allowed; their nonpositive normal products improve the
lower bound. Individual foot inequalities are summed against the SAME x.

The [independent mathematical source review](../results/review_summed_actual_foot_kernel_source_20261021.md)
(`504c12b7`) reviews the prepared predecessor. Its only subsequent proof
change replaces `exact entries` by `simpa only [sumN_add] using entries`
at the final finite-sum equality. The [independent delta review](../results/review_summed_actual_foot_kernel_sum_linearity_delta_20261021.md)
(`4535f25e`) verifies the archived predecessor, exact one-line difference
and unchanged mathematical scope. The [preserved failed compile](../results/summed_actual_foot_compile_attempt1_failure_20261021.json)
binds that predecessor and exact failure output; the successful receipt
binds the successor. This elaboration repair changes no theorem or premise.

## Whole A20 application and remaining boundary

The accepted written whole-A20 argument uses the twelve indices5..16
of unit A20 normals in21 Euclidean coordinates. Its reference has
`w=sum beta_i omega_i`, `w.n_i=beta_i`, `beta_i>=0`, and
`H_i=w.omega_i`. At a strong index choose `u_i=0,T_i=H_i`.
At an isolated dotted index choose `u_i=gamma_i^-2` and

```text
T_i = sqrt(D_i-u_i)-sqrt(u_i(1-u_i)),  D_i=||omega_i||^2.
```

The actual dotted-pair foot supplies the required non-strict pressure.
This instantiates the structural formula in equation(12) of the
[whole-A20 written proof](../results/arithmetic_whole_a20_summed_actual_feet_exclusion_candidate_20261021.md)
(`457debf9`), supported by its
[full outside-author review](../results/review_arithmetic_whole_a20_summed_actual_feet_outside_author_20261021.md)
and [written acceptance receipt](../results/whole_a20_and_two_end_couplings_written_acceptance_20261021.json)
(`92784b48`).

Still outside Lean are the hyperbolic/farthest-point coordinate realization
and `||x||^2<1`; actual wall incidence, acute feasibility and the dotted
foot construction; the complete strong-or-isolated-dotted classification,
including arithmetic first-type and finite-root arguments; and the
concrete reference realization and nonnegative coefficients. The displayed
rational heights, `E=1.007658`, strict comparisons with radical strong
floors, internal total loss `<.005511`, favorable endpoint contributions
and final numerical contradiction `||x||^2>1.002147>1` are also outside
this module. Its compilation does not instantiate these inputs or prove
the whole-A20 exclusion, an arbitrary-spectator theorem or a dimension bound.

## Reproducing this single-module verification

Run from the repository root, serially with no other proof job. The
unchanged runners require40GiB free disk and6GiB available RAM before
each launch, retain the shared exclusive lock and owner checks, and
abort below8GiB disk or3GiB RAM. Lean uses2048MiB and one thread;
the outer wrapper additionally enforces59s and1GiB process-tree RSS.
Keep its process-monitor access enabled. First inspect the cache/gate:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.SummedActualFootKernel --check-only
```

Use fresh, unused ledger filenames in the two commands below. The first
compiles only a missing/stale target; with the current bound object it
checks cached provenance. If any imported module needs compilation,
obtain a new serial batch admission before proceeding.

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --ledger results/summed_actual_foot_replay_compile_NEW.json -- scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.SummedActualFootKernel
```

Only after that command is terminal/PASS, check the fresh gate and lock
again, then run the separate actual printed-axiom audit:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --ledger results/summed_actual_foot_replay_audit_NEW.json -- scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.SummedActualFootKernel --audit-axioms
```

Preserve new ledgers/raw output before another audit can overwrite the
operational cache records. The runners remove unused intermediates and
retain needed `.olean` imports. Local import objects are hashed; there
is no separately hashed transitive Std object closure beyond the pinned
toolchain/shared-library/configuration bindings.
