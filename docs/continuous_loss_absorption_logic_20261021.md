# Continuous loss absorption: written proof and Lean boundary

21 October 2026. The new
[shorter-A proof](../results/arithmetic_a15_a17_improved_loss_summed_feet_20261021.md)
has full outside-author written review. It proves component energies
above .8348, .8703 and .9032 for arithmetic A15, A16 and A17 with
arbitrary finite spectators. It removes the ordinary-T2 spectator
beside A16 and the ordinary nonshort A1 beside A17. No global ceiling17
follows.

The analytic step treats a continuous unknown parameter, rather than
checking finitely many sample values. For an isolated or foreign-cleared
cap, its actual foot has pressure

    x.(omega_i-u_i n_i)>=sqrt(D_i-u_i)-sqrt(u_i(1-u_i)),
    0<u_i<1.

Given a positive reference w=sum beta_i omega_i with w.n_i=beta_i,
sum these inequalities BEFORE estimating the energy. The shifted vector
W=w-sum beta_i u_i n_i has its full overlap with w retained.
Nonpositive products between different old normals bound its square,
even for adjacent feet. It suffices to absorb the loss by

    2[a u+sqrt(u(1-u))-ell] < beta(2u-u^2).

The written proof supplies the three rational parameter rows:

| Component | a upper bound | beta positive lower bound | ell | r upper bound for sqrt(1+a²) |
| --- | ---: | ---: | ---: | ---: |
| A15 | .1982 | .022 | .601 | 1.0195 |
| A16 | .1955 | .0205 | .6 | 1.019 |
| A17 | .194 | .018 | .6 | 1.019 |

## What the new kernel proves

[ContinuousLossAbsorptionKernel.lean](../formal/CoxeterProofs/ContinuousLossAbsorptionKernel.lean)
works over an arbitrary linearly ordered commutative ring. It uses no
division, rational proxy for real vectors, real-analysis library or
assumed dimension conclusion. Its seven declarations prove:

1. The circle identity from s²=u-u².
2. The exact two-square identity for the loss a*u+s.
3. Its continuous upper bound 2(a*u+s)<=a+r from r>=0 and r²>=1+a².
4. The bound2s<=1 for the nonnegative square root.
5. Adjusted loss absorption by splitting2u<=1 and2u>=1. The abstract
   coefficient tests are a+1<2ell and4(a+r-2ell)<3b, with beta>=b>=0.
6. Nonnegativity of the corresponding weighted energy correction.
7. The final reference-energy inequality from actual-foot comparisons,
   primitive unit/obtuse normals, nonnegative weights and duality.

The last theorem uses the previously kernel-checked
[summed actual-foot theorem](../formal/CoxeterProofs/SummedActualFootKernel.lean).
The loss maximum and interval absorption are derived here, not supplied
as conclusion-shaped fields. The pointwise scalar absorption parameter
of the final theorem can be discharged using the preceding theorem.

Zero weights need no scalar loss estimate. The final theorem still
accepts a foot comparison at every index. At an omitted zero-weight
slot, set u=0 and loss=w.omega+ell-x.omega; the foot comparison is
then equality and its energy contribution is zero. This dummy completion
does not assert a geometric floor at the omitted index. Unit/obtuse/dual
conditions still concern the actual complete old-root system.

## What remains written or externally supplied

The kernel does not identify its scalar type with a concrete real
coordinate model. An application supplies that realization and the
nonnegative square-root value s. The actual all-wall feasibility and
farthest comparison, complete changed-threshold cap inventory,
first-type field/conjugate arguments, golden time-table packing,
reference coefficients, rational parameter comparisons, and endpoint
H4/graft construction remain the written arguments. The concrete
component bounds and ambient dimension conclusions are not Lean
theorems merely because these generic implications compile.

## Actual execution and preserved failed attempt

The first capped compilation failed because the Std-only environment
does not include the tactic `by_contra`. Its original source and actual
compiler diagnostics are preserved. Only the two tactic blocks were
replaced by `by_cases` with the same polynomial contradiction; statements
and hypotheses did not change. Failed output containing `sorryAx` is
not credited as verification.

The repaired source compiled in .82 seconds and its separate actual
axiom audit ran in .89 seconds. Four unchanged dependencies were reused
from capped objects; there was one new module and no aggregate build.
Both invocations used one thread and the 2048MiB compiler cap, with
start gates40GiB disk/6GiB RAM and abort gates8GiB/3GiB. The ledgers
record global resources, not a measured process-tree peak RSS.
All seven actual axiom lists contain exactly propext, Classical.choice
and Quot.sound, with no sorryAx. One unused-section-variable warning
remains; no cosmetic rerun was made. Temporary products were removed
by the existing signal-safe runner; its necessary olean is retained.

The source, compile and actual-output bindings are in the
[verification receipt](../results/continuous_loss_absorption_logic_verification_20261021.json).
To inspect the single-module replay without running it:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py \
  CoxeterProofs.ContinuousLossAbsorptionKernel --check-only
```

For a necessary fresh replay, replace `--check-only` with
`--audit-axioms`; retain the same resource gates and serial scheduling.
