# Lean coverage of the chamber-period closure

20 October 2026. The reusable
[ChamberPeriod.lean](../formal/CoxeterProofs/ChamberPeriod.lean) module
kernel-compiles and passes its seven-theorem axiom audit. It verifies the
group-theoretic final implication used in the completed H9 host
exclusions. It is not an instantiation of the actual isometry groups or
a geometric Lean proof of any host exclusion.

The source SHA256 is
`7bb2d31edcc18d15b5801588890751aa094b7f43d64930ff73144da5b9f25288`.
The [compilation ledger](../formal/.proof-batches/CoxeterProofs_ChamberPeriod.json)
records0.35 seconds; the
[separate axiom-audit ledger](../formal/.proof-batches/CoxeterProofs_ChamberPeriod_axiom_audit.json)
records0.29 seconds. All seven prints pass; there is no sorryAx, and the
only printed logical dependencies are propext and Quot.sound. The module
imports only Std. One induction proof was corrected after a rejected
development compile; both final records bind the corrected source above.

## Exact statement and checked logic

The conclusion constructs a map from Nat into the actual coset type with
pairwise distinct values. This is an infinite-index witness once the
coset type is interpreted as the isometry-group cosets gW of the full
reflection subgroup W.

|Step|Kernel-checked theorem|
|---|---|
|Group left cancellation|`left_cancel`|
|Natural-power addition|`power_add`|
|Strict trace excess forbids a positive period|`no_positive_period`|
|No positive period makes every pair of powers distinct|`powers_distinct`|
|All powers preserve the chamber|`powers_stabilize`|
|Chamber stabilizers have distinct cosets unless equal|`stabilizers_have_distinct_cosets`|
|Combine these into an injection Nat → cosets|`infinitely_many_cosets_from_trace`|

The proof does not assume distinct powers, distinct cosets, infinite
order or infinite index as its final hypotheses. It proves those steps
from the separate inputs below.

## Explicit inputs and remaining obligations

`GroupActionModel` lists the group and action laws. `TraceInterface`
assumes the standard matrix fact that a positive finite order forces
trace at most the ambient matrix dimension. The strict excess for the
specific isometry is a separate premise. No real field, eigenvalue
theory or exact host matrix is silently encoded as a natural number.

`ChamberCosets` assumes that equal gW cosets imply h⁻¹g∈W, and that a
reflection-group element stabilizing its own chamber is the identity.
These are the standard coset and chamber facts. The map g↦gW, the actual
subgroup and its chamber action have not been constructed in this
Std-only project.

The remaining premise is that the specific isometry preserves the full
chamber. In these host proofs this follows from the entire integral-mirror
arrangement and two FULL actual inward corner germs. It is an internally
proved geometric bridge, still supplied as a premise to this module.
The primitive norm inventory and complete edge profiles also remain
written proofs. The exact scripts discharge matrix/finite-arithmetic
interfaces; their outputs are not imported as unproved Lean axioms.

|Host|Written/script input supplying the premises|Still outside this Lean module|
|---|---|---|
|Literal square-root-of-three host|[Completed proof/reproduction](compact_h9_sqrt3_first_edge_reproduce_20261020.md); period trace5+3sqrt3>9|Actual facet transfer, full lattice and group interpretation, complete edges, real trace inequality and chamber-germ bridge|
|Literal golden discriminant59 host|[Completed proof/reproduction](compact_h9_golden_discriminant59_reproduce_20261020.md); period trace40+55phi>10|Actual base facet and full lattice/group interpretation, complete edges, real trace inequality and chamber-germ bridge|
|Literal quadratic13 host|[Far-edge and period proof](../results/compact_h9_sqrt13_a9_far_end_and_period_20261020.md); period trace268eta+80>10|Actual lattice/group interpretation, primitive alphabet, four complete edges, real trace inequality and chamber-germ bridge; the new final worker is separately tracked in its reproduction document|
|Literal discriminant19 host|[Complete proof/reproduction](compact_h9_discriminant19_reproduce_20261020.md); period trace10+14phi>10|Actual lattice/group interpretation, full primitive alphabet and inherited H8 calibration, complete edges, real trace inequality and chamber-germ bridge|

The earlier immutable evidence packets retain their historical
“not Lean verified” status. This later shared module adds conditional
coverage of their internal final logic without modifying either packet
or upgrading the whole geometric result to Lean verified.

## Reproduction and resources

From the repository root:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py \
  CoxeterProofs.ChamberPeriod --audit-axioms
```

The runner holds the shared proof lock, starts only above40GiB disk and
6GiB available RAM, uses one thread and a2048MiB Lean memory cap, and
compiles one module before its separate audit. Both successful starts
had179.26GiB free disk and over8.2GiB available RAM. Temporary build
outputs and the failed development log were removed by the runner;
only the needed `.olean`, source and compact ledgers remain. No other
scientific job ran concurrently.
