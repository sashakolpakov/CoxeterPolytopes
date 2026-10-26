# Conditional Lean verification

The project uses **Lean 4.32.2 and Std**, pinned by [lean-toolchain](lean-toolchain).
The kernels prove selected internal implications. Geometry, classification,
actual-field interpretation and named computation outputs remain explicit
parameters; printed foundational axioms do not enumerate those parameters.

Read the **[detailed formal guide](https://sashakolpakov.github.io/CoxeterPolytopes/formal_logic.html)**
for module coverage, actual compilation and separate axiom-audit records,
assumption maps, cache provenance and failures. The full predecessor inventory
is retained in [detailed records](../docs/formal_verification_records_20261025.md).

Only when a named replay is needed, use separate unused ledgers:

```sh
python3 -B scripts/run_proof_batch_20261024.py lean-compile \
  CoxeterProofs.CompactDimensionSeventeenComposition --ledger results/local_compile.json
python3 -B scripts/run_proof_batch_20261024.py lean-audit \
  CoxeterProofs.CompactDimensionSeventeenComposition --ledger results/local_audit.json
```

Run from the repository root. Use one monitored process and the shared lock;
Lean has one thread and a 2048 MiB compiler cap. Admission requires 2% free
disk and 15% available RAM; stop below 1% free disk or 6% available RAM. The outer watchdog retains its own
59-second/1 GiB bound. Never run a raw whole-project build or rerun successful
unchanged modules to refresh a date.

Sources, actual receipts and required legacy restoration archives are retained.
Rebuildable compiler products are removed after a work item. A missing target
needs compilation before its separate audit; restoring bytes establishes no
new source-to-object proof. See the [cleanup record](../results/tex_lean_repro_only_cleanup_20261025.json).
