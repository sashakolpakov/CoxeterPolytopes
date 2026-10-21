# Anchor matching and paid-foot logic

The source is [AnchorAndPaidFootKernel.lean](../formal/CoxeterProofs/AnchorAndPaidFootKernel.lean).
Its six declarations have a [source-logic review](../results/review_anchor_and_paid_foot_kernel_source_20261021.md).
All six now have **successful bounded compilation and actual axiom audits**:
[verification receipt](../results/anchor_paid_foot_kernel_verification_20261021.json),
[compile record](../results/AnchorAndPaidFootKernel_compile_20261021.json),
[audit record](../results/AnchorAndPaidFootKernel_axiom_audit_20261021.json),
and [exact679-byte output](../results/AnchorAndPaidFootKernel_axiom_stdout_20261021.txt).
The first two declarations are axiom-free; the other four use only
`propext`, `Classical.choice` and `Quot.sound`, with no `sorryAx`.
This is conditional logic verification, not a geometric dimension theorem.

A first
bounded compiler invocation found a reserved binder name and a failed
rewrite of a partially applied definition. The
[failed source/output record](../results/anchor_paid_foot_initial_elaboration_failure_20261021.json)
is preserved. Renaming the binder and exposing the definition before
rewriting changes no statement or assumption. The repaired source's
retry was refused before a child launched at5.4GiB available RAM.
The6GiB launch gate was not lowered. The later successful compile and
audit started with6.55/6.58GiB available RAM and156.42GiB disk, in one
serial one-thread2048MiB batch. They took2.15s/.44s. These ledgers record
global resources and the compiler cap, not measured peak process-tree RSS.
The predecessor map is byte-preserved in
[this archive](../results/source_archive/anchor_paid_foot_logic_before_successful_verification_20261021.md.gz).

| Declaration | What the proof derives | Application inputs still required |
| --- | --- | --- |
| `same_k_reflection` | For the concrete19 indices, equal i(19-i) means equal or reflected indices | None beyond the finite arithmetic; kernel `decide` is used, not `native_decide` |
| `weak_index_unique` | Necessary matches plus the reflected obstruction leave one weak index | Actual cap extraction, verified matching output and written Cartan-integrality obstruction |
| `remaining_weight_pressure` | Summing nonnegative retained weights preserves the original reference floors | Actual floors and the chosen retained weights |
| `paid_reference_pressure` | The paid old-coefficient subtraction and one actual foot give the new pressure | Pointwise payment identity, old reference energy/height and feasible-foot pressure |
| `paid_reference_energy` | One finite-coordinate completed square gives the full graft energy | Actual payment, flat-tail and squared-tail identities; SAME original vector throughout |
| `paid_reference_contradiction` | The derived energy above1 contradicts squared norm below1 | The concrete strict scalar comparison and interior norm bound |

For A18, the weak predicate on `Fin 19` must be FALSE outside the central
indices5..14. The census supplies no statement about weak caps at the
other raw height indices. The module does not silently extend that domain.

The script propositions are necessary matching statements, never claims
of field isomorphism, a field-square or geometric realization. The Cartan
obstruction, actual field reduction, complete finite-white inventory,
feasibility of all auxiliary points, concrete reference construction and
numerical bounds remain written inputs. No parameter assumes the final
dimension bound, and no new geometric formalization is claimed.

For a future reproduction with the unchanged40GiB disk/6GiB RAM gate,
the existing signal-safe runner reuses unchanged capped dependency
objects and audits the six actual printer commands:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.AnchorAndPaidFootKernel --audit-axioms
```

The receipt above records the actual successful run; the command is a
reproduction instruction. It retains the exclusive lock, one-thread2048MiB
cap and cleanup. The owned temporary work directory was removed.
Printed axiom names alone would not replace the application assumptions
in this map. The earlier ten selected kernels/87 declarations retain
their scope. This module raises the total to eleven selected kernels
and93 distinct audited declarations; none is an end-to-end geometric proof.
