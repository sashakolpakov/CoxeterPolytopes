# Exact arithmetic for the large inert parent financing proof

**Bounded compilation and all four actual axiom outputs PASS.**
The [compile snapshot](../results/large_inert_parent_compile_20261021.json)
records0.55s with6.03GiB available memory; the [separate audit](../results/large_inert_parent_axiom_audit_20261021.json)
records0.41s with6.58GiB. Both use one thread and a2048MiB cap.
The retained [375-byte actual output](../results/large_inert_parent_actual_axioms_20261021.txt)
matches both ledger hashes. Its [complete parsed summary](../results/large_inert_parent_axiom_summary_20261021.json)
has exactly the four requested declarations, all axiom-free. The ledgers
record the memory cap and admission resources, not a measured peak RSS.
No existing successful kernel was rerun and no catalogue selection changed.

Source: [LargeInertParentArithmetic.lean](../formal/CoxeterProofs/LargeInertParentArithmetic.lean),
SHA-256 `a62fe70788b10ee6b1f47d046ab732b79b3ad2ccc5c77b069efc53d42f919c7d`
(72 lines / 2618 bytes). Its only import is `Std`; it imports no local
geometry, scientific worker, saved table or earlier arithmetic module.

The target is Section 5 of the written
[parent26–40 all-port financing proof](../results/finite57_large_inert_parent_26_to40_all_port_financing_20261021.md),
SHA-256 `6679d964f8e6f5e99cdf54209928bca6fdb260c4af1edabbc9deb7fb1928bf2f`.
The complete proof was read before preparing the source. This map records
the arithmetic interface; it does not replace that proof's independent
review or assert its acceptance.

Put `M=40-p` in the integer range `0,...,14`. The definitions compute

\[
\begin{aligned}
\binom n4&=0\quad(n<4),\qquad
\binom n4=n(n-1)(n-2)(n-3)/24\quad(n\ge4),\\
A_M&=30464-\binom{30-M}4+\binom{13-M}4-\binom{M-10}4,\\
V_M&=17(41-M)(M+1),\\
\Gamma_{2,M}&=(84-M)V_M-227(41-M)(M+1)\\
&\qquad-34(256-5M)(M+1)-10A_M,\\
D_{2,M}&=\begin{cases}
64(41-M)(M+1),&M\le4,\\
4(41-M)(M+1),&M\ge5,
\end{cases}\\
R_{2,M}&=\Gamma_{2,M}-D_{2,M}.
\end{aligned}
\]

Here `linkVertexCount` computes `V_M`, the link-vertex count, not
hyperbolic volume. Also `Gamma2_M=2*Gamma_M` and `R2_M=2*R_M` in the
written application. The split at M4/M5 is exactly the change from p36 to p35,
with source debit eight versus one half per eligible face. The last
binomial term is retained: it contributes one at M14. At M0 the same
formula includes the rank-zero factor convention.

| Kernel-checked theorem | Exact content |
|---|---|
| `coefficient_table` | All fifteen `A_M`, in increasing M, computed from the binomial expression |
| `canonical_surplus2_table` | All fifteen doubled canonical surpluses from the displayed signed formula |
| `remaining2_table` | All fifteen doubled remaining surpluses after the complete piecewise debt |
| `remaining2_lower_bound` | For every `M : Fin 15`, `173 <= R2_M` |

Each proof is `by decide`, requesting kernel reduction of a closed
finite proposition; no `native_decide`, added axiom, `sorry` or supplied
table-correctness premise occurs. Each theorem has a `#print axioms`
command. The three `List.range 15` equalities compute the formulas rather
than merely repeat an assumed list. The lower bound covers all fifteen
indices, not a sampled minimum. All numeric variables are exact
combinatorial integers; real heights, norms and curvature have not been
replaced by integers.

Still outside this module are the coefficient-extraction interpretation,
the link h-polynomial and curvature identity, canonical geometric charge
envelopes, convex-order concentration and its direction, the comparison
between the concentrated and actual link-vertex counts, positive-margin
scaling, complete finite forest classifications, the cheap/high source
split, actual face incidence and injection, all-port counts, once-only
charge reassignment and compatibility with older allocations. In
particular the module does not prove the selected-family financing
theorem, a global finite maximum, ceiling56 or any polytope existence.

The existing bounded runner reproduces this module in two separate batches:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.LargeInertParentArithmetic
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.LargeInertParentArithmetic --audit-axioms
```

The linked immutable receipts record both actual executions. The second
requires terminal compile PASS and a fresh admission. Keep the unchanged shared
lock, 40GiB disk/6GiB available-memory start gate, 2048MiB cap, one thread
and process monitor; root owns execution and subsequent immutable output
bindings. The six paper sources/PDFs, existing kernels and catalogues
are unchanged. The prepared map144b67fd is preserved byte-for-byte at
`results/source_archive/large_inert_parent_arithmetic_map_before_actual_compilation_20261021.md.gz`.
This execution adds four exact arithmetic facts, not a geometric Lean theorem.
