# Residual Bellman certificate: assumptions map

The module [ResidualBellmanCertificate.lean](../formal/CoxeterProofs/ResidualBellmanCertificate.lean) addresses the saved-table clause of Proposition `dp-certificate` in [finite_volume57.tex](../manuscripts/finite_volume57.tex). It imports only `Std`. Its SHA-256 is `cd2a05a693cf9dd8c93b5a143bb90ed2ae2932f2ca7d872e5a1139fc05daa3d9` (202 lines, 8274 bytes). **Bounded compilation and all eight actual axiom prints pass.** Each printed declaration uses only `propext` and `Quot.sound`; there is no `sorryAx`. The [verification receipt](../results/publication_logic_verification_20261020.json) binds the source, compiled object, actual output and execution ledgers under Lean 4.32.2. The catalogue now selects this abstract logic kernel for finite57; no concrete scientific table or geometric theorem is thereby formalized.

The finite item alphabet is `Fin n`. A `Certificate n cutoff` supplies positive natural item ranks, arbitrary integer item weights and `memo : Nat → Option Int`, with exactly these local obligations:

- `base`: rank zero stores `some 0`.
- `transition`: every item appended to every stored predecessor within the cutoff has a stored successor whose value dominates predecessor value plus item weight. This quantifies over **all** items and predecessors, not only a selected attaining path.
- `tight`: every stored positive rank within the cutoff has an item and a stored predecessor with the exact rank and weight equations.

`list_upper` uses list induction to derive reachability and the upper bound for every item list within the cutoff. `memo_attained` uses strong rank induction to derive an attaining list for each stored value. Strict rank positivity makes every predecessor smaller. `memo_exact` and `some_iff_exactMaximum` combine these directions; no maximum or attainment conclusion is a certificate field. `none_iff_no_list` derives absence of all item lists. For positive ranks, `none_iff_no_predecessor` derives the equivalent absence of any stored predecessor transition; zero is correctly treated separately as the base.

`totals_perm` proves invariance of both additive totals under list permutation, so the exactness statement also applies to residual multisets. `length_le_rank` and `attained_within_58` give an attaining witness of at most 58 items when the explicit additional input is `cutoff ≤ 58`. Empty alphabets and negative weights are allowed; `none` remains distinct from every integer, including zero.

The application must still supply the complete allowed alphabet, its actual ranks and scaled-64 **integer** residual weights, the context split removing all pair interactions, and a concrete memo table satisfying all local obligations. Those obligations are finite: at most `cutoff + 1` states and `n` item transitions per state, with tight witnesses checked at stored positive states. The module does not implement a table generator or an enumeration checker, instantiate the actual table, prove geometric domination or context completeness, or reinterpret real heights/tensor weights as integers. It therefore supplies the abstract certificate implication and witness bound, not an end-to-end formalization of the finite-volume theorem.

The [parsed actual axiom output](../results/ResidualBellmanCertificate_axiom_summary_20261020.json) covers `list_upper`, `memo_attained`, `memo_exact`, `some_iff_exactMaximum`, `none_iff_no_list`, `none_iff_no_predecessor`, `totals_perm` and `attained_within_58`, with all eight distinct qualified names checked. The saved 802-byte [Lean output](../results/ResidualBellmanCertificate_axiom_stdout_20261020.txt) has the same hash as both the compilation and axiom-audit ledger output. Compilation took0.54s and the separate audit0.41s, each under the one-thread/2048MiB runner. There are no project-local imports, so the local dependency-object map is empty. `Std` is supplied by the pinned toolchain; the runner's compiler/configuration/shared-library fingerprints are not a separately hashed transitive closure of Std objects.

The [source audit](../results/review_residual_bellman_certificate_source_batch_guard_20261020.md)
accepted predecessorc3876078. The [header-change receipt](../results/residual_bellman_publication_logic_preparation_20261020.json)
byte-preserves that source/map and verifies that the entire Lean proof
body is unchanged in currentcd2a05a6. Only the introductory status comment
now refers to the durable compilation receipt instead of saying that
compilation has never occurred. That preparation receipt and the earlier
source-only review remain historical records. The actual compilation and
printed-axiom results are recorded separately in the current verification
receipt; no Lean proof statement changed for this status update.
