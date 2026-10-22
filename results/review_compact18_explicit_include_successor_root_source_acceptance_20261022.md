# Root source acceptance: compact eighteen explicit-include repair

22 October2026. **Corrected source accepted for a fresh separately admitted
compile, followed by a separate actual axiom audit. No successful compile
is supplied by this source acceptance.**

Root reread the entire224-line successor
`formal/CoxeterProofs/CompactDimensionEighteenComposition.lean`, SHA-256
`441bc21675cb6e2ff707c2267468ac667f06984009d9d118b9ab2d6c0d341287`,
10,051bytes, the updated entire map
`docs/compact_dimension_eighteen_composition_logic_20261022.md`, SHA-256
`587a918c72432c176a3e1773c188a180f3834c126205d43c0843edbb3b80fb29`,
11,537bytes, and the full128-line other-author report
`results/review_compact_dimension_eighteen_composition_explicit_include_successor_source_20261022.md`,
SHA-256 `26e39eb99dd0adbc4d4e647342784d1ba9bc0808ad4f1a83b7245c1448d585ed`.

The actual first attempt failed elaboration. Root preserved its diagnostics,
wrapper stdout and byte-verified old source/map snapshots in
`results/lean_compact18_composition_verification_20261022/attempt1_failure_and_preservation.json`.
The printed sorryAx placeholders belong to that failed attempt and give
zero accepted axiom-audit credit. Both root's original source check and
the earlier other-author review missed the omitted section context.

Root made the sole17-byte Lean edit: `include written` and a blank line
between the record variable and the four conditional theorems. This
explicitly includes the shared written record and its dependent parameters
which were previously referenced only in theorem bodies. All23 input
fields, five conclusions, proof bodies, imports, types, names and print
requests are unchanged. The independent review confirms the same full
logic, actual-vertex quantifiers, native/typed alternatives and external
geometry boundary. Root accepts that report at source-review scope.

There is still no geometric instance and no input dimension<=18 or whole
A19/B19/D19 impossibility. The intended implication remains conditional
on the23 named written theorems. A future successful compiler/audit run
must use this new hash; the predecessor's acceptance cannot be relabelled
as execution evidence. All successful older modules remain untouched.

The next launch remains one module at2048MiB/one thread under the existing
shared lock, with40GiB disk/6GiB available-RAM admission and8GiB/3GiB abort
floors. A fresh failure or resource refusal must be retained. Compilation
and the five actual printed-axiom outputs must be separate jobs.
