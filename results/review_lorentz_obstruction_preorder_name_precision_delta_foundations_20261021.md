# Lorentz kernel: final preorder-name precision review

21 October2026. **Static/source mathematical PASS** at final source
`formal/CoxeterProofs/LorentzObstructionKernel.lean`, SHA256
`5c0f8e57c9a1fb2dd1ca740d2424a5359169cde199c351d782188fe328c0a6d1`
(4765 bytes). No Lean/compiler/test/import job was run by this reviewer.
Actual execution evidence is a separate requirement.

The second guarded compile found name resolution, not a new algebraic
defect: the unqualified le_trans was unavailable. I inspected the complete
Std source namespace this time. Init/Grind/Order.lean line14 opens
**Lean.Grind.Order**, so the earlier review/message's suggested
Lean.Grind.le_trans was misqualified. That static API identification in
`review_lorentz_obstruction_explicit_transitivity_delta_foundations_20261021.md`
is corrected here; it is not evidence of successful elaboration.

The final source instead uses the class field directly:

    Std.IsPreorder.le_trans _ _ _ cauchy upper

Its three points are 2q², `(q+c)*((q+c)+second+q)` and
2q²-2c². It has exactly the required premise/result types and needs
only the already available preorder inherited from Std.IsLinearOrder.
The explicit c²>0 and substituted Cauchy premise from the preceding
repair remain. The resulting contradiction is the same linear order
contradiction checked in the full review. All seven theorem statements,
hypotheses and requested #print commands are unchanged.

I read the exact delta from preserved8a34: besides the explicit
class-field name, only the header changes from pending wording to a
timeless pointer to the logic map, expressly retaining that preparation
alone is not verification. No geometry, lattice, finite-order bridge,
Schur-rank interpretation, axiom or theorem conclusion was added.

The predecessor4690-byte source is preserved at
`results/source_archive/lorentz_obstruction_transitivity_name_failure_20261021/LorentzObstructionKernel.lean`,
SHA256 `8a34ae404b1bdb4cfcc8767c88d59a6b57d254fc8cecf0009f7af9775999bf3b`.
Its failed child output matches
`834f0ef7230611cb76d23fd8772e55cfd7a2cbed615be8fd39a54b8c7bed33ce`;
the preserved outer output matches
`9319ca8e586befb502265c93039bae13b10d019dc646423d9c29316bdebf65f5`.
Receipt:
`results/lorentz_obstruction_transitivity_name_failure_preservation_20261021.json`.
Together with the preserved initial8713 attempt, these failures have
ZERO verification/actual-axiom credit.

The full seven-statement mathematical and assumption review remains
`results/review_lorentz_obstruction_kernel_source_foundations_20261021.md`,
SHA256 `70a73b90e0c7b8dbc82453888b30e17d5bec9e1f9c5e88730d9d4b533d0a4c2a`.
The corrected map remains
`9a51f5b4af9cb82b558962aeae79d66ebdbeb00615a8d64ad86adffcc387ffaa`.
This final source binding corrects the API-name issue without expanding
any mathematical/formal scope. Separate actual terminal compile and
axiom audit results must determine execution status.
