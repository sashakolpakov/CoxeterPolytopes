# A12 adjusted-loss kernel: preserved two-step tactic repair

23 October 2026. **SOURCE-ONLY DELTA PASS.** I independently compared the
entire current source with the decoded predecessor and checked both repaired
proof obligations. This is clearance of the exact repair, not a successful
Lean compilation or axiom audit. I did not execute Lean, import a worker,
run tests, or change a reviewed source.

## Exact bindings and preservation

- Current `formal/CoxeterProofs/A12AdjustedLossKernel.lean`:
  SHA256 `f02dededf06fa12f53811d1fa4be5da982747ee42b78629050290916c4ac6952`,
  10,219 bytes / 242 lines.
- Decoded predecessor: SHA256
  `5ab0165dd843e9599f7a2c4ccb68701f2bcc9f5d701d3715ba5a7eecf9f4858d`,
  10,075 bytes / 237 lines. Its gzip is
  `results/source_archive/a12_adjusted_loss_before_explicit_inner_scaling_and_empty_case_20261023/A12AdjustedLossKernel.lean.gz`,
  SHA256 `5cb90e721c623d1a4dd1fd2f103cf3e079ede2ebc663e6ea437d0577fd105237`,
  2,804 bytes. I independently decoded it and recovered these exact
  predecessor bytes, matching the preservation entry in
  `results/lean_a12_adjusted_loss_verification_20261023/attempt_2_failed_compile.json`.
- Actual failed-attempt payload:
  `results/lean_a12_adjusted_loss_verification_20261023/compile_attempt_2_lean_failure.txt`,
  SHA256 `1b9ea49427d3bed97c30574e7d3ae30a0321d038d2c0885dfe539ffcea6bd7af`,
  4,298 bytes. Its two errors are the inner-role final `grind` and the
  empty-list nonnegativity goal. The payload also prints `sorryAx` in
  failed descendants; those prints are failed-compilation evidence and
  confer no successful axiom-audit credit.

## Complete delta and algebra

The byte comparison has exactly two hunks, both inside existing proofs.
No import, declaration, statement, hypothesis, structure, constant,
requested print, or other proof body changes.

1. In the inner role of `individual_adjusted_loss`, the source adds
   `scaled_beta : 2000 * beta * (2*u-u^2) = 61 * (2*u-u^2)`, proved by
   `rw [inner.1]`, before the unchanged final `grind`. This is precisely
   multiplication of the existing equality `2000*beta=61` by the indicated
   ring element; it needs no sign or extra geometric premise. Subtracting
   it from the existing inequality
   `4000*loss < 820*u + 4000*b` gives

       4000*loss - 2000*beta*(2*u-u^2)
         < 698*u + 61*u^2 + 4000*b < 2402.

   The desired conclusion follows by positive integer scaling by two.
   The added equality supplies the nonlinear product missing from the
   old tactic's recorded arithmetic context.
2. The `finite_gain` nil case changes only from
   `simp [sumCredit, sumDebit]` to unfolding these two definitions with
   `simp only`, followed by `constructor <;> grind`. Both sums are zero.
   The first conjunct is `0 <= 0-0`; the implication in the second has
   the contradictory premise `[] != []`. These are the correct empty
   base obligations, with no positivity assumed for an empty list.

The repaired arguments are valid in the same explicitly assumed ordered
commutative-ring setting. Their actual elaboration and the absence of
`sorryAx` in a new successful object remain for the root's separate guarded
compile and axiom audit. Attempt 2 remains FAILED; this review does not
change that history or assert an attempt-3 result.

I authored/reviewed some of the inherited geometric supporting arguments,
but did not author this Lean implementation or this two-hunk repair. This
report concerns the implementation delta and its exact existing interfaces,
not an independent reproof of the external geometry or a new dimension
bound.
