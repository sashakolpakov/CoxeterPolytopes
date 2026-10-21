# Lorentz kernel: explicit-transitivity repair review

21 October2026. **Static/source mathematical PASS for the repair; no
successful compilation or actual axiom credit is asserted here.** Current
source `formal/CoxeterProofs/LorentzObstructionKernel.lean` has SHA256
`8a34ae404b1bdb4cfcc8767c88d59a6b57d254fc8cecf0009f7af9775999bf3b`.
Its seven theorem statements, hypotheses and requested axiom commands are
unchanged from the full source review
`results/review_lorentz_obstruction_kernel_source_foundations_20261021.md`
(SHA256 `70a73b90e0c7b8dbc82453888b30e17d5bec9e1f9c5e88730d9d4b533d0a4c2a`).

The first guarded compile exposed an actual tactic failure in the final
grind of second_orbit_gap: its rewritten upper inequality and the
unrewritten Cauchy premise were not connected automatically. The failure
does not invalidate the manual algebra, but it confirms why the earlier
static review explicitly gave no elaboration/axiom credit. That source
and its failed output are preserved, not silently treated as successful.

I read the exact unified delta. It only makes c²>0 explicit by rewriting
c²=c*c and applying OrderedRing.mul_pos, substitutes W=q+c in the
Cauchy premise, and composes that premise with the already rewritten
upper bound by le_trans. Their middle expressions now agree exactly:

    2q² <= (q+c)*((q+c)+second+q) <= 2q²-2c².

The result contradicts 0<c² by linear order cancellation. The generic
transitivity API was checked source-only in the pinned Lean toolchain:
`Lean.Grind.le_trans` in Init/Grind/Order.lean takes exactly these two
inequalities under Std.IsPreorder. No new geometric assumption, cutoff
premise, custom axiom or admitted conclusion is introduced. The other
six proof bodies are unchanged. A bounded recompile and separate actual
axiom audit remain the execution requirement.

The exact predecessor4523-byte source remains at
`results/source_archive/lorentz_obstruction_initial_tactic_failure_20261021/LorentzObstructionKernel.lean`,
SHA256 `8713ee5909934f8b9263ff8cb77c2a0fe7752de603d36edd0168090524964690`.
The preserved actual failed-child and outer outputs respectively match
`7b5b7d16447475a510b40a961d8316a0bede19a4b550be988a9c731139fc7ca8`
and `9319ca8e586befb502265c93039bae13b10d019dc646423d9c29316bdebf65f5`.
The preservation receipt is
`results/lorentz_obstruction_initial_failure_preservation_20261021.json`.
Recovery prints containing sorryAx from failed elaboration have ZERO
verification credit.

The corrected squared-norm map remains SHA256
`9a51f5b4af9cb82b558962aeae79d66ebdbeb00615a8d64ad86adffcc387ffaa`;
its explicit geometric/lattice/rank boundaries are unchanged. I ran no
Lean, compiler, test, mathematical import or scientific job and changed
neither the kernel nor that map in this repair review.
