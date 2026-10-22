# Root retry gate: explicit unfolding in arithmetic16 frontier proofs

22 October2026. Root has read the complete predecessor failure and the
complete successor diffs, and the ENTIRE113-line different-author repair
review. **Corrected SOURCE accepted for a capped retry; no successful Lean
compile or actual axiom-audit credit is granted by this document.**

| Object | SHA256 |
|---|---|
|formal/CoxeterProofs/ArithmeticSixteenFrontierRefinement.lean|f00fa9541ee0960096ddfdcb33a5729ae775089b5474c2937d4a221be78a9114|
|docs/arithmetic_sixteen_frontier_refinement_logic_20261022.md|739fa74848eae7c3911c2e17f20a4b4db370997613b0c8048fce92ff9a684312|
|results/review_arithmetic_sixteen_frontier_explicit_unfolding_repair_foundations_20261022.md|f5fa9c78be247d8ee7e80a232d592a44bd5b44ccb2a66871ba71d15a23017646|
|results/source_archive/arithmetic16_frontier_before_survives_decide_repair_20261022/manifest.json|a1b4806a48a913e27df5c6d746ba4ec3b26219d7c8da70d28d0be173b854f6b4|
|results/lean_arithmetic16_frontier_thirteen_group_verification_20261022/first_compile_failure.txt|515c50f5fb4e779502fc5d32d17f6aab0fd9e166a0d2e3a0ef31386b6e618a43|

The first ACTUAL2048MiB/one-thread compile exited1. Its retained `decide`
branches did not unfold `Survives` for decidability; typed branches also
left a free meaning parameter. All six printed entries in that failed
output, including three downstream `sorryAx` entries, earn ZERO accepted
compile/audit credit. No successful object was published. The previous
static gates did not discover this implementation issue and are historical.

The exact source diff is the explanatory header and sixteen retained
membership proofs replaced by

    simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]

This expands the conjunction, eliminates the unused typed parameter through
the group map and simplifies explicit constructor membership. Every datum,
theorem type, binder, excluded-case contradiction, nine written interfaces,
same-point witness and six print requests is unchanged. It therefore fixes
the diagnosed tactic problem without changing the mathematical conclusion:
nineteen original groups refine to thirteen, while A14 retains four cases.

The source/map and both old static gates were byte-preserved before editing.
Root created and checked those four gzip round trips; the different-author
review independently rechecked their hashes/lengths and decoded hashes.
The map now records the actual failure and binds the corrected source.
All geometry and arithmetic assumptions remain external, as in the complete
source gate6839. No global15 or geometric instance is inferred.

Proceed only through the unchanged signal-safe2048MiB/one-thread runner,
one module then a SEPARATE actual axiom audit, each with fresh40GiB/6GiB
admission and8/3 abort floors. Successful older modules are not rerun.
