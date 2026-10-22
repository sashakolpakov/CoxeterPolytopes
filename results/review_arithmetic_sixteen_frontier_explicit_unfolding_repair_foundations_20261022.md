# Static delta review: explicit frontier-membership unfolding after failed compile

22 October2026. **Corrected-source/delta/assumption-map PASS; no correction
requested. No successful compilation or axiom-audit credit.** No Lean,
compiler, mathematical import, worker, test, enumeration or TeX job ran
by this reviewer. Only source and preserved metadata were read.

Reviewed ENTIRE repaired source/map and their complete predecessor diffs:

* `formal/CoxeterProofs/ArithmeticSixteenFrontierRefinement.lean`,
  SHA256 **f00fa9541ee0960096ddfdcb33a5729ae775089b5474c2937d4a221be78a9114**,
  12,157 bytes/286 lines;
* `docs/arithmetic_sixteen_frontier_refinement_logic_20261022.md`,
  SHA256 **739fa74848eae7c3911c2e17f20a4b4db370997613b0c8048fce92ff9a684312**,
  14,454 bytes/213 lines.

This supplements my full predecessor review
`results/review_arithmetic_sixteen_frontier_thirteen_group_precompile_full_scope_foundations_20261022.md`,
SHA256 **cdc6b2c0414b74e1642cb4ea33bea5977213eefa403e60737ed05868170476c6**.
I did not author the formal source or this repair; I authored some external
written geometry premises, as disclosed there. This reviews the logical
composition and exact repair, not an independent self-review of those
geometric premises or an assertion that Lean has accepted the successor.

## Actual predecessor failure, preserved without retrospective credit

The actual full compiler output is
`results/lean_arithmetic16_frontier_thirteen_group_verification_20261022/first_compile_failure.txt`,
SHA256 **515c50f5fb4e779502fc5d32d17f6aab0fd9e166a0d2e3a0ef31386b6e618a43**,
4,663 bytes. I read all of it and freshly matched its hash/length.
The preserved parent stdout also matches9319ca8e586befb502265c93039bae13b10d019dc646423d9c29316bdebf65f5
and records the child-exit1 exception. The archive manifest records the
2048MiB/one-thread compile as FAILED, with no successful object.

The errors concern exactly the retained branches of
`classified_forest_refinement`: `decide` did not unfold `Survives` enough
to synthesize its Decidable instance; the two retained typed branches
also retained a free meaning parameter. Here "opaque" means opaque to
that tactic's decidability synthesis, not a new Lean `opaque` declaration.
My earlier static review missed this real elaboration issue and therefore
did not establish compilability. Its source/assumption analysis remains
historical; the failed run is NOT retrospectively relabeled PASS.

The failure stdout prints all six requested names, but the downstream
three include `sorryAx`. Even the first three printed entries are NOT
credited as a successful module compile or separate actual axiom audit.
The whole failed batch earns zero accepted outputs. Those error-recovery
placeholders are diagnostic evidence, not source proofs or accepted axioms.

## Exact code change and why it addresses the diagnosed failure

The ENTIRE source diff has only the explanatory header and SIXTEEN
retained branch replacements:

    decide
      --> simp [Survives, ForestCase.group,
                survivingGroups, survivingA14Spectators]

Those are the twelve retained non-A14 branches (including both typed
branches) and four retained A14 spectator branches. No inductive type,
constructor, list, definition, theorem statement, binder, written premise,
excluded-case contradiction, selection witness, or print command changed.

In each replacement, unfolding `Survives` exposes the conjunction;
`ForestCase.group` removes the irrelevant typed parameter, and the lists
reduce the group membership to explicit constructor equalities. For A14,
the match additionally exposes precisely its four-label membership. The
lists contain each retained constructor, so simplification proves the
goal without asking for a Decidable instance on the unexpanded predicate.
This directly addresses BOTH recorded diagnostics. The explicit A14 list
argument is unnecessary in some branches but does not weaken the goal or
introduce a mathematical premise; no zero-warning result is predicted.

The closed original-enumeration `decide` proofs remain unchanged. The
separate derived ALL-THREE-case D14 theorem and all external contradiction
applications are unchanged. The final existential still returns exactly
the same `family, forest` produced for the same supplied p. Thus the
nine external interfaces,19-to13 refinement, A14 five-to-four refinement
and actual typed meanings are exactly those of the complete prior review.
No `Survives` conclusion or whole-D14 conclusion has been added as input.

## Complete preservation and truthful current map

The failure/preservation manifest is
`results/source_archive/arithmetic16_frontier_before_survives_decide_repair_20261022/manifest.json`,
SHA256 **a1b4806a48a913e27df5c6d746ba4ec3b26219d7c8da70d28d0be173b854f6b4**.
Its current hash and all four gzip hashes/byte lengths match. Independent
decoding recovers the complete named predecessors with exact SHA256:

| Object | Decoded SHA256 |
|---|---|
| source11dd |11dd035a5556b7a8627889b98dd9ba6f835844503302007b05a57c03c9714e33|
| mapf41b |f41bf48414e500400fee889475ed637fea186e9c305eb3c848c30c4b61c6d328|
| full outside static reviewcdc6 |cdc6b2c0414b74e1642cb4ea33bea5977213eefa403e60737ed05868170476c6|
| root static gate6839 |68391402af0c111ceef7e061332f7bcc2d5df795ce5f02860600e921c53a43f5|

The map diff updates the exact source binding/status, records the actual
failure and its full stdout, explains both diagnostics and this repair,
and changes the proof description to explicit unfolding. ALL nine
mathematical bindings, current written root gates, application meanings
and the earlier18cf/815e archive chain are unchanged. Old static gates
remain source-bound historical records and are not promoted to this source.

The current map correctly distinguishes one ACTUAL FAILED compile from
six REQUESTED prints and zero ACCEPTED actual outputs. It claims no
successful current object or axiom list. All geometry, field arithmetic,
classification, energy and actual-vertex instantiation remain external;
the conditional refinement proves no new global bound or new example.

**Verdict:** the exact repair is logically sound and addresses the actual
elaboration errors without changing mathematics or assumptions. It is
ready for root's source gate and separately admitted capped retry/audit.
Only such actual successful evidence can close the remaining formal gate.
