# Compact eighteen: prepared conditional composition kernel

**CORRECTED SOURCE; NOT YET SUCCESSFULLY COMPILED.** The first admitted
compile failed because Lean omitted section hypotheses used only in the
theorem bodies. The successor adds the explicit `include written` before
the four conditional theorems. Its23 premises, five conclusions and proofs
are otherwise byte-identical. The failed stdout and Lean diagnostics are
[preserved](../results/lean_compact18_composition_verification_20261022/attempt1_failure_and_preservation.json),
including the exact old source/map in byte-verified gzip snapshots. The
failed elaboration's placeholder `sorryAx` prints are not an accepted
axiom audit. Full successor source/map review and separate admitted compile
and audit are still required. Existing successful modules and six papers
are unchanged.

The new [Lean source](../formal/CoxeterProofs/CompactDimensionEighteenComposition.lean)
is **10051 bytes / 224 lines**, SHA-256
`441bc21675cb6e2ff707c2267468ac667f06984009d9d118b9ab2d6c0d341287`.
It imports **Std only**, with no local module imports, added axioms,
`sorry`, native evaluation, real/number-field library or geometric model.
Its target is the complete written
[global composition](../results/compact_dimension18_global_composition_and_whole_d19_exclusion_candidate_20261022.md),
SHA-256 `9a8558607afaac0889b1d359a361ec1ca42942e2dc4b725f320bf126bb382ed7`
(15462 bytes). Source hashes identify written premises; they do not
instantiate them as Lean proofs. This preparation makes no acceptance,
execution or global-bound promotion beyond its conditional logic.

## What the source derives

All public theorem names have namespace
`CoxeterProofs.CompactDimensionEighteenComposition`.

| Declaration / future axiom request | Derived conclusion |
|---|---|
| `weak_end_count_cases` | Two Boolean weakness flags give exactly zero, one or two weak sites |
| `whole_a19_exclusion` | The zero/one/two written exclusions exhaust whole A19, at the same vertex |
| `whole_b19_exclusion` | The actual two-terminal reduction and its two exclusions exhaust whole B19 |
| `whole_d19_exclusion` | The actual head reduction and all six branch exclusions exhaust whole D19 |
| `compact_dimension_at_most_eighteen` | Prior GLOBAL19 plus the complete exact-dimension19 frontier and its exclusions imply dimension at most18 |

The dependency graph is acyclic: the Boolean lemma feeds A19; terminal
and head witnesses feed B19 and D19; these three conclusions feed the
nine-group forest dispatch; only the final theorem uses prior GLOBAL19.
Under the contrary `dimension > 18`, that prior bound gives the exact
equality19 by integer order reasoning. The frontier then supplies ONE
actual farthest vertex and ONE forest witness. No downward inheritance
of existence, Coxeter-face assertion or fresh optimizer is invoked.
There is no input `dimension <= 18`, whole dimension19 nonexistence or
whole-A19/B19/D19 impossibility. The three whole exclusions are derived.

## Typed meanings and quantifiers

The types are finite **labels for external mathematical interfaces**.
Their constructors do not prove that geometric forests exhaust these
labels. That theorem is the explicit `F19` premise. `Forest19` has nine
grouped constructors: A16, A17, A18, whole A19, B18+actualI, whole B19,
D17+typed7, D18+actualI, whole D19. A16 retains all four spectator
families, including both typed-seven meanings; A17 retains ordinary T2
and either typed7; A18 retains ordinary nonshort A1 and actualI.

`TypedSeven` separately names OLD A2 / NEW isolated7 and OLD I2(7) /
NEW isolated3. Ordinary T2 means ordinary irreducible I2(m), m>=3; it
is not assigned to a typed defect. ActualI is the frontier's ACTUAL
isolated short-A1 opposite, not an arbitrary old A1. The precise native
opposites, ordinary/nonshort convention, fixed original vector, root
normalizations and spectator incidences remain part of interpretation
and the written proofs. No arithmetic or field-list hypothesis is added.

The two Boolean functions mean weakness at the SAME indices5 and15,
with the SAME nonstrict .61 convention. Their geometric definition is
external; the kernel proves their Boolean count is exhaustive. The
`A19one` premise quantifies over either orientation, rather than assuming
symmetry of the candidate. Equality belongs to weak, as in all three
written sources.

`B19Terminal.shortThree/shortFour` mean exactly the two actual surviving
last-long alternatives of B19r (own3 and the respective short contact).
The kernel does not derive that reduction from a diagram or skip its
actual-wall hypothesis. `D19Head` denotes only the post-DH isolated
native5/7/8/9, finite-label>=10 and strictly dotted classes. Its finite
label carries `10 <= label`; dotted represents the ENTIRE actual
continuous magnitude interval greater than1. Excluding own0/3/4/6,
parallel and nonisolated cases is DH's EXTERNAL geometric theorem.

All predicates concern ONE fixed compact candidate and its dimension.
The witness p from F19 is passed unchanged through every branch. The
written input fields and the final theorem quantify over arbitrary p,
labels and typed meanings. They do not choose different points to add
energies. No uniqueness of a fork cap or equality of fork parameters is
assumed. There is no assertion that declaring a predicate gives a metric
or Coxeter instance.

## Exactly23 written input fields

`WrittenInputs` is a Prop-valued record with EXACTLY the following23
fields, in the composition's literal source-table order. All23 current
source hashes were streamed and matched while preparing this map. Paths
below are under `results/`; the actual source paths, hashes and scopes
are also in the linked composition. Only source bytes were inspected;
no transitive mathematical replay is claimed.

| Field | Written source | SHA-256 |
|---|---|---|
| G19 | `compact_dimension19_and_arithmetic_dimension17_global_composition_20261021.md` | `0ac802108e0ef8a9ca821c3091fca82f7bf9331a0037639d93888f2fabb37a75` |
| F19 | `compact_dimension19_complete_forest_frontier_candidate_20261021.md` | `c443832caaf78411d85e718b8d2813e1bba76d9b3f066ae388bc5dc6fa107f70` |
| A16 | `compact19_a_family_continuation_candidate_20261021.md` | `a707755e7a6e7f3172f524bf1cf378fb4a2d934b30265a12f81968ba8ec94a43` |
| A17o | `compact19_a17_ordinary_rank_two_endpoint_foot_exclusion_20261021.md` | `4faf7073434c7b9bdd9815ac80c4718f169dd7c51f9c38e2a556e63527f1e431` |
| A17t | `compact19_a17_typed_seven_endpoint_coupling_candidate_20261021.md` | `bc1b2acdcddd885a0bf350505ac3c8aa3d39216cbf61e47c151943a0e4cc35ea` |
| A18o | `compact19_a18_ninety_one_hundredths_endpoint_foot_candidate_20261021.md` | `c093ab0189dd8198bf842c66de98689eab8885761d728169c4edf6be0d9cbe01` |
| A18I | `compact19_a18_isolated_i_two_end_paid_foot_exclusion_candidate_20261021.md` | `b040a9b8e488154d94beb589cb0e3c11ca79a0eb34366a233b1ee803728fd820` |
| B18I | `compact19_b18_isolated_i_full_endpoint_foot_exclusion_candidate_20261021.md` | `d7e6fa43104c001c0e0263e9b3e6e98d5fdc00c56783fbab394383a7f0ef9cdb` |
| D17t | `compact19_d17_typed_seven_fresh_centre_exclusion_candidate_20261021.md` | `bef4e450bf91109579888661c2d9385df14afbe3434bf2c63cc36575fd7a62dc` |
| D18I | `compact19_d18_isolated_i_complete_centre_and_head_exclusion_candidate_20261021.md` | `8b8e4b90f4059a85778fe0a9c033d80c10d02a932d877b0d993dc588973bed73` |
| B19r | `compact_whole_b19_two_terminal_wall_reduction_foundations_20261022.md` | `f62d79fdbcf76b9eaa90d14e90b929487bc87de8dacdd681ea7157a07fba9664` |
| B19four | `compact_whole_b19_short_four_terminal_fresh_index17_exclusion_foundations_20261022.md` | `4c7e28403ba97c0f58a64a8b4daff79289ab03c66941937df601dfb7e3701ad8` |
| B19three | `compact_whole_b19_short_three_terminal_fresh_index17_exclusion_foundations_20261022.md` | `a960457719cf3da6ebdf55db36ac5aa69c8232279dec2deb633f4aafe234f76a` |
| A19two | `compact_whole_a19_two_weak_ends_exact_cross_paid_feet_candidate_20261022.md` | `2f6e2083921476c89d9fc082457b282dc12d22cc514497e8a6b095eb4ce44c81` |
| A19one | `compact_whole_a19_one_weak_end_joint_isolated_endpoint_candidate_20261022.md` | `1437bf2fd6f873fff7b17979dc4199e87aa5261fe51ebf1d8ce72a10ea17b445` |
| A19zero | `compact_whole_a19_zero_weak_ends_endpoint_foot_exclusion_candidate_20261022.md` | `9a0f93aee72cd788d5b6e28e1117cb4a7388b0f2df6acb6fea5f385b71d75333` |
| DH | `compact_whole_d19_long_head_complete_reduction_foundations_20261022.md` | `d65a3daa56ea3a5dfcc72e95e9c23ec9051caa7c2f655cbd511bcf9e20bb863c` |
| D5 | `compact_whole_d19_native_five_head_adjacent_vertex_exclusion_foundations_20261022.md` | `ee9a9a3f9cc8b79cea6286da18efaa53ce5f3712ba4d1919ff82a34e4ce3ce3c` |
| D7 | `compact_whole_d19_native_seven_head_complete_exclusion_foundations_20261022.md` | `b0f20a6942ae853c966d41475d0c4a6c4a368101f9144972fec0689138e88bf4` |
| D8 | `compact_whole_d19_native_eight_head_complete_exclusion_foundations_20261022.md` | `abcfffae1c07659c60f86ace8154db3b8a7dcd430fccceb78589bb7400472931` |
| D9 | `compact_whole_d19_native_nine_head_complete_exclusion_foundations_20261022.md` | `30c055cc50079e14fa3ed1cdc353a636554939bae45f51ba9ad7bf7f0f8cb00e` |
| Dlarge | `compact_whole_d19_finite_long_head_at_least_ten_exclusion_foundations_20261022.md` | `3f14f2b7071e14e808d1facb0b8ff34fee3855ce7263366c739712372205bc8b` |
| DD | `compact_whole_d19_uniform_dotted_head_complete_exclusion_foundations_20261022.md` | `27e28380f347c0481f9e49d97027fb15572efd99b125dae926437a7eecdf9c60` |

G19 supplies a prior dimension bound valid for ALL higher dimensions;
F19 supplies existence and COMPLETE nine-group coverage at EXACT19.
A16 handles its entire typed spectator type. A17o/t and A18o/I remain
separate; B18I, D17t and D18I are their complete actual-spectator
exclusions. B19r and DH return actual witnesses for the SAME p. The
remaining fields exclude exactly those terminal/head/count alternatives.
A18o incorporates the written .91 plus ordinary>.1 contradiction at
its exact interface; the kernel does not independently prove either
real inequality. No numeric forest enumeration is represented as
geometric classification by a finite `decide` computation.

## Geometry and verification left external

Compactness, acute simplicity, existence of a farthest vertex from a
fixed strict interior point, finite/affine/native diagram classification,
actual opposite and reciprocal endpoints, distinct-wall native gaps,
all-wall future feasibility, positive finite-Gram inversion and the
original farthest-vector comparison are not formalized here. Neither
are the component energies, exact ranks and typed defect propagation,
complete changed-threshold inventories, continuous radical bounds,
omission/chord losses, one-vector payments, signed fork checks or
adjacent-q slack and monotonicity. The earlier global tail and its
published/script/companion inputs remain external written mathematics.

The author manually reviewed every new declaration, dependency and
quantifier against the complete written composition. This is an AUTHOR
source check, not independent review or elaboration. The expected proof
is elementary logic plus `omega`/Boolean `decide`; no claim is made
about actual compiled axioms before the separate compile/audit occurs.
Actual source/object/compiler/standard-library bindings and five real
printed outputs require future execution receipts. Even a future PASS
would check this conditional composition only, not instantiate its23
geometric premises or formally certify a new polytope/bound application.
