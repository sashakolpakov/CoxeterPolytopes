# Lean verification of the internal proof logic

This is a real Lean4 project, pinned by `lean-toolchain` to the locally
installed Lean4.32.2. It currently uses only Lean's standard library, so
no unpinned mathematical library or network download is required.

For a new named module, run separate compilation and actual audit commands
from the repository root. Choose unused ledger filenames:

```sh
python3 -B scripts/run_proof_batch_20261024.py lean-compile \
  CoxeterProofs.CompactDimensionSeventeenComposition --ledger results/local_compact17_compile.json
python3 -B scripts/run_proof_batch_20261024.py lean-audit \
  CoxeterProofs.CompactDimensionSeventeenComposition --ledger results/local_compact17_audit.json
```

The public entry monitors one owned process tree under the exclusive lock.
Lean uses one module, one thread and a2048MiB compiler cap; the outer
research watchdog retains59seconds and1GiB RSS. The current percentage
policy admits work at2% free disk and15% available physical RAM and stops
below1% disk or6% RAM. Exact measured byte ratios determine the gates.
Do not use a whole-project build or an uncapped direct Lean command. Reuse
successful unchanged modules and dependencies.

The read-only `--check-only` command reports cache provenance. Existing
objects without a guarded-build ledger are labelled
`LEGACY_CACHED_UNCAPPED`; their reuse does not claim a new capped verification.
New objects have streaming hashes of sources, imported local objects,
compiler/configuration and the resulting object recorded in
`formal/.proof-batches/`. Changed inputs prevent publication; a mismatched
ledger prevents cache reuse. Output retention is bounded to 512KiB.

The runner removes temporary Lean/C/IR/LSP and TeX outputs after every
batch, including failures and entirely cached runs. Objects needed by the
next separate audit remain during that sequence. After a completed work
item, [reproduction-safe cleanup](../results/tex_lean_repro_only_cleanup_20261025.json)
removes regenerable objects and preserves exact legacy objects with their
trace/hash provenance in a local byte-verified archive. Its manifest records
original paths, hashes and nanosecond timestamps. Sources, certificates,
compilation/audit evidence, final PDFs and pinned configuration remain.
A fresh checkout can give sources newer timestamps: the runner must then
refuse stale cache reuse. Do not change times to force acceptance or infer
a source-to-legacy-object link from restoration. Only a necessary new replay
should rebuild a missing target, followed by its separate actual audit. Below2% free disk it
writes `.proof-batches/cleanup_needed_percentage_20261025.json`: the coordinating agent must
start a cleanup agent before further jobs. The subprocess cannot spawn a
Codex agent itself. Frozen evidence and active application data must be
preserved. Mathematical work takes priority; Lean is final verification
of settled statements.

The aggregate check also uses the capped runner for its axiom audit:

```sh
python3 scripts/check_proof_translation_20261018.py --lean-only
```

The root module prints the axioms used by its audited theorems. The current
audited build uses only Lean's standard `propext`, `Quot.sound` and
`Classical.choice`. No custom axiom, `sorry`, `admit`, `native_decide`, or unsafe
proof shortcut has been introduced. Standard foundational axioms are not
the same thing as a published external mathematical assumption.

## Current coverage

`QuadraticDyadicTripleKernel` has [bounded compilation and a SEPARATE
fourteen-declaration actual axiom audit PASS](../results/lean_quadratic_dyadic_triple_verification_20261023/manifest.json).
Compile2.35seconds; audit.50seconds reused the object. All14 actual prints use
only propext and Quot.sound, with no sorryAx. Its [input map](../docs/quadratic_dyadic_triple_kernel_logic_20261023.md)
keeps its immutable preparation status/hash; the execution receipt supplies
the later PASS. Primitive Nat valuation/parity constraints yield the selected
dyadic and bounded-norm deductions. Actual fields, maximal-order strips,
K-square exceptions, geometry and the prior global theorem remain external.
One unused time_parity warning is preserved. This is not the full geometric
arithmetic15 formalization and supplies no actual number-field instance.

`B16IsolatedPairScalarKernel` has [repaired-source compilation and a SEPARATE
twelve-declaration actual axiom audit PASS](../results/lean_b16_isolated_pair_scalar_verification_20261023/manifest.json).
Compile2.09seconds; audit1.66seconds reused the object. The first print uses
propext/Quot.sound; the other11 also use Classical.choice. No successful print
uses sorryAx. Both literal scalar-prefix contradictions and the cubic root
bound compile. The [logic map](../docs/b16_isolated_pair_scalar_kernel_logic_20261023.md)
retains actual ideal/field facts, original geometry, census, primitive pressures
and real instances as external premises.

The first admitted B16 compile failed on reserved local names. Its complete
source/map/diagnostics are [byte-preserved with ZERO audit credit](../results/lean_b16_isolated_pair_scalar_verification_20261023/failed_compile_source_preservation.json).
The exact two-name [outside delta review](../results/review_b16_isolated_pair_scalar_reserved_identifier_and_map_delta_compact_20261023.md)
and [root clearance](../results/review_b16_isolated_pair_scalar_reserved_identifier_repair_root_clearance_20261023.md)
preceded the successful retry. Two unused-section-variable warnings remain.
Neither new audit recompiled an already successful module; both imported only
Std. The unchanged40GiB/6GiB admission,2048MiB/one-thread named runner and
runtime monitoring applied. No process-tree RSS peak is inferred where the
Lean receipt records none. Local objects remain cache provenance, never
publication assets; older selected-module scopes are unchanged.

For a fresh checkout, named compilation and audit are separate commands:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.QuadraticDyadicTripleKernel
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.QuadraticDyadicTripleKernel --audit-axioms
```

`B18TerminalPressureKernel` has full different-author source/map review,
[root precompile clearance](../results/review_b18_terminal_pressure_kernel_root_precompile_acceptance_20261023.md)
and an [actual bounded compile and separate ten-declaration axiom audit PASS](../results/lean_b18_terminal_pressure_verification_20261023/manifest.json).
Compilation took1.71seconds and the separate audit1.16seconds. The audit
reused the object; no old successful module ran. All ten actual prints
contain only propext, Classical.choice and Quot.sound, with no sorryAx.
Two unused-variable warnings are preserved; owned temporary work is absent.

The [logic map](../docs/b18_terminal_pressure_kernel_logic_20261023.md)
exposes actual geometry, inventories, future feasible points, original
farthest-vector pressure, terminal sums and the real instance as ordinary
explicit premises. The code derives the continuous unbounded polynomial
seal and subcritical radical/convex-sum contradiction. It does not verify
the full geometric near3/own6 theorem, wholeB18 or later role exclusions.
Its source/map preparation labels remain immutable history; actual PASS
comes from the execution receipt. The runner reports the2048MiB Lean cap,
one thread and disk/RAM monitoring; no process-tree RSS peak is recorded
or inferred for this package.

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.B18TerminalPressureKernel
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.B18TerminalPressureKernel --audit-axioms
```

`A13IEndpointGraftScalarKernel` has full different-author source/map review
and [root precompile acceptance](../results/review_a13_i_endpoint_graft_scalar_kernel_root_precompile_acceptance_20261023.md),
and an [actual compile PASS in1.33seconds](../results/lean_a13_i_endpoint_graft_verification_20261023/pre_audit_manifest.json),
and now has a [separate actual eleven-declaration audit PASS](../results/lean_a13_i_endpoint_graft_verification_20261023/manifest.json).
The audit took1.28seconds at6.34GiB RAM; no successful module was recompiled.
Its earlier audit admission was refused at5.3GiB RAM. Its first compile admission was refused
at4.7GiB RAM. The next admitted compile [failed at two tactic scalings](../results/lean_a13_i_endpoint_graft_verification_20261023/attempt_2_failed_compile.json).
Two explicit identity rewrites have [outside/root clearance](../results/review_a13_endpoint_graft_two_identity_rewrites_root_preclear_20261023.md);
the repaired retry was [refused before Lean at5.8GiB](../results/lean_a13_i_endpoint_graft_verification_20261023/attempt_3_resource_deferral.json).
The failed diagnostic prints, including error-propagated sorryAx, have zero
verification credit. The successful compile admitted6.49GiB RAM and retained
its682968-byte ignored object. All eleven actual audit prints use only
`propext`, `Classical.choice` and `Quot.sound`; none uses `sorryAx`.
The [logic map](../docs/a13_i_endpoint_graft_scalar_kernel_logic_20261023.md)
separates the continuous square/factor gains and case dispatch from actual
geometry, inventories, individually legal prices and weak-foot absorption.
Those strong geometric vector interfaces are ordinary explicit arguments.
This is conditional scalar verification, not a full geometric group9 theorem.

`A12AdjustedLossKernel` has [bounded compilation and a separate
ten-declaration actual audit PASS](../results/lean_a12_adjusted_loss_verification_20261023/manifest.json).
It proves continuous individual-loss and finite-list paid-square logic from
primitive ordered-ring data. Compilation took2.37seconds and the separate
audit1.51seconds. All ten prints use only `propext`, `Classical.choice` and
`Quot.sound`; no `sorryAx` occurs. Its first admission was refused below6GiB
RAM. Its second admitted compile failed on two tactics;
the [failure evidence](../results/lean_a12_adjusted_loss_verification_20261023/attempt_2_failed_compile.json)
and predecessor are preserved. The exact two-step implementation repair has
[outside/root clearance](../results/review_a12_adjusted_loss_two_tactic_repairs_root_preclear_20261023.md),
and its third compile and separate audit passed. The earlier failed output
retains zero verification credit. Two unused-variable/section warnings remain;
no cosmetic rerun was made. Owned temporary proof work is absent.
The [logic map](../docs/a12_adjusted_loss_kernel_logic_20261023.md) retains
the geometry, actual fields/inventories, legal pressures and spectator floors
as premises. It does not formalize the later universal A12/A13 theorems or
the later arithmetic frontier refinements. The source header and logic map
are immutable preparation records; the execution receipt supplies the later
PASS. Older successful modules were not rerun or relabelled.

`HeavyCentreForestScalarKernel` has [bounded compilation and a separate
eleven-declaration actual audit PASS](../results/lean_heavy_centre_forest_verification_20261023/manifest.json).
Its [map](../docs/heavy_centre_forest_scalar_kernel_logic_20261023.md) separates
ordered-ring cone algebra, finite-chain induction from two local bridges and
incidence arithmetic from the actual geometric and graph instances. Six prints
use `propext`, `Classical.choice`, `Quot.sound`; the other five use `propext`
and `Quot.sound`. No `sorryAx` occurs. Compilation took2.4s and the separate
audit1.93s. The initial RAM refusal and unused-variable warning are preserved;
owned temporary work is absent. No old successful module was rerun.

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.HeavyCentreForestScalarKernel
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.HeavyCentreForestScalarKernel --audit-axioms
```

`A18TwoDottedEndpointScalarKernel` has [bounded compilation and a separate
ten-declaration actual axiom audit PASS](../results/lean_a18_two_dotted_endpoint_verification_20261023/manifest.json).
It derives the continuous radical gap, endpoint gains and one legally funded
completed square over a genuine ordered ring. All ten actual axiom sets contain
only `propext`, `Classical.choice` and `Quot.sound`; none uses `sorryAx`.
The [immutable preparation map](../docs/a18_two_dotted_endpoint_scalar_kernel_logic_20261023.md)
exposes actual geometric pressures, reference inventories, legal prices,
orthogonality and the real interpretation as premises. The runtime receipt
supplies the later status; preparation headers are not silently rewritten.
It supports the two-dotted A18 branch, not the full whole-A18 theorem.
Compilation took3.25s and the separate admitted audit1.77s. The earlier
admission refusal and one unused-variable warning are preserved; owned work
was removed. No successful earlier module was rerun.

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.A18TwoDottedEndpointScalarKernel
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.A18TwoDottedEndpointScalarKernel --audit-axioms
```

The later written compact18 frontier has ONLY wholeD18, after the
[complete written B18 exclusion](../results/review_compact18_whole_b18_complete_composition_and_single_d18_frontier_root_acceptance_20261023.md).
This concerns a selected farthest vertex, not every vertex or an example.
Neither the whole-A18 nor whole-B18 geometric theorem is formalized by
this two-dotted scalar module.

The [GLOBAL arithmetic15 proof](../results/review_arithmetic16_whole_d16_manual_dyadic_exclusion_and_global15_root_acceptance_20261023.md)
now closes the last whole-D16 family. Its earlier frontier, complete actual
quadratic head/fork reduction and dyadic splitting argument have full written
reviews. The new Std-only kernel above prepares a conditional discrete part;
it does not formalize those actual fields, all original walls, coefficient
tables or geometric prices. Older compiled/audited16 and precursor kernels
retain their historical scopes and acquire no new geometric instance.

The [fresh arithmetic15 necessary frontier](../results/review_arithmetic15_fresh_complete_twenty_two_group_frontier_root_full_source_acceptance_20261023.md)
has a full written review of its original22 grouped rows. The
[complete whole-B14+I proof](../results/review_arithmetic15_whole_b14_actual_i_all_degree_manual_root_acceptance_20261023.md)
removes original group19. The
[complete B13 rank-two proof](../results/review_arithmetic15_b13_ordinary_and_typed_seven_whole_group17_root_acceptance_20261023.md)
also removes original17. Ordinary spectators
use only the fresh6..10 body; actual typed-seven spectators have freshly
proved11/12 estimates in K7 before E12 is used. The
[complete D13 theorem](../results/review_arithmetic15_d13_both_actual_typed_seven_and_whole_group18_root_acceptance_20261023.md)
removes original18 too. Its separate
typed proof uses a full maximal order, exhaustive coefficient strips,
both signed fork lengths>6 and one reference omitting11. The new
handwritten field/fork argument has no Lean instance.
The [A13 typed-seven theorem](../results/review_arithmetic15_a13_both_actual_typed_seven_whole_group14_root_acceptance_20261023.md)
also removes original14, leaving18 necessary groups at that stage. Its new endpoint
and two-block reference arguments have full written reviews; the old
central .574 theorem remains an explicit input. There is no new Lean
instance or endpoint theorem for arbitrary extra spectators.
The [four-family extension](../results/review_arithmetic15_four_typed_seven_a_families_full_foreign_root_acceptance_20261023.md)
retains all extra spectator costs and removes original2,5,8,11, leaving14
necessary groups. The [two-A1 theorem](../results/review_arithmetic15_a13_two_ordinary_a1_root_full_source_acceptance_20261023.md)
narrows original12 to ordinary irreducible T2. Both are complete reviewed
written proofs with no new Lean instance or numerical14 theorem. The
[compact spectral tool](../results/review_spectral_frog_canonical_lorentz_gram_band_root_full_source_acceptance_20261023.md)
likewise has a full written proof of B³-B positive semidefinite for an
appropriately scaled FULL wall Gram. Unknown-wall control and a local
dimension application remain open; this has no formalization credit.
These new proofs and the whole-B14 joint scalar and
all-degree closure, small-tail H4 projection and fresh partition exhaustion
have no new Lean instance;
the older compiled frontier modules retain their stated dimensions.

The new [complete native-six finite budgets](../results/review_finite57_all_single_native_six_finite_budgets_root_acceptance_20261023.md)
cover all finite vertices containing native6, including Context and old
tokens. The same-rank replacement certificate, strict modular endpoint
and complete finite remainder bounds have full written review. These
and the weighted incidence/minima/conditional ideal comparison have no
new Lean instance. The [selected old ideal debts](../results/review_finite57_native_six_all_old_ideal_debts_root_full_source_acceptance_20261023.md)
are covered in a separately reviewed written argument. The residual
bare small-factor worker now has [actual generation, separate full replay
and complete ideal-budget acceptance](../results/review_finite57_native_six_complete_weighted_ideal_budgets_root_acceptance_20261023.md):
All 6,626 cases are strictly positive across 862 profiles. Every one of the 6,656
canonical record was regenerated with strict EOF, under two separately
admitted guarded units. The exact selected weighted ideal ledger is
closed by its conjunction with the old-debt theorem. This is same-
implementation execution, not a new Lean instance or silent additional
premise of any older compiled module. The k0 finite-heavy and global
face-coverage obligations remain.
The [new weighted strict-counting proof](../results/review_finite57_native_six_weighted_global_extrema_strictness_root_acceptance_20261023.md)
is complete written mathematics, with no new Lean instance. It closes
the global positive-correction step while retaining every face-payment,
ledger-conservation and finite/ideal cap premise.
The new B14+actual-I
[head](../results/review_arithmetic15_b14_all_finite_pure_heads_positive_foreign_root_acceptance_20261023.md),
[short-cap](../results/review_arithmetic15_b14_complete_short_white_isolated_dotted_root_acceptance_20261023.md)
and [cubic norm](../results/review_arithmetic15_b14_cubic_terminal_even_norm_foreign_four_root_acceptance_20261023.md)
reductions likewise have full written reviews and no new formal instance.
The later complete whole-B14+I proof removes group19 as written mathematics;
none of these results enlarges any older kernel's scope.
The [complete zero-exchanged far-row8 refund](../results/review_finite57_zero_exchanged_far8_all_corner_refund_one_unit_root_acceptance_20261023.md)
also has full written review. Its actual finite/affine classification,
fifteen Gram/native-gap cases and simultaneous full-ledger proof have
no new Lean instance; no earlier kernel is credited with this theorem.
The [positive-exchange far8 extension](../results/review_finite57_all_far8_positive_exchange_normalized_refund_root_acceptance_20261023.md)
shares one receiver port across the zero/positive union. The
[broader pure-D17 theorem](../results/review_finite57_whole_d17_large_inert_d_far8_refund_root_acceptance_20261023.md)
retains arbitrary inert D-factors of rank>=9 and total40, conditional
on an actual qualifying centre, including the old D40+D17 load.
Both have full written reviews and no new Lean instance or56 theorem.
The [entire D18 isolated/near3 family](../results/review_compact18_d18_whole_isolated_near_three_family_root_acceptance_20261023.md)
is also excluded as written mathematics: actual future-core and
A16+2A1-neighbor arguments retain all heads, unequal owns, both reverse
coordinates and one norm budget. It does not enlarge any older
compiled kernel or prove compact17.
The [complete isolated/far3 finite-mutual reduction](../results/review_compact18_d18_isolated_far_three_complete_finite_reduction_root_acceptance_20261023.md)
has full written source review and root acceptance. Actual B18/H3
points and the finite7+ continuous core initially left native3 or dotted mutuals,
and only dotted mutuals under native5 head. No old scalar kernel
formalizes that complete geometry or its new simultaneous S2 bound.
The [subsequent native3 theorem](../results/review_compact18_d18_isolated_far_three_native_three_two_dual_root_acceptance_20261023.md)
now leaves dotted mutuals only for every head. Its actual-A18 factor8/7,
joint reverse duals, fresh .71 inventory and single funded reference
are written and fully reviewed. They enlarge no existing Lean scope.
The later [whole isolated/isolated exclusion](../results/review_compact18_d18_whole_isolated_isolated_pair_root_acceptance_20261023.md),
[isolated native5 long-head conjunction](../results/review_compact18_d18_isolated_native5_head_complete_exclusion_root_acceptance_20261023.md)
and [six central bridges](../results/review_compact18_d18_central_bridge_six_configurations_root_full_source_acceptance_20261024.md)
are likewise complete reviewed written mathematics with their own original-vector
and all-wall boundaries. The new finite57 [row5](../results/review_finite57_positive_row_five_two_corner_refund_root_full_source_acceptance_20261023.md)
and [row11](../results/review_finite57_whole_positive_row_eleven_refund_root_full_source_acceptance_20261023.md)
refund theorems retain their complete geometric/ledger premises. None of these
new geometric results is added to an older kernel's formal scope. Their current
research frontier is summarized in the [comparison guide](../docs/current_bounds_and_verification.md).
The older five-case and six-case modules below retain their original verified
scope; they do not acquire the newer written exclusions automatically.

`GeneralEighteenEndpointRefinement` has a [bounded compile and separate
eleven-declaration actual audit PASS](../results/lean_general18_endpoint_refinement_verification_20261023/manifest.json).
It derives the five compact18 forest cases at its preserved checkpoint at the SAME vertex
from eleven explicit written interfaces. Its arithmetic section proves
thirty Taylor-table integer inequalities, finite Machin-sum bounds and
all endpoint comparison signs by ordinary kernel reduction. Eight prints
are axiom-free; the remaining three use only propext/Quot.sound. The
[logic map](../docs/general_eighteen_endpoint_refinement_logic_20261023.md)
exposes the geometric and analytic bridges still external. It reuses the
earlier capped six-case dependency without rerunning it, and proves no
whole-A18 exclusion, global compact17 bound or geometric instance.

Reproduce that module with two separate guarded commands:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.GeneralEighteenEndpointRefinement
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.GeneralEighteenEndpointRefinement --audit-axioms
```

`ArithmeticDimensionSixteenComposition` has bounded compilation and a
SEPARATE actual three-declaration axiom audit
[PASS](../results/lean_arithmetic16_composition_verification_20261022/manifest.json).
It derives whole-A17 exclusion, dispatches the fourteen original forest
groups, and combines exact-dimension17 exclusion with the prior global17
tail. Its [logic map](../docs/arithmetic_dimension_sixteen_composition_logic_20261022.md)
exposes27 geometric, arithmetic and computed-data premises; there is no
geometric instance or end-to-end paper formalization. The first two
declarations are axiom-free; the final Nat implication uses only propext
and Quot.sound, with no sorryAx. The initial pre-child RAM refusal and
precompile source repair are preserved. The successful audit reused the
compiled object. Frozen preparation headers describe their historical stage;
the actual receipt records the later execution.

To reproduce only this module in two serial batches, from the repository root:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.ArithmeticDimensionSixteenComposition
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.ArithmeticDimensionSixteenComposition --audit-axioms
```

The existing aggregate root does not automatically import this new module.
Do not rebuild already successful modules simply to refresh their status.

`ClassicalSeventeenCensusKernel` has bounded compilation and a SEPARATE
actual four-declaration axiom audit
[PASS](../results/lean_classical17_census_verification_20261022/manifest.json).
It derives the strict head gain, fork square comparison and1.00184
same-vector scalar contradiction. Its
[logic map](../docs/classical_seventeen_census_kernel_logic_20261022.md)
keeps actual geometry, field reductions, complete census, radical
interpretation and the legal paid-reference inequality explicit outside
Lean. Only standard axioms occur; this is not the complete D17 geometry.
The audit reused the new compilation, and no old module was rebuilt.

`D15TypedSevenCentreKernel` has bounded compilation and a separate actual
six-declaration axiom audit
[PASS](../results/lean_d15_typed_seven_centre_verification_20261022/manifest.json).
It derives the continuous31.59 radical gap and the conditional1.00357
and1.00069125 funded margins from primitive scalar premises. Only standard
axioms occur, with no `sorryAx` in the successful audit. Its
[logic map](../docs/d15_typed_seven_centre_kernel_logic_20261022.md)
keeps actual geometry, complete classification, field interpretation
and reference applicability outside Lean. A failed elaboration and its
sole tactic correction are preserved; the corrected module was compiled
once and reused for the separate audit.

`D18CentreScalarKernel` now has bounded compilation and a separate actual
five-declaration axiom audit [PASS](../results/lean_d18_centre_verification_20261022/manifest.json).
`CubicB16ShortCapKernel` likewise has six actual audited declarations and
its own [PASS receipt](../results/lean_cubic_b16_short_cap_verification_20261022/manifest.json).
Each imports only Std, uses standard axioms and has no `sorryAx`. Both
audit runs reused the successful compilation without rebuilding it; one
unused-variable/section warning per module is preserved. Their source
headers and [D18](../docs/d18_centre_scalar_kernel_logic_20261021.md)/[B16](../docs/cubic_b16_short_cap_kernel_logic_20261022.md)
logic maps keep the preparation-stage bytes. These later receipts record
the successful verification. Actual feasible-point geometry, field/basis
and radical bounds, complete cap coverage and the funded-reference bridge
remain explicit application premises.

`LorentzObstructionKernel` has bounded compilation and a separate actual
seven-declaration axiom audit
[PASS](../results/lean_lorentz_obstruction_verification_20261021/manifest.json).
It derives the second-orbit cutoff, the exact golden comparison, and the
restricted large-D Schur numerator/sign and rank contradiction. Only
`propext`, `Classical.choice` and `Quot.sound` occur; no `sorryAx` is used.
The [logic map](../docs/lorentz_obstruction_kernel_logic_20261021.md) retains
actual finite-order averaging/Cauchy realization, full chamber germs,
original-lattice preservation, arithmetic covolume and Gram/rank
interpretation as written inputs. It does not formalize the entire host
exclusion or change a dimension bound. Two failed predecessors are preserved
with zero credit; successful dependencies were not rerun.

`EndpointTriangleNormKernel` has bounded compilation and a separate actual
ten-declaration axiom audit PASS. Its
[receipt](../results/lean_endpoint_triangle_verification_20261021/manifest.json)
binds the final source, both serial ledgers and actual output; only standard
axioms occur, with no `sorryAx`. It proves cycle cancellation, the explicit
integral Bezout expression, physical/conjugate inequalities, finite norm
product contradictions and the rational positive-integer case. The
[logic map](../docs/endpoint_triangle_norm_kernel_logic_20261021.md) leaves
actual walls, cyclic integrality, field-norm interpretation and geometric
classification as application inputs. Two failed tactic predecessors have
zero verification credit and remain preserved. Two unused-section-variable
warnings are retained; no successful old module was rerun.

The selected `RegularIdeal24CellHeightKernel` and
`ClassicalDimensionComposition` now each have bounded compilation and a
separate actual printed-axiom audit PASS. The
[completion receipt](../results/lean_publication_completion_20261021/manifest.json)
binds both unchanged sources, four individual ledgers and four actual stdout
records. Regular24 has five actual declarations, using only `propext`,
`Classical.choice` and `Quot.sound`; Classical has four, using only `propext`
and `Quot.sound`. Neither uses `sorryAx`. The four runs were strictly serial,
with 2,048 MiB/one thread and unchanged resource gates; successful dependencies
were not rerun. Two Regular24 unused-section-variable warnings are retained;
Classical has no warnings.

Regular24 derives the explicit opposite-normal heights and coordinate-center
consequence; its actual-cell isometry and geometric interpretation remain
inputs. Classical derives the quartic integer norm contradiction and conditional
19/17 compositions; actual forests, geometry, field reductions and the other
written branch exclusions remain premises. Their
[Regular24](../docs/regular_ideal_24cell_height_logic_20261021.md) and
[classical](../docs/classical_dimension_composition_logic_20261021.md) maps
give the precise boundaries. Historical preparation headers in the Lean files
are retained to preserve the compiled source hashes.

The new per-manuscript interfaces have an
[explicit assumption register](../docs/manuscript_assumptions.md).
`ManuscriptIdeal4` now compiles in a capped batch and its actual axiom
prints contain only the standard axioms above. It derives the finite
defect-count implications, restricted-cut induction, chain counts and
second-volume equality conclusion **conditional on its named premises**.
The actual geometric objects, spectrum reduction and cutting/recognition
interfaces remain written proofs. The source-handoff register predates
this successful compilation; the
[final run](../results/ideal_ra4_unified_core_verification_final_20261018.json)
records the current status. `CompactArithmetic20`, `ManuscriptDataBounds`,
`ManuscriptCompact20` and `ManuscriptFinite57` have since compiled in
separate capped modules. Their two main composition audits list only
propext and Quot.sound, with no sorryAx. The
[compact20 core run](../results/compact20_core_verification_20261019.json)
passes all ten declared jobs, and the
[finite57 audit](../results/lean_execution_ledgers_snapshot_20261021/CoxeterProofs_ManuscriptFinite57_axiom_audit.json)
binds its fresh source. Actual whole-star/factor coverage, arithmetic
geometry and certificate interpretation remain typed premises. The
sharp right-angled interface has since passed its capped compilation and
12-theorem printed-axiom audit, as recorded in the table below.
The old root module does not automatically import this new paper interface.

|Module|What is proved|What it does not establish|
|---|---|---|
|CertificateBounds|Last-item inequalities imply an upper for every allowed finite forest, by induction; finite defect-sum lower bounds; soundness of the executable table check|Completeness of a geometric forest domain or the correctness of its weights|
|OptionalTableBounds|The same universal forest induction with explicit unreachable rank values; a `none` entry excludes every forest of that rank|The saved optional tables, their complete domains or geometric interpretation|
|Generated/Residual58|A frozen 149-item, rank0..58 residual table; upper378 at rank58 and a formal attaining witness|Other fixed-component/context branches, the whole finite maximum, or a dimension theorem|
|FiniteTableArithmetic; Generated/Finite59|All 5148 frozen residual tables at their used rank caps, including all 285643 required Bellman inequalities; the induction gives an upper for every permitted residual forest|Complete fixed/context domain coverage, the context-cover argument, or the geometric meaning of the weights|
|ForestDomination; Generated/TypeDomination59|All180 serialized actual/representative base pairs and32400 cross-cost pairs; induction proves ordinary-cost domination for every finite forest|The two guarded contextual costs, geometric interpretation, fixed-domain/context decomposition, or the whole finite bound|
|CompactHVector|Strict weighted h-vector inequality from paired incidence identities, nonnegative masses and positive endpoint mass; exact coefficient checks for dimensions21–29|The geometric f/h identities, local arithmetic rules, their complete domains or field bounds|
|CompactFieldBudget|The explicit sqrt-two component costs imply a vertex budget2n−2 for every rank n≥9, by a deficit proof with no rank cutoff|Identifying the explicit formulas with actual weighted diagrams and proving every triangle/quadrilateral payment|
|IdealProfileArithmetic; Generated/Ideal59|All1890 terminal ideal-link profiles: complete enumerated terminal domain, exact polynomial multiplication, vertex/index sums, rational envelope costs, normalized margin≥1/340 and absolute surplus≥91/2|The all-rank concentration theorem, affine-factor payment geometry, and cusp/global interpretation|
|IdealDomainCompleteness|Every ordered rank58 profile with at most one low factor below9 occurs in the checked data; the recursive partition enumerator is proved complete by induction|Concentration of an arbitrary actual cusp into that terminal form and geometric payment interpretation|
|ShapeCurvatureArithmetic|Each extremal h-vector contributes slack at least870; two extrema give scaled correction8700; finite curvature constants and the final payment/shape contradiction|The geometric index/incidence identities and identification of the two distinct extrema|
|IdealRA4Gap|The complete integer case split forces W>=12; excludes3,6,9; converts the gap to volume numerator at least8|Geometric construction of the selected-facet partition, cusp-shape normal argument, and the external facet classification|
|ManuscriptIdeal4|Conditional gap/equality composition; finite defect-list reduction; strong induction from one-step Oct/D cuts to cell counts; two-cell chain implies the double. Fresh capped compilation and explicit axiom prints pass.|Instantiation of the published/script premises and actual geometric interfaces; universal decomposition outside the Oct/D class|
|CompactArithmetic20; ManuscriptCompact20|Exact paired coefficients, endpoint550 surplus, scaled contradiction, and conditional composition of all13 fields at21 with the higher tail. Fresh capped compilation and printed-axiom audits pass.|Whole-star coverage, same-hub geometry, D5/H4 consumers and interpretation of the admitted script outputs|
|ArithmeticEightField19|Rank20 paired coefficients, endpoint500 surplus, ten-family quarter-unit cost reduction, elementary face-payment implications and conditional dimension19 composition for eight actual ground fields. Capped compilation and separate seven-theorem axiom audit pass.|The full actual-face consumer proof, role separation, complete forest coverage, actual charge/h-vector realization and earlier universal20 theorem remain explicit premises. No field construction or end-to-end geometric formalization.|
|ArithmeticK7Field19|All thirteen quarter-unit family budgets, vertex sum, endpoint500 surplus and conditional actual-K7 dimension19. Capped compilation and separate six-theorem axiom audit pass; printed dependencies only propext/Quot.sound.|Actual forest classification, geometric consumer payments, role counts, h-vector realization and earlier universal20 remain explicit premises. Imports the checked eight-field dependency and retains legacy CompactHVector provenance; no end-to-end geometric formalization.|
|RightAngledMixed11|Complete polygon/E/R coupling algebra, strict opposition lower184460/1619, exact cross-product gap5861741, contradiction with upper1420561/12500 and conditional dimension11 composition. Capped compilation and separate six-theorem audit pass.|Actual P3/links/Q8/incidence geometry, the strict weighted row, the full real-valued positive-part transport proof and published12 tail remain explicit premises. No actual analytic row is asserted integer-valued.|
|RightAngledMixed10|New primitive-count balance, defect/cusp eliminations, strict3722B>159195f+2794E, upper4B<=168f+3E, contradiction and conditional dimension10. One capped module and separate eight-theorem audit PASS; only propext/Quot.sound, exact comparison axiom-free.|Actual three-face/P3 geometry, weighted/incidence/error inputs and cleared aggregate analytic transport8B<=336f+3H are explicit premises. The eleven-face interface includes face heredity and published compact4/ideal6; no earlier internal ceiling is assumed.|
|RightAngledSharpDimensions|Current sharp8/general and4/entire-ideal conditional logic. Exact cubic/Gram coefficients, count27-to26 endpoint, three actual projection steps and both dimension implications. Capped compilation and12-theorem axiom audit PASS; only propext/Quot.sound in non-axiom-free dependencies.|Actual Lorentz projection, real height bounds, code interpretation, universal real tensor-square/equality argument, higher-face inheritance and published Dufour/Musin/Kolpakov inputs remain explicit. No integer proxy for a real height or norm; no geometric end-to-end Lean claim.|
|ManuscriptFinite57|Finite/ideal exact outputs, actual domination and strict4205 correction imply no rank58 charge model; the earlier58 tail then gives57. Fresh capped compilation and printed-axiom audit pass.|Actual Coxeter geometry, universal routing/concentration and interpretation of the index and certificate rows|
|DimensionComposition|The scaled ambient59 contradiction and tail composition; the compact three-field composition|Internal global geometric inequalities, local transfers, field bounds or their full certificates|
|IdealRA5|All-24-cell facet incidence identities contradict the dimension5 Euler equation|A general exclusion of ideal right-angled5-polytopes|
|FiveFaceCarrier|The ideal5 incidence/Euler identities force total κ48 and a positive carrier; mixed finite/ideal incidence inequalities combine correctly|Actual face/cusp incidence maps, the geometric RA4 classification, or a dimension exclusion|
|DoublingKappa|Exact κ recurrence, the16-cell path's counts, and the exponential lower bound under further doublings|The geometric halfspace-union lemma or validity of an actual reflection construction|
|GuardedRank58|Exact rank-capacity guards, the eight lowered and six increasing forest costs, the4205/36 contradiction and conditional57 tail composition|The D13 geometric theorem, the complete rank58 finite and ideal certificates, or their actual-polytope interpretation|
|RightAngledDensity|The56/170 ideal5 and227/3859 ideal6 consequences; four five-face/vertex-ratio eliminations and the3/40 six-face fraction|The actual-facet/cusp injection, geometric carrier theorem, or published average/incidence premises|
|RightAngledCuspCapacity|The reviewed59/180 ideal5 and240/4980 ideal6 integer consequences, exact outside-neighbor deficit, strict ideal5 incidence sum, updated4/43 and5/71 densities and mixed vertex ratios. Fresh capped compilation and axiom audit pass.|The geometric at-most-two-neighbor lemma,Figure7/large-facet bridge,actual injection,box/Euler interpretation,published minima and universal transport remain explicit inputs. No dimension exclusion.|
|RightAngledCuspDefect|Positive cusp-defect count logic and the stronger62/189 ideal5,252/5504 ideal6 consequences, with2/23 and1/15 density bounds. Fresh capped compilation and axiom audit pass.|The spinor polygon lemma,actual good-opposite geometry,triangle/outgoing-edge bijection and geometric interpretation of the count interfaces remain written premises,along with inherited carrier and transport inputs. No dimension exclusion.|
|RightAngledCuspPropagation|The ideal5 facet-propagation integer logic and62/206,269/5875 consequences; updated15/242 six-face density. Fresh capped compilation and axiom audit pass.|The exact outside-vertex identity,spinor polygon counts and9G≤2N0 actual-opposite incidence remain explicit geometric inputs,along with inherited carrier and transport premises. No dimension exclusion.|
|ContextCover59|The complete0/1/2-reserved-A1 presence-guard cover preserves both contextual charges, with no witness multiplicity|Fixed rank-domain completeness, item identities and the connection of contextual choices to the saved table indices|
|RestrictedCarrierCounts|Cusp-count equations force divisibility by6; the orthogonality graph inequality forces30 carriers and66 doubles; the457/2397 arithmetic assignment is checked|The width and decorated-box geometry, graph clique bound, or any geometric realization|
|Concentration58|All28 low8 transfers recomputed from the two uniform distributions, with nonnegative absolute-moment difference for each of58 possible remaining-degree shifts; the signed moment stays nonnegative for arbitrary nonnegative convolution coefficients|Iteration to canonical ranks, rational margin/product propagation and the actual cusp interpretation|

Among the earlier verified successors are these three kernels, with capped
compilation and actual printed-axiom audits covering28 selected declarations:

| Module | Additional verified logic | Explicit remaining inputs |
| --- | --- | --- |
| RightAngledCubicScalarKernel | Coordinate Gram positivity, cubic sum of squares, count bound and strict endpoint; constructs the older cubic interface. | Actual code/tensor realization and its two contraction identities; Lorentz projection, heights, face inheritance and the named published theorems. |
| ConcentrationDescent58 | Mass/potential/product laws,1596-step termination, terminal recognition and sorted rank57 profiles; normalized-to-absolute surplus lifting. | Real convex order and normalized comparison; actual cusp interpretation; exact terminal enumeration and margins. |
| CompactFarthestComposition23 | Real local energies imply the scaled charge budget; dominant uniqueness, A/B exclusions, sole D23/I frontier and conditional ceiling23. | Actual component estimates/caps, rank/debit extraction, the coupled geometric obstruction and published29. |

Only `propext`, `Quot.sound` and `Classical.choice` occur in the audited
declarations; none uses `sorryAx`. The
[verification receipt](../results/publication_logic_verification_20261020.json)
binds current sources, objects, import/compiler ledgers and parsed actual
axiom output. See the [logic map](../docs/publication_logic_formalization_20261020.md)
for the precise boundaries and serial reproduction commands. This is
selected internal logic, rather than complete geometric formalization.

The finite-forest checker uses `by decide`, which produces a proof checked
by the Lean kernel; it does not invoke the Python producer or trust native
evaluation as a new axiom. The generic induction is proved separately.

Reproduce the imported residual certificate from the repository root:

```sh
python3 scripts/export_lean_residual58_20261018.py
python3 scripts/export_lean_finite59_20261018.py
python3 scripts/export_lean_ideal59_20261018.py
python3 scripts/export_lean_type_domination59_20261018.py
```

This checks an exact export against the frozen source hash and manifests.
`--write` regenerates the Lean text and source metadata. The Lean source
then checks the actual numerical inequalities independently. The metadata
also records precisely which original type each item represents.

The complete ideal-profile kernel build passed in688 seconds on the
current environment. It uses ordinary `decide` and recomputes the actual
polynomials rather than trusting the saved coefficients or Python output.
Cached builds are much quicker; a clean build must repeat that work.

## Assumption boundary

The user permits published external theorems to be assumed. The registry
in [assumptions.json](assumptions.json) names those statements separately
from unformalized internal facts. At present the Lean modules use explicit
arguments rather than introducing geometric axioms. For example,
`defect_gap_of_partition` proves its conclusion from a selected/rest defect
partition; the new geometric construction of that partition still needs a
Lean proof. Likewise the field-composition theorem does not by itself
prove the internally new field-specific dimension bounds.

Consequently **none of the complete geometric dimension/volume theorems is
yet labelled end-to-end Lean-verified**. Successful compilation is genuine
verification of the stated module theorems, with their explicit hypotheses;
it is not a substitute for proving the remaining internal bridges.

The manuscripts and per-result lemma maps document the next obligations.

## Compact17: actual compilation and separate axiom audits

The general compact17 written theorem assumes no arithmeticity. Its new
[scalar source](CoxeterProofs/D18UniformDottedHeadScalarKernel.lean) and
[composition source](CoxeterProofs/CompactDimensionSeventeenComposition.lean)
have complete different-author static reviews and
[root pre-execution clearance](../results/review_compact17_scalar_composition_lean_sources_root_preexecution_clearance_20261024.md).
They use Std only, with13 public theorems and ten selected axiom-print roots.
No custom axioms or proof placeholders are declared. Both now have
[actual bounded compilation and separate actual axiom audits](../results/compact17_targeted_lean_completion_20261025/manifest.json):
seven scalar declarations and three composition declarations. The earlier
[resource deferral](../results/compact17_lean_source_clearance_and_resource_deferral_20261024.json)
is preserved as history. The scalar's first compile failed on the reserved
binder `mutual`; an exact three-occurrence alpha rename to `mutualProduct`
preceded its successful retry. Two unused-section-variable warnings remain.

The [logic map](../docs/compact_dimension17_conditional_logic_map_20261024.md)
lists calculus/radical estimates, actual-wall points, coefficient funding,
complete geometric inventories and the prior global18 tail as explicit
parameters. Scalar theorems derive their conclusions from lower-level
Taylor and value/derivative bounds; the composition derives whole-D18
exclusion and17 from a single coherent geometric record. The scalar and
geometric records are not automatically instantiated or wired together.
Printed axioms do not list these ordinary theorem parameters. Older
successful modules retain their original scope and have not been rerun.
