# Lean verification of the internal proof logic

This is a real Lean4 project, pinned by `lean-toolchain` to the locally
installed Lean4.32.2. It currently uses only Lean's standard library, so
no unpinned mathematical library or network download is required.

Run these commands from the repository root:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py --check-only
python3 -B scripts/run_lean_batches_signal_safe_20261018.py --audit-axioms
```

The runner holds an exclusive lock and compiles one module at a time with
one thread, a 2048MiB Lean memory cap and a 15-minute timeout. It requires
40GiB free disk and 6GiB available memory before each child starts, and
stops below 8GiB disk or 3GiB available memory. Do not use unrestricted
`lake build` or direct uncapped Lean commands. A named module can be supplied
to verify only that module and its local dependencies.

The read-only `--check-only` command reports cache provenance. Existing
objects without a guarded-build ledger are labelled
`LEGACY_CACHED_UNCAPPED`; their reuse does not claim a new capped verification.
New objects have streaming hashes of sources, imported local objects,
compiler/configuration and the resulting object recorded in
`formal/.proof-batches/`. Changed inputs prevent publication; a mismatched
ledger prevents cache reuse. Output retention is bounded to 512KiB.

The runner removes temporary Lean/C/IR/LSP and TeX outputs after every
batch, including failure and entirely cached runs, while preserving proof
objects, sources, certificates and final PDFs. Below 20GiB free disk it
writes `.proof-batches/cleanup_needed.json`: the coordinating agent must
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
