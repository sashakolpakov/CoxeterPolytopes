# Verified logic after the six source audits

20 October2026. These seven small kernels now compile under the bounded
Lean runner. Their **62 actual axiom prints pass**: ten for the cubic,
ten for concentration, eight for compact23 accounting, eight for the
residual Bellman certificate, eight for projection/height descent and nine
for finite-volume tail arithmetic, plus nine for dimension58 shape counting.
The former six-kernel/53-print receipt is preserved. Only
`propext`, `Quot.sound` and `Classical.choice` occur; there is no
`sorryAx`. The [verification receipt](../results/publication_logic_verification_20261020.json)
binds current sources, objects, compiler/import ledgers and parsed actual
axiom output. The catalogue now includes these modules in verified
conditional coverage. Existing frozen sources remain unchanged.
The explicit inputs below remain assumptions; this is not an end-to-end
formalization of the geometric manuscripts or a new dimension bound.

All seven kernels also have directly inspectable raw axiom output and
fourteen immutable compilation/audit snapshots. The
[two-output recapture receipt](../results/publication_logic_two_raw_output_recaptures_20261020.json)
restores the concentration and compact23 bytes missing from the original
package: both match their original immutable audit hashes exactly. These
were cached audit runs, with no recompilation or additional declarations.

| Verified kernel | Internal logic checked | Remaining inputs |
| --- | --- | --- |
| [RightAngledCubicScalarKernel](../formal/CoxeterProofs/RightAngledCubicScalarKernel.lean) | Derives finite Gram positivity from coordinate squares, then the cubic SOS/count/equality arguments and the older kernel's code-consequence fields. | Unit six-coordinate vectors and their strict pair gap; rank-two/three coordinate arrays with two explicit contraction identities. Actual wall/radial realization, face inheritance and published dimension inputs remain separate; the new projection kernel derives the scalar projection and height steps. |
| [ConcentrationDescent58](../formal/CoxeterProofs/ConcentrationDescent58.lean) | Mass/potential/product laws,1596-step termination, terminal recognition and sorted rank57 profile with unchanged link product; normalized-to-absolute gap lifting. | Actual cusp interpretation; convex-order/independence yielding normalized comparison; exact terminal enumeration and margin. |
| [CompactFarthestComposition23](../formal/CoxeterProofs/CompactFarthestComposition23.lean) | Derives the integer charge budget from real local inequalities, then dominant-component accounting, the A/B exclusions, sole D23+I frontier and conditional23. | Actual orthogonal energy decomposition, local energy bounds/caps and rank/debit extraction; exact D23/I obstruction and published29. |
| [ResidualBellmanCertificate](../formal/CoxeterProofs/ResidualBellmanCertificate.lean) | Derives list upper bounds, attainment and exact maxima, unreachability equivalences, permutation invariance and a rank58 witness-length bound from local transitions and tight predecessors. | Complete allowed residual alphabet and integer weights, interaction-free context split, concrete memo certificate and actual geometric domination. |
| [RightAngledProjectionHeightKernel](../formal/CoxeterProofs/RightAngledProjectionHeightKernel.lean) | Derives raw/normalized Lorentz coordinate identities, strict wall signs, exact orthogonal heights, maximal squared recurrence, finite descent and scalar radial gap. | Positive normalizers, actual finite unit walls/interior and maximum choices, future-sheet/facet interpretation, ridge-neighbor completeness, face inheritance and radial/tensor realization. |
| [FiniteVolumeTailKernel](../formal/CoxeterProofs/FiniteVolumeTailKernel.lean) | Checks all 105 literal finite-cap comparisons by kernel `decide`; derives scaled payment subtraction and near/far Lorentz scalar contradictions. | Actual charges and shape correction, complete finite/cusp table interpretation, framing walls and Lorentz geometry, near/far loss estimates and the compact-facet alternative. |
| [FiniteVolume58ShapeKernel](../formal/CoxeterProofs/FiniteVolume58ShapeKernel.lean) | Derives minimum uncut-edge assignment, both 841 gains and the distinct-extrema 8410 correction, finite 754 and the tagged global shape bound from primitive index splits/signs/incidence. | Actual links/index distribution and uncut-edge incidence, positive extremal h-coefficients, degree balance and common-denominator interpretation. The old `shapeLower` record is not automatically instantiated. |

## Sharp right-angled cubic

The frozen `RightAngledSharpDimensions` interface assumes
`StrictCodeSixConsequences.tensorCount` and its equality-endpoint
consequence. Its existing12-theorem audit checks their use. The new
source derives both from earlier scalar data instead of introducing
them again as hypotheses.

The scalar type is any linearly ordered commutative ring supported by
Lean's ordered-ring interface. It is not restricted to rational points.
For a code of size `n`, `ScalarCode` supplies its Gram entries `g`, unit
diagonal, `4g<1` off the diagonal and the ordinary four-point positivity
inequality, without an `n=27` condition. `TensorContractions` supplies
kernels `A,T`, nonnegative full pair sums for `g,A,T`, and precisely

```text
6 A_ij = 6 g_ij² - 1
8 T_ij = 8 g_ij³ - 3 g_ij.
```

The source derives `sum(16F)=n²+6sum(g)+12sum(A)+16sum(T)`.
The diagonal is27 and each off-diagonal cubic is nonpositive, yielding
`n²<=27n`. At equality each such cubic vanishes; its strict-domain
factorization forces `2g+1=0`. Four distinct indices then contradict
the unconditional four-point positivity inequality. The old Nat
endpoint argument supplies `n<=26`. No desired SOS, cardinality bound,
endpoint implication or dimension conclusion is an input field.

The principal declarations are `cubic_sum_of_squares`, `tensor_count`,
`endpoint_half_products`, `no_endpoint_twenty_seven`,
`strict_code_count_le26` and `cubic_geometry_from_scalar_bridge`.
Real tensor construction, its contractions and their geometric
interpretation are still obligations. The generic scalar theorem does
not instantiate Lean real numbers or formalize Lorentz geometry.
Musin and Kolpakov remain ideal-branch inputs; this cubic argument uses
neither. Initial compilation exposed three elaboration errors in casts,
a redundant tactic and an implicit scalar argument. Those were repaired
without changing the mathematical statements. The current module and
its ten declaration audits pass.

The preferred successor interface is now `CoordinateWitness`. It gives
actual finite coordinate arrays rather than PSD conclusions.
`sumN_swap` proves finite Fubini, and `coordinate_gram_square` rewrites
the full Gram sum as the sum of coordinate squares of the vector sum.
`coordinate_gram_nonneg` therefore derives all three positivity fields.
`coordinateDot_four_square` derives the four-point condition from the
unit six-coordinate vectors, including repeated indices. The rank-two
and rank-three arrays have36 and216 coordinates; only their two
pointwise contraction identities remain as algebraic/script inputs.
`coordinate_code_count_le26` and
`cubic_geometry_from_coordinate_bridge` use this stronger route. No
PSD, SOS, endpoint or cardinality premise is supplied to that route.
The earlier scalar interface remains as its internal reusable target.

## Lorentz projection and height descent

The [projection/height kernel](../formal/CoxeterProofs/RightAngledProjectionHeightKernel.lean),
source `2fac8028…` (282 lines,12333 B), now passes bounded compilation
and eight actual axiom prints. Its finite Lorentz form is defined from
the imported coordinate sums. From `q.q=-1`, `n.n=1`, `q.n=-t`, it
derives that `v=q+t*n` has `v.n=0` and `v.v=q.v=-1-t*t`.
Positive `lambda` satisfying `lambda^2*(1+t*t)=1` gives the derived
unit norm and negative pairing of `lambda*v`. Primitive strict old-wall
and nonpositive inter-wall products give strict new-wall products.
For orthogonal neighbors, the exact new height is `lambda*h`; an old
maximum `h<=t` yields the cleared squared recurrence.

The finite induction proves `j*a_j<1` from the local recurrence
`(1+a_j)*a_(j+1)<=a_j`, using no division or square roots. A generic
positive parameter `p` turns `(p-1)*h_i^2<1` and the unit-wall/all-wall
equations into the strict radial gap `p*z<1`. Thus the one-/three-step
scalar height and radial comparisons are derived rather than supplied
as final projection/code conclusions.

The [assumptions map](right_angled_projection_height_logic_20261020.md)
retains actual wall/incidence and future-sheet interpretation, positive
normalizer existence, ridge-neighbor completeness, inherited face geometry
and real radial/tensor realization. The older abstract geometry records
are not yet instantiated. The [actual output summary](../results/RightAngledProjectionHeightKernel_axiom_summary_20261020.json)
checks all eight qualified names, each using only `propext`,
`Classical.choice` and `Quot.sound`, with no `sorryAx`. Its1057 bytes
match both runner output hashes. First compilation took3.47s and the
initial separate audit3.15s. The current six-job sharpRA run repeats
its eight-print audit with identical child bytes. The publication logic receipt
binds immutable compilation/audit records;
operational cache paths remain metadata. The complete project-local
object closure is bound; standard
library objects remain within the pinned toolchain boundary.

## Low-eight concentration

Ranks are natural-number lists. Zero entries are padded factors: their
mass and potential contribution is zero and their link-size factor is
one. A `Step` moves a pair `1<=a<=b<8` to `(a-1,b+1)`, up to arbitrary
permutation. It is not limited to one chosen greedy strategy.

The potential is the sum of squared ranks. Each step increases it by
at least two. Mass57 gives `57<=potential<=57²=3249`, hence
`chain_length_at_rank57` bounds **every** finite chain by1596.
Strong induction on `3249-potential` constructs a terminal state; the
proof uses the successor potential bound to control truncated Nat
subtraction. `terminal_iff_lowCount` now identifies no outgoing `Step`
with at most one positive rank below8, counting repeated factors with
their multiplicities. `terminal_rank57_sorted_profile` removes neutral
zeros and constructs a permutation of one optional low factor followed
by sorted ranks at least8. Its total mass is57 and its link product is
unchanged. This is the terminal domain's mathematical recognition;
complete membership in the recorded finite enumeration is still a
separate input.

`absolute_surplus_lift` uses the denominator-cleared comparison

```text
vp*kq + vq*cp <= vq*kp + vp*cq,
```

where `v` denotes link size and `k,c` are curvature and charge after
clearing one common positive denominator. Together with positive
terminal size, `vq<=vp`, and `cq+gap<=kq`, it derives
`cp+gap<=kp`. It uses no truncated charge subtraction and assumes no
integer-valued real margin. `chain_absolute_surplus` obtains the volume
inequality from the pair chain. The normalized comparison is a
whole-chain premise; the real convex-order argument has not thereby
been formalized. A filter-normalization mismatch found at compilation
was repaired without changing its theorem. The current module and its
ten declaration audits pass.

## Compact23 accounting

The ledger's integer charges are certified rational **lower bounds**
multiplied by700. Actual projection energies remain arbitrary real
numbers. `RealAccounting.RealLedgerWitness` now takes local scaled
lower inequalities and a strict real energy sum below one.
`dominant_lower_sum` sums the local component inequalities, and
`strict_scaled_budget` derives the integer budget below700 using
order-preserving Nat casts. Repeated isomorphic blocks retain separate
energy positions. `RealLedgerWitness.toLedger` constructs the older
integer ledger without rounding any actual energy, and
`ceiling23_from_real_component_bounds` reuses its accounting proof.
The ledger records ordinary paid rank, at most one
unfunded short-A1 or rank-two-seven defect, and the dominant A/B/D
blocks. Type-J pairs are already paid in ordinary rank.

From the budget, `at_most_one_dominant` and
`unique_dominant_of_large_dimension` derive uniqueness/existence for
`n>=24`. The common debit accounting yields
`charge+50*spectatorRank<750`. Local caps use the explicit published29
input. The A half and five-sevenths floors force A20 with at least four
spectator ranks, contradicting its stronger charge. The B floors and
prefixes exclude all its cases. `single_D_frontier` derives exactly
dimension24, a sole D23 block and one short type-I A1 with no other paid
or seven-defect rank.

Only that derived configuration is fed to the final coupled geometric
obstruction. The bridge asks for a ledger only in hypothetical
dimensions at least24, where the high-dimensional component lemmas
apply; lower dimensions are handled directly. Neither a dominant-type
classification nor a dimension23 conclusion is supplied as a premise.
The declarations `ceiling23_from_local_component_bounds` and
`compact_dimension_le23` express the final conditional implication.
Root's source reading accepted the accounting and requested these two
scope refinements before compilation. The singleton witness needed one
elaboration repair. The current module and its eight declaration audits
pass.

## Residual Bellman certificate

The additional [ResidualBellmanCertificate](../formal/CoxeterProofs/ResidualBellmanCertificate.lean)
source `cd2a05a6…` (8274 B,202 lines) now passes bounded compilation and
all eight actual axiom prints, each using only `propext` and `Quot.sound`.
Its
[assumptions map](residual_bellman_certificate_logic_20261020.md) separates
the additive integer recurrence from alphabet/context completeness and
actual table validation. The finite57 wording now makes its certificate
premises precise: tight predecessors at every positive stored finite
entry, and a no-predecessor characterization only at positive ranks;
rank zero is the reachable base. Its [publication preparation receipt](../results/residual_bellman_publication_logic_preparation_20261020.json)
preserves sourcec387 and establishes that the current successor changes
only header comments; the proof body is unchanged. The [actual axiom
summary](../results/ResidualBellmanCertificate_axiom_summary_20261020.json)
binds all eight qualified declarations and the802-byte saved Lean output,
which matches both execution-ledger output hashes. Compilation took0.54s
and the separate audit0.41s under the one-thread/2048MiB runner.
The local dependency-object map is empty because the sole import is
`Std`, supplied by the pinned Lean4.32.2 toolchain. Compiler/configuration/
shared-library fingerprints are recorded; a separately hashed transitive
Std object closure is not claimed. Concrete alphabets, context separation,
actual memo-table checks and geometry remain application inputs. No
scientific checker was executed by this metadata update.

## Finite-volume tail arithmetic

The [FiniteVolumeTailKernel](../formal/CoxeterProofs/FiniteVolumeTailKernel.lean),
source `b8f03fd4…`, has actual bounded compilation and nine parsed axiom
prints in its [execution record](../results/finite_volume_tail_kernel_execution_20261020.json).
Its literal table matches all finite caps for 59–163 and proves their
comparisons with the parity-dependent shape coefficient by kernel `decide`.
Separate finite/cusp cost lists and a positive correction yield payment
subtraction; an additional natural scale clears individual denominators.
From the cleared Lorentz identity and loss estimates it derives
`9p>=M(N-37)` or `9p>=M(N-43)`, then contradicts `p<14` at N≥163,
including the boundary 20M=21. No dimension ceiling is supplied as a premise.
The [assumptions map](finite_volume_tail_kernel_assumptions_20261020.md)
and [independent source review](../results/review_finite_volume_tail_kernel_source_composition_20261020.md)
retain actual geometric/table bridges and the compact-facet branch as inputs;
this is not the complete n≤58 geometric tail.

The [dimension58 shape kernel](../formal/CoxeterProofs/FiniteVolume58ShapeKernel.lean),
source `40d12b57…`, has actual 0.59s compilation and a 0.55s nine-print audit;
its [summary](../results/FiniteVolume58ShapeKernel_axiom_summary_20261020.json),
[execution](../results/finite_volume58_shape_kernel_initial_execution_20261020.json)
and [assumptions map](finite_volume58_shape_kernel_assumptions_20261020.md)
retain the eight harmless lint warnings and explicit geometric/normalization
boundary. It replaces primary finite57's ambient59 `ShapeCurvatureArithmetic`
slot, preserving that source/evidence and historical atlas selection;
the 11 core/123 full job counts remain unchanged in the
[current plan](../results/publication_shape_counting_verification_plan_20261020.json).
The [original Khovanskii audit](../results/review_finite57_published_shape_inputs_20261020.md)
checks index/uncut-incidence scope. Separate [facet/acuteness](../results/review_sharp_ra_published_inputs_dufour_20261020.md)
and [ideal-branch](../results/review_sharp_ra_published_inputs_ideal_branch_20261020.md)
audits check the sharpRA inputs: Dufour's original is paywalled, with the
stated mixed 27 consequence checked in Alexandrov; Musin, Kolpakov and
Vinberg were checked in their originals.

Separately, the [two ancestor-forest replays](../results/finite57_ancestor_forest_replay_execution_summary_20261020.json)
actually checked 8114 tables and 10,195,241 inequalities in 5.106s,
with 318,799,872 B maximum process-tree RSS. They add concrete coarse and
first increasing-family coverage under their
[contract](../results/finite57_ancestor_forest_replay_contract_20261020.md),
not whole ancestor regeneration. The
[finite57 composition](../results/review_finite57_theorem_composition_adversarial_20261020.md)
and [sharpRA foundations](../results/review_right_angled_sharp_foundations_adversarial_20261020.md)
reviews independently check the remaining written mathematical interfaces.

## Reproduction

All six current PDFs match their completed sources, as recorded in the
[PDF receipt](../results/six_paper_pdf_completion_20261020.json).
Finite57's revised source `c300db69…` now has a48-page PDF (`1f3c020c…`),
with21 underfull spacing warnings and no overflow or unresolved references.
The44-page `89bd7658…` source/PDF and41-page phase remain preserved.
Compact23 source `2718192d…` / PDF `39e665d8…` remains36 pages with zero
warnings after its [one-phrase correlation correction](../results/compact23_correlation_precision_20261020.json)
to "at most1/3"; the strict demand is unchanged. These TeX builds add no
script replay or Lean verification. [Selected finite57 pages5,16,24,25](../results/finite57_publication_delta_visual_review_20261020.json)
passed visual inspection; no all-page or new compact23 visual inspection
is claimed.
The [bridge receipt](../results/finite57_publication_bridge_completion_20261020.json)
records the reviewed changes and actual bounded build.
Recheck the logic one module per batch with the unchanged signal-safe
runner, requesting its printed axioms:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.RightAngledCubicScalarKernel --audit-axioms
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.ConcentrationDescent58 --audit-axioms
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.CompactFarthestComposition23 --audit-axioms
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.ResidualBellmanCertificate --audit-axioms
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.RightAngledProjectionHeightKernel --audit-axioms
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.FiniteVolumeTailKernel --audit-axioms
```

These commands have passed for the current source hashes. The same
40GiB disk/6GiB available-RAM start gate, exclusive lock,
one-thread/2GiB Lean cap and cleanup apply to reproduction. Actual
axiom output was parsed and checked against the three standard axioms,
and source/project-local-import/object hashes were checked before coverage
advanced, with the toolchain import boundary described above.
The [current verification receipt](../results/publication_logic_verification_20261020.json)
records those bindings. The [earlier bridge receipt](../results/publication_logic_bridge_completion_20261020.json)
records the previously uncompiled versions, not the current status. The
[earlier preparation receipt](../results/publication_logic_and_replay_preparation_20261020.json)
is preserved; its three original drafts have byte-verified gzip archives
under `formal/source_archives/`.
