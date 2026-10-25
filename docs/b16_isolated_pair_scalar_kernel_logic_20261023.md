# B16 isolated pairs: scalar logic and external input map

23 October2026. **NO successful compilation or actual axiom audit is
credited yet.** The first admitted named compilation failed because two
local proof names used the reserved token `prefix`. The successor changes
only those names to `middle_gain`; it requires different-author delta
review and root clearance before a new named bounded compilation and
a SEPARATE actual axiom audit. The mathematical whole-B16 exclusion has
its own accepted written/computational boundary below. Its acceptance
does not depend on treating this uncompiled source as verified.

Module: `formal/CoxeterProofs/B16IsolatedPairScalarKernel.lean`,
12,112 bytes/282 lines, SHA256
`aa1f20bed7ae1dffc50423c72d5828eaa96cd9742dbc20a132020edca6207138`.
It imports ONLY `Std`. The pinned toolchain is `leanprover/lean4:v4.32.2`,
`formal/lean-toolchain` SHA256
`2bdc48adfa58d0017e538a0ad117c5d73d35deec879978f909406a80c8037273`.
No successful older module, dependency or frozen paper was changed.

Root authored this implementation and supplied the final quadratic
reference strategy. Replay authored the written cubic proof, census and
whole-family composition. Their written outside reviews are not a
different author's review of this NEW Lean source. The cubic-only
uncompiled predecessor was byte-preserved in
`results/source_archive/b16_scalar_before_cubic_quadratic_consolidation_20261023/manifest.json`;
it earned no execution credit. The consolidated source covers both
scalar branches within the existing arithmetic paper family.

## Literal mathematical and computation inputs

| Input | SHA256 | Imported boundary |
|---|---|---|
|results/arithmetic16_whole_b16_cubic_isolated_pair_dyadic_square_and_single_reference_exclusion_candidate_20261023.md|c20d7d91a35945bcbe65e5b327dc46204d0c2a53eabca325a51f85085b87582c|Entire dyadic argument, small-conjugate cubic lemma, original edge pressures and one reference|
|results/review_arithmetic16_whole_b16_cubic_dyadic_square_single_reference_root_acceptance_20261023.md|bf4dbd13fed5da4a367bdc7360d255406b9a7b3c270dbe50f87edd906120ffad|Accepted written cubic branch and full different-author review at its declared interface|
|results/arithmetic16_whole_b16_complete_quadratic_census_and_actual_edge_reference_exclusion_candidate_20261023.md|5328fa8e04606a209242c1dd9e254dc026f04586e3fa4fc030b22de1cff9243c|Entire whole-family composition and the sole quadratic pair's original-edge pressures|
|results/review_arithmetic16_whole_b16_census_runtime_and_complete_exclusion_full_scope_compact_reviewer_20261023.md|47486689bc0359756b4fcacea0b4dcb9479654d3f92531fbcf565d4ca3f25fe5|Full written/runtime outside review; not review or execution of this Lean module|
|results/review_arithmetic16_whole_b16_complete_census_and_reference_root_acceptance_20261023.md|fcb3ae95f196e654a67f5eedbed30368afbcacd2f8bffe6f6fc0354b2d71a0fe|Accepted wholeB16 exclusion; ONLY original group15 removed, ceiling16 unchanged|
|results/b16_isolated_quadratic_pair_census_verification_20261023.json|7db61f691ced96ffb9113dc98feb7f1b06e73a64e52aa0ab50a37136c45c4d70|Four actual guarded units, complete43-field/273-pair domain, ONE necessary pair; not an empty census|
|results/arithmetic16_whole_b16_complete_short_white_and_isolated_dotted_field_reduction_candidate_20261023.md|4fb202752c01f8ace7f38763403ef93dc13f937ddf1ee31cc656cdf80c27a514|Complete original head/short inventory, actual degree2/3 and unconditional strong6|
|results/review_arithmetic16_whole_b16_complete_short_white_isolated_dotted_field_root_acceptance_20261023.md|c987452b13d48fe7876b772785471c4c3645253c5592c13d347efcaf136e735f|Accepted original geometry/field reduction before any isolated-pair census|
|results/compact_farthest_vertex_mixed_replacement_inequalities_20261020.md|3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84|Actual all-wall comparison and one original vector with squared norm below1|

Hashes identify the exact mathematical interfaces; they do not prove
their claims in Lean. The module neither imports the census JSON nor
constructs actual mirrors, number fields, ideals or a real geometric
instance. Its scalar conclusions are conditional on the premises below.

## Types and external facts kept visible

Metric scalars use a generic linearly ordered commutative ring `R` with
`Lean.Grind.OrderedRing`. They are not integer proxies for real lengths.
Only the parity argument uses `Int`, because ideal valuations really are
integers. `CubicData R` is ordinary data passed as a theorem parameter,
not a typeclass or an automatically supplied instance.

The external number-field layer supplies the dyadic valuation identity
v=(e-1)+1-2e-(e-1) and, for a square, v=2k with integer e,v,k.
The kernel derives their inconsistency. The existence of the prime,
norm2 ideal, ramification values, ideal valuations and the reduction to
this exact identity remain external; oddness is not a substitute for
those field arguments.

For the cubic lemma, the actual integral W=U/2 has physical root w>1
and two distinct conjugates0<a<b<1. Its midpoint m and scalar trace S,
norm A and shifted norm B satisfy

    2m=a+b, S=w+a+b, A=wab,
    B=(w-1)(1-a)(1-b).

Positive integral norms supply the explicit alternatives A=1 or A>=2
and B=1 or B>=2. Integral trace supplies S<=5 or S=6 or S>=7.
These discreteness statements, separability, total reality and actual
field interpretation remain EXTERNAL ordinary fields of `CubicData`.
Neither w>5, positivity of a polynomial interval, an assumed minimal
polynomial coefficient nor the conclusion of the cubic lemma is an input.

For each prefix contradiction, actual geometry supplies scalar x1,S6,S16,
energy<1, the literal dot-product identity and ONE raw nonnegative norm
square. It also supplies the following primitive original-vector bounds:

| Branch | Head | Middle prefix | Short prefix |
|---|---|---|---|
|Cubic|1000*x1>618|100000*S6>201815|10000*S16>35276|
|Quadratic sole pair|100*x1>61|100000*S6>201815|100*S16>356|

Deriving these from original edge cosh lengths, the actual field/census,
the cubic w>5 conclusion, radical comparisons and strong6 is external.
The two scalar branch contradictions are separately derived; there is
no automatic field-case wrapper or claim that the Lean source realizes
the original quadratic/cubic alternatives. In particular the census
correctness/completeness is not asserted by a hidden Lean instance.

## Algebra actually implemented

All names below belong to `CoxeterProofs.B16IsolatedPairScalarKernel`.

| Declaration | Derived assertion |
|---|---|
|`odd_dyadic_valuation_contradiction`|The specified odd integer valuation cannot equal twice an integer|
|`unit_interval_cube_gap`|For0<m<1, m^2-m^3>0 from a positive product|
|`polynomial_factorization`|The cubic defined from S,A,B equals (x-w)(x-a)(x-b)|
|`midpoint_polynomial_positive`|The midpoint is in(0,1) and the polynomial there is strictly positive|
|`norm_A_large_contradiction`|Under S<=6, A>=2 contradicts the positive midpoint|
|`norm_B_large_contradiction`|Under S<=6, B>=2 contradicts the positive midpoint|
|`low_trace_contradiction`|With both norms1, S<=5 is impossible|
|`cubic_physical_root_gt_five`|All these branches imply w>5, including the remaining trace6 case|
|`cubic_literal_prefix_cost`|The cubic scaled reference norm equals40074077390|
|`cubic_prefix_reference_contradiction`|The three primitive cubic pressures contradict one energy budget|
|`quadratic_literal_prefix_cost`|The quadratic scaled reference norm equals40256297070|
|`quadratic_prefix_reference_contradiction`|The three primitive quadratic pressures contradict one energy budget|

The two auxiliary unprinted lemmas prove negativity of the whole
quadratics -5m^2+6m-2 and -5m^2+4m-1 using exact squares, rather than
sampled points. The midpoint factorization supplies positivity. Under
w<=5 the trace alternatives give S<=6. The two norm branches force
A=B=1; the low-trace branch forces S=6. The polynomial at5 is then
-1, contradicting its nonnegative factored value when w<=5. This proves
the full continuous cubic lemma without enumerating cubic fields.

For both reference lemmas use scale200000 and put T=200000*x.w.
The SINGLE primitive square has the form

    0<=40000000000*energy-400000*T+N.

The cubic literal vector is (123600,56006 five times,30189 ten times),
with N=40074077390 and

    T=67594*x1+25817*S6+30189*S16,
    400000*T-80148154780
      =27037600*(1000*x1-618)
       +103268*(100000*S6-201815)
       +1207560*(10000*S16-35276).

The quadratic literal vector is (122000,56326 five times,30837 ten times),
with N=40256297070 and

    T=65674*x1+25489*S6+30837*S16,
    400000*T-80512594140
      =262696000*(100*x1-61)
       +101956*(100000*S6-201815)
       +123348000*(100*S16-356).

Each coefficient is proved positive before the pressure terms are
combined. The literal squared norm is computed in Lean, not imported
as a final gain. The contradiction follows from that one norm square
and energy<1. The physical costs are1.00185193475 and1.00640742675.
No independently optimized norms or different farthest vectors are added.

## Execution and audit boundary

The original complete static review and root clearance bind predecessor
`b92f33dc...`, not the repaired source. Its actual failed compile and
all diagnostic output are byte-preserved in
`results/lean_b16_isolated_pair_scalar_verification_20261023/failed_compile_source_preservation.json`,
with verified gzip round trips for the predecessor source and map.
The failure produced no successful object or audit. Its diagnostic
`sorryAx` prints have zero verification credit. The successor changes
exactly two bound local names; no statement, hypothesis, tactic,
coefficient, mathematical input or requested declaration changes.

Exactly TWELVE footer `#print axioms` requests correspond to the twelve
table entries. A named compile and a SEPARATE actual audit are required;
successful older modules/dependencies must not be rerun. Actual admission
remains40GiB disk/6GiB RAM, one shared-lock process,2048MiB/one thread,
with the existing abort and cleanup behavior. A failed/refused attempt
must retain zero PASS credit and its exact source identity.

At this successor snapshot no successful compilation or separate actual
axiom audit is credited. Static review, literal
source hashes, mathematical acceptance and successful scientific output
are distinct from Lean compilation or actual axiom prints. Later execution
receipts may record a change in status without rewriting these historical
preparation bytes. Foundational axiom prints do not enumerate ordinary
theorem/structure premises; those are exposed explicitly above.
