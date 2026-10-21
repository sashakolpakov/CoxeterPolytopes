# Dimension ten: proof, assumptions and reproduction

20 October 2026. Current result:
[finite-volume right-angled dimension at most10](../results/right_angled_mixed_dimension10_20261020.md).
The new argument excludes mixed dimension11 and uses face heredity for
the higher tail. It is unconditional as a mathematical theorem and has
conditional, rather than end-to-end geometric, Lean verification.
The consolidated [TeX](../manuscripts/right_angled_dimension10.tex) and
[12-page PDF](../manuscripts/build/right_angled_dimension10.pdf) are
complete: three visually checked vector figures, zero warnings and
two full manuscript reviews. The older11-page ceiling11 paper is
preserved as a predecessor.

|Step|Written source|Verification boundary|
|---|---|---|
|Actual P3 recognition, fixed componentwise null direction|[Assembly](../results/right_angled_mixed_dimension10_assembly_candidate_20261020.md), geometry shared with the [predecessor](../manuscripts/right_angled_dimension11.tex)|All-wall feasibility, actual finite/cusp links and orientation are written geometry.|
|Strict weighted row and coface/Euler identities|[Weighted-face proof, dimension11 section](../results/right_angled_mixed_weighted_triangles_research_20261018.md)|DS, standard truncation/links and published82 facets in dimension8. Earlier exact DS checker has fresh bounded replay.|
|One-cusp exclusion, R7 balance and exceptional budget|[Local proof](../results/right_angled_one_cusp_ten_facet_exclusion_research_20261020.md), [balanced counts](../results/right_angled_mixed_eleven_three_face_corner_balance_20261020.md)|Complete own facet/edge argument and full root review; no 3D classification imported.|
|H<=2E without D3 recognition|[Full assembly review](../results/review_right_angled_mixed_dimension10_assembly_finite57_20261020.md)|Direct three-unit excess, Euler and integer cusp capacity.|
|Structured scalar contract, k<=6|[Moment source, Section7](../results/right_angled_mixed_k10_subset_moment_reduction_20261020.md)|Elementary Cauchy and complementary-triple convexity.|
|k7,8,11|[Moment and rational-interval proof](../results/right_angled_mixed_eleven_structured_small_and_full_subsets_20261020.md)|Every interval reconstructed with exact rational square-root rounding; all11 case is a single rational comparison.|
|k9|[Global analytic proof](../results/right_angled_mixed_eleven_nine_direction_global_proof_20261020.md)|Full human inequality; exact polynomial identity replay. Independent earlier interval proof also passes.|
|k10|[Moment reduction](../results/right_angled_mixed_k10_subset_moment_reduction_20261020.md), [rational envelopes](../results/right_angled_mixed_k10_rational_envelopes_20261020.md)|Complete corrected moment proof; all12550 Bernstein coefficients positive on840 closed intervals.|
|Deletion, ordinary reciprocal summation and contradiction|[Assembly](../results/right_angled_mixed_dimension10_assembly_candidate_20261020.md), [review](../results/review_right_angled_mixed_dimension10_assembly_finite57_20261020.md)|One fixed null vector per actual component, paired reverse edges; no maximum-height substitute.|
|Dimension tail|Same full review|Eleven-faces inherit finite volume and right angles; published compact4 and ideal6 force them mixed.|

## Reproduce exact scalar checks

The catalogue provides one command for the four declared exact checks,
the conditional kernel/axiom audit and a fresh temporary PDF build:

```sh
python3 -B scripts/verify_manuscript.py mixed10 --mode core --pdf
```

Each job is serial and individually resource capped. A successful run
has the stated conditional scope; it does not formalize actual geometry.
The [actual six-job run](../results/mixed10_core_verification_20261020.json)
passes, including the separate axiom audit and a fresh temporary PDF build.

Run the following commands sequentially. Each uses the shared lock and
the59-second/1GiB process-tree watchdog; the normal40GiB/6GiB start
gates apply. A replay writes no certificate and launches no optimizer.

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- scripts/check_mixed_eleven_balanced_and_subset_rows_20261020.py --verify
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- scripts/check_mixed_eleven_nine_direction_global_20261020.py --verify
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- scripts/check_mixed_k10_rational_envelopes_20261020.py --verify
```

The six matching generation/replay ledgers begin
`mixed_eleven_balanced_subset_`, `mixed_eleven_nine_direction_global_`
and `mixed_eleven_k10_envelope_` in `results/`. The optional earlier
nine-direction interval checker is
`scripts/check_mixed_eleven_nine_direction_subset_20261020.py`.

## Conditional kernel

[RightAngledMixed10.lean](../formal/CoxeterProofs/RightAngledMixed10.lean)
imports only Std. It checks the count logic from typed integer identities
and local incidence/error bounds. The cleared aggregate real-analytic
transport comparison is an explicit structure field; the module does
not pretend that individual real ratio rows are integers.

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.RightAngledMixed10 --audit-axioms
```

Actual compilation and a separate eight-theorem axiom audit pass.
Only propext and Quot.sound occur, with no sorryAx; the integer
comparison has no axioms. Ledgers are
`formal/.proof-batches/CoxeterProofs_RightAngledMixed10.json` and
`formal/.proof-batches/CoxeterProofs_RightAngledMixed10_axiom_audit.json`.
This verifies the explicit logic and its assumptions, not universal
hyperbolic geometry. The script results, written geometry and typesetting
retain separate verification records.
