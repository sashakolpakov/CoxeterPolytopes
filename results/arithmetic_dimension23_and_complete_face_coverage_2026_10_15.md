# Arithmetic dimension bound and complete face coverage

**Internally reviewed computer-assisted proof candidate: every finite-volume
arithmetic hyperbolic Coxeter polytope has dimension at most 23.** This lowers
the published ceiling of 29 within the repository proof programme. The complete
face classifications, independent exact weight checks, and
[all-field composition](fresh_arithmetic_dimension_composition_20261015.md)
pass. External mathematical validation and literature priority remain open.

The completed finite search has **5,172 states, no pending cases and 2,070
singular stopping configurations**. It is an exhaustive local face-star
calculation through the already known dimension ceiling 29. A separate
implementation checks its arithmetic and all attachment domains. The dotted
families have their own complete, independently reviewed parameter and
attachment calculations. No new finite-volume polytope or example dimension
has been found.

## The dimension improvement

The nonexistence argument first concerns compact arithmetic polytopes over
Q(√5). The independently reviewed
[dimension-28/29 rule](fresh_sqrt5_parameter_weight_stress_20261015_cross_review.md)
[full-catalogue rules for 25–27](fresh_sqrt5_complete_face_weights_20261015_review.md),
and the [new dimension-24 rule](fresh_sqrt5_dimension24_i25_cross_20261015_review.md)
give the following exact bounds on the total angle weight at every vertex:

| Dimension | Vertex upper bound | Strict lower threshold from face counting |
|---|---:|---:|
| 24 | 267,527/5,000 = 53.5054 | 54 |
| 25 | 59 | 60 |
| 26 | 31,017/500 = 62.034 | 65 |
| 27 | 67,009/1,000 = 67.009 | 71.5 |
| 28 | 230/3 | 77 |
| 29 | 84 | 84 |

Every triangular face must receive at least two units and every quadrilateral
at least one. Larger faces require only nonnegative weights. Double-counting
and the symmetric nonnegative h-vector of a simple polytope force the average
vertex weight strictly above the final column. The strictness also excludes
dimension 29, where the two displayed numbers agree. The [complete coverage
composition](fresh_sqrt5_dimension23_coverage_20261015.md) establishes the
face inequalities across every branch of the classification, including
branches outside the optimization's original catalogue.

The new weights use the full spherical diagram at a vertex. In dimensions
25–27, connected-pair bonuses depend on the exact number a of isolated A1
components. A marked A1-to-component cross pair receives u_a/a; because each
other root forms exactly a such pairs, its total cross contribution is u_a.
This permits an exact additive enumeration of all vertex forests for fixed a.
All bonuses are nonnegative and preserve the earlier triangle and short-path
payments. Numerical optimization only proposes coefficients; rational checks
and a separate graph implementation certify them.

In dimension 24, the rule also counts entire I2(5) components. With `b` such
components, each eligible I2(5)-to-root pair receives `v_a/(2b)`. Its target
component excludes A1 and I2(5), so these payments do not overlap the A1
cross payments. Existing connected and A1 bonuses also distinguish whether
`b` is zero. Exact enumeration checks all 157 possible `(a,b)` contexts and
65,300 common spectator cases. The independent maximum is 53.5054.

The earlier A1-only rule has an exact model barrier of 55.5 at dimension 24.
Adding I2(5) crosses without the presence distinction still gives a barrier
above 54.19. These are obstructions to those specific weight models, not
existence theorems. The successful distinction crosses the required 54.

## Exhaustion and proof dependencies

The [finite queue](fresh_sqrt5_nonsingular_rank29_fast_20261015.md) starts from
all 311 nonsingular first extensions of the reviewed minimal quadrilateral
cores. Its 51,599 profile trials retain all singular stopping choices,
including expanded states. The 31-node boundary is forced by dimension:
a quadrilateral star in dimension n has n+2 roots. It is not a search-depth
guess. The [independent full-queue review](fresh_sqrt5_nonsingular_rank29_fast_20261015_review.md)
passes both generation and a fresh saved replay: all 5,172 states, 51,599
profiles and 10,085 seed/child maps agree. The primary terminal generation
and structural database/map checks also pass. The independently reviewed
initial batch remains preserved. A nonsingular sealed core is excluded by the full-star rank
condition; a singular one requires a direct weight payment.

The [dotted-parameter review](fresh_sqrt5_dotted_parameters_20261015_review.md)
reduces all 24 residual dotted families to 532 exact integral candidates and
32 seeds. Their [complete attachment review](fresh_sqrt5_dotted_rank29_exhaustion_20261015_review.md)
checks 438 states and 264 singular stopping choices. Distances are not sampled
or truncated numerically.

The final binder checks 2,661 zero-touch rows, of which 1,897 fit dimension 24,
in 60,735 common spectator cases. These include the earlier sealed cores, all singular first
extensions, the 51-state singular frontier and its 17-state closure directly
under the new weights. Orthogonal spectator components are allowed, with
every dimension-feasible A1 count and, for dimension 24, every I2(5) count. Structural sealing is retained; old
numerical payments are not assumed to survive changed weights.

Two earlier exploratory rules are not theorem inputs. The original all-A1
cross calculation omitted rank-two rooted contributions because of the
A2/I2(3) and B2/I2(4) naming aliases; the
[corrected certificate and erratum](fresh_sqrt5_all_a1_cross_corrected_20261015.md)
repair this, leaving its overall maximum unchanged. Subsequent partial
rules failed some late finite-queue cases. The complete-catalogue rules above
were reoptimized with all those cases and independently checked.

## Arithmetic scope

Vinberg's published theorem restricts compact arithmetic ground fields in
dimensions at least 22 to Q(√2), Q(√5), and Q(cos(2π/7)). The previous complete
repository candidates exclude dimensions above 20 over
[Q(√2)](fresh_sqrt2_quad_exhaustion_20261015_review.md) and
[the cubic field](fresh_compact_global_bound_20261015_review.md).
Consequently the completed Q(√5) bound of 23 gives the same bound for all
compact arithmetic reflection groups within this proof candidate. The
[composition certificate](fresh_arithmetic_dimension_composition_20261015_bound23.json)
checks all eighteen field/dimension branches for dimensions 24–29. The noncompact arithmetic case already
has the published sharp ceiling 21, with dimension 20 excluded.
[Vinberg, Theorems 1–2](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt),
[Esselmann](https://doi.org/10.1006/jnth.1996.0141).

Every arithmetic hyperbolic lattice generated by reflections is of the
simplest type. Thus this is already the arithmetic first-type reflection
problem; it is not a bound on arbitrary first-type arithmetic lattices.
Vinberg's root enumeration applies to the number-field lattice setting.
Allcock's published edgewalking algorithm is specifically formulated for
Z-lattices and requires a starting vertex; its geometric idea alone does not
provide a complete number-field or nonarithmetic search algorithm.
[Vinberg's lemma as stated in Belolipetsky, Lemma 2.1](https://arxiv.org/pdf/1506.03111),
[Allcock, Introduction](https://web.ma.utexas.edu/users/allcock/research/edgewalking3.pdf).

The unrestricted compact ceiling remains the published 29. The separate
unrestricted finite-volume **163 proof candidate** is unchanged in the
[preceding checkpoint](all_quadrilateral_domains_and_two_endpoint_bound_2026_10_15.md).
Internal exact verification, external mathematical validation and establishing
novelty in the literature are distinct; no priority claim is made here.

## Next research step

The same hub-presence model has an [exact restricted-model obstruction in
dimension 23](fresh_sqrt5_dimension23_i25_cross_20261015.md): its vertex bound
cannot reach the required 49.5 under the supplied conservative face constraints.
This does not prove that dimension 23 is attainable or that a stronger bound
is impossible. Further progress needs richer whole-diagram contexts, different
weights, or geometric exclusions of the obstructing principal configurations.
The constructive search for a new finite-volume example remains open.

## Constructive lane and archive

A [new E8-based arithmetic H9 lattice prefix](fresh_silver_e8_h9_20261015.md)
has 20 exact walls and a strict infinite-volume escape. It remains an open
host; this prefix does not prove nonreflectivity or give a new example.
The older H9/H20 geometric constructions and excluded lattice hosts retain
their preceding checkpoint status.

The [archive checker](../scripts/audit_arithmetic_dimension23_checkpoint_20261015.py)
and [integrity certificate](arithmetic_dimension23_and_complete_face_coverage_2026_10_15_audit.json)
bind the new results and preserve the preceding archive. Two historical
batch snapshots retain explicitly nonbinding active-database metadata; their
immutable bytes and exact certificates remain bound. Current terminal
databases have no such exemption. Archive integrity does not replace the
mathematical reviews.
