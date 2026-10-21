# Complete compact Q(√5) dimension-27 coverage

The [binder source](../scripts/fresh_sqrt5_dimension27_coverage_20261015.py) and [complete certificate](fresh_sqrt5_dimension27_coverage_20261015.json) compose the frozen, independently reviewed arithmetic queues with direct exact face payments. Generation and saved replay pass. Its `--stage` mode remains a partial diagnostic that writes no data.

The certified conclusion is precise: every **compact arithmetic hyperbolic Coxeter polytope with ground field exactly Q(√5)** has dimension at most 27. The known compact dimension bound 29 reduces this claim to excluding dimensions 28 and 29. Other ground fields and noncompact chambers are outside this binder. The argument is a composition of reviewed face reductions and exact weights, not an enumeration of global polytopes or a claim of priority.

## A fixed rule on full vertex diagrams

Use the frozen [cross-weight rule](fresh_sqrt5_parameter_weight_stress_20261015_cross.md) and its [separate implementation review](fresh_sqrt5_parameter_weight_stress_20261015_cross_review.md). A pair at distance one or two in a spherical vertex component has base weight one. An entire H4 also receives a unit on its 5-edge and on its endpoint pair. The nonnegative connected-pair bonuses depend on the entire component and the number of entire A1 components in the full vertex, capped at two. One additional disconnected-pair bonus, 7/9, joins the unique entire A1 to a specified rooted E7 orbit, only when the full vertex contains exactly one A1.

The unique-A1 condition makes this last contribution additive over E7 components: each root in the specified orbit makes exactly one pair with that A1. There is no uncounted interaction between components. The separate review reconstructs all rooted-node multiplicities, all 92 finite component types through rank 29, the exact primal/dual certificate, and a different exact forest dynamic program. The binder recalculates the component totals and the full forest maxima from the literal rule.

The rule preserves the **weaker H4-only triangle baseline**. It does not require the optional extra unit on an I2(5) or H3 5-edge from the earlier triangle discussion. All 76 minimal Lannér triangle markings are checked directly with the H4-only baseline and receive at least two. Existing paths remain in any spherical extension, and an H4 cannot be a proper connected subdiagram of a larger spherical diagram. These are the extension-stability facts in the [independent triangle review](fresh_sqrt5_triangle_weights_20261015_review.md). Nonnegative new bonuses cannot reduce those payments.

## Exhaustive reduction of an unpaid quadrilateral

A marked distance-one or distance-two corner already pays one, persistently under further white extensions. Hence an unpaid quadrilateral satisfies the no-short-path premise of the [complete 5,440-support census](fresh_sqrt5_quad_exhaustion_20261015.md) and its [independent traversal](fresh_sqrt5_quad_exhaustion_20261015_incidence_review.md). It is unnecessary to impose every unrestricted 8,298-class union as a new weight constraint: the current rule has never decreased its short-path baseline.

For completeness, the geometric reduction uses the full compact polytope. The common white support W is spherical, and W plus either adjacent black pair is a finite corner. Each opposite-black extension is indefinite. A spherical extension would force the opposite facets to meet; an affine extension is excluded by compactness. Each chosen minimal indefinite piece contains both of its opposite black nodes and belongs to the complete arithmetic Lannér catalogue. The enumeration includes all shared whites and all finite cross contacts. It does not assume that the union is the entire face star.

The 1,180 marked union classes pass through the reviewed [principal arithmetic filters](fresh_sqrt5_quad_filters_20261015_review.md): 784 fail the field or signature conditions, leaving 362 finite-label cores and 34 dotted families. The [contact review](fresh_sqrt5_quad_contacts_20261015_review.md) then accounts for all 396 survivors:

| Disposition | Classes | Remaining obligation |
|---|---:|---|
| Persistent entire-H4 payment | 84 | Direct baseline check |
| Certified affine subdiagram | 2 | Compactness exclusion |
| Singular finite cores with no unpaid contact | 112 | Entire-component payments |
| Nonsingular finite cores with no unpaid contact | 44 | Full-star rank exclusion |
| Nonsingular finite cores with contacts | 130 | All 458 first extensions |
| Dotted families | 24 | Complete parameter and attachment queues |

“No unpaid contact” is conditional sealing. A first touching white root either creates a persistent short-path/H4 payment or occurs in the complete enumerated domain. If all unpaid profiles are excluded, every remaining white root is orthogonal to the core. Its marked components are then entire; arbitrary orthogonal white components remain possible. This is why every payment below is checked with all dimension-feasible spectator A1 counts. A sealed nonsingular core plus positive-definite orthogonal whites stays nonsingular, contrary to the full quadrilateral star's one-dimensional relation.

## The finite branches

The [complete first-extension containment census](fresh_sqrt5_sealed_subcore_closure_20261015_review.md) partitions the 458 accepted first extensions into 86 singular rows containing an original sealed subcore, 61 further singular rows, and 311 nonsingular rows. The binder reconstructs all 458 literal extensions and their partition. It checks every saved sealed-subcore map and directly evaluates **all 147 singular first-extension rows** with the current rule; it never assumes an old numerical payment transfers to changed weights.

The 61 singular rows have 51 marked representatives. Their [independently reviewed complete contact domains](fresh_sqrt5_extended_forest_weights_20261015_review.md) give 46 sealed states and five states with seven possible touching continuations. The binder checks all 61 representative maps and all 51 zero-touch payments under the current rule. The [17-state closure and review](fresh_sqrt5_singular_attachment_closure_20261015_review.md) close the seven touching continuations, together with sixteen auxiliary longer-chain starts. Six states contain known sealed subcores and eleven are newly sealed. All seventeen zero-touch rows, all seven continuation maps, and all saved subcore maps are checked again. Thus the earlier 174-row catalogue is used through explicit current-rule evaluations and structural sealing, not through an assertion that its old weight table is unchanged.

The remaining 311 nonsingular first extensions seed the [dimension-specific finite queue](fresh_sqrt5_nonsingular_rank29_exhaustion_20261015.md). A separately [accelerated exact continuation](../scripts/fresh_sqrt5_nonsingular_rank29_fast_20261015.py) starts from the immutable independently reviewed 1,433-state batch; the original implementation remains unchanged. Its final export and [separate complete independent review](fresh_sqrt5_nonsingular_rank29_fast_20261015_review.md) both pass generation and fresh saved replay. The terminal queue contains 5,172 states and 2,070 singular zero-touch obligations, all checked here. The earlier batch review is explicitly partial provenance, not a claim of exhaustion. Every singular state is a zero-touch obligation, including expanded states and states already paid under older rules. Every nonsingular no-touch state is rank-excluded. A touching state continues through an exact complete one-white domain unless it creates a persistent baseline payment or violates a necessary arithmetic condition.

The boundary of 31 star nodes is geometric for dimensions at most 29. A full n-dimensional quadrilateral star has n+2 roots and rank n+1. A 31-node nonsingular core cannot fit dimension 29; a singular core at that boundary is retained and its zero-touch payment is checked. This is not an experimental cutoff. The queue orders touching nodes by their first contact with the growing core; components orthogonal to it are spectators. Intermediate singular states are not omitted merely because they admit further touching roots.

## The dotted branches

The [independent parameter review](fresh_sqrt5_dotted_parameters_20261015_review.md) covers all 24 residual families: the exact physical upper bounds, conjugate interval and arithmetic integrality reduce them to 532 tested integer candidates and 32 admissible parameters. No floating-point distance sample is used. Their [independently reviewed dimension-29 queue](fresh_sqrt5_dotted_rank29_exhaustion_20261015_review.md) has 438 complete states and 264 singular zero-touch obligations: 250 sealed, four expanded and ten at the dimension boundary.

All 264 literal rows are evaluated directly under the current rule, including their dimension-feasible spectator cases. The dotted parameter is retained in each origin record. Although weights depend only on corner diagrams, parameter-specific arithmetic coverage is inherited from the exact parameter and queue certificates. No sealed-only subset is substituted for the full singular obligation list.

## Exact weight contradiction

In a zero-touch row with k star nodes, the displayed corner rank is k−2. At ambient dimension n, the external white rank is n−k+2. The binder checks every possible capped A1 count obtained by assigning some of that rank to A1 components and leaving either zero or at least two dimensions for non-A1 components. Thus the boundary cases do not assume that all three contexts are always realizable.

The exact maximum total vertex weights are:

| Dimension | No A1 | Exactly one A1 | At least two A1 | Overall |
|---|---:|---:|---:|---:|
| 28 | 230/3 | 230/3 | 230/3 | 230/3 |
| 29 | 239/3 | 84 | 236/3 | 84 |

With all quadrilateral cases paid, every k-gonal two-face receives at least 5−k: triangles receive two, quadrilaterals one, and larger faces need only nonnegative weights. Double-counting gives total weight at least `5 f₂ − binom(n,2) f₀`.

For a simple compact polytope, the symmetric nonnegative h-vector writes `f₂/f₀` as an average of the paired coefficients `(binom(n−i,2)+binom(i,2))/2`. Their minima are 91 at n=28 and 98 at n=29. The ratio is strictly greater than these minima because h₀=hₙ=1 and the extreme paired coefficient is larger. Hence the average vertex weight is strictly greater than 77 at dimension 28 and strictly greater than 84 at dimension 29. Both contradict the displayed upper bounds. This is the same published weighted-face/h-vector framework used in the [independent cubic-field proof review](fresh_compact_global_bound_20261015_review.md).

The complete certificate checks 471 literal maps and 2,661 zero-touch rows in 13,901 dimension-specific payment cases. It combines complete branch coverage, direct current-rule payments and this strict inequality. Generation and saved replay pass, with 60 recursive manifests and 477 source-hash references verified.
