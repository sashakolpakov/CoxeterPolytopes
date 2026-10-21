# Independent scalar and composition review of dimension71

**Completed independent exact reconstruction; fresh replay passed.** This review supports the [unrestricted dimension71 theorem](fresh_ideal5_weighted_dimension71_20261017.md). It imports no producer or project arithmetic code. Its [standard-library source](../scripts/fresh_ideal5_weighted_dimension71_20261017_review.py) and [saved certificate](fresh_ideal5_weighted_dimension71_20261017_review.json) independently reconstruct all4,913 scalar profiles, all13,034 finite-budget endpoint cases, and the changed order13/order17 composition interfaces.

The primary's local geometric theorems have also received a full source and mathematical read in the scopes described below. The inherited complete core, profile and induced-embedding censuses are frozen dependencies; this review does not claim to independently repeat those censuses or every local Gram calculation.

## Independent scalar calculation

The normalized factor-cost envelope is

|Rank|Cost|
|---|---:|
|1 through11|r/2|
|12|75/13|
|13|87/14|
|14|97/15|
|15|111/16|
|16|120/17|
|17 and above|163/23|

The independent enumeration chooses the optional high rank, optional low rank and optional rank12 first, then recursively fills the remainder with11,13,14,15,16,17. The producer enumerates the repeated parts first and obtains the optional parts from its remainder. Exact profile-set equality verifies both completeness and absence of duplicates.

Low parts concentrate toward11 without changing their linear cost. Two12 parts become11 and13, increasing cost by16/91. Two parts at least18 concentrate until one reaches17, preserving their common constant cost. These moves decrease normalized curvature by the inherited convex-order result. Each nontrivial move strictly increases the integer sum of squared parts, bounded by the square of their fixed total. Thus termination and the stated reduced domain follow. Repeated15 parts are retained: their old elimination has the wrong cost direction after the rank16 reduction.

For each profile, sliding-window convolution constructs the coefficients of `H(t)=Π(1+...+t^r)`. Independently, the bounded-composition numerator `Π(1−t^(r+1))` gives a cumulative-count identity

`G=floor(d/2)V−2 Σ_s a_s binom(floor(d/2)−s+k,k+1)`.

This agrees exactly with `Σ_j h_j min(j,d−j)`. The curvature `(5G−dV)/(2V)` also agrees with the absolute-deviation formula `3d/4−(5/2) E|Σ U_i−d/2|`. Integer convolution and rational comparisons therefore have two distinct exact reconstructions in addition to comparison with the primary output.

|Ambient dimension|Profiles|Minimum ideal margin|Minimizing ranks|Finite slack|
|---|---:|---:|---|---:|
|72|1530|6875/49152|11,15,15,15,15|18|
|73|1635|41401/145152|11,11,11,11,13,15|28|
|74|1748|135455/304128|10,11,11,11,15,15|37|

Every primary profile row matches, including its full h-polynomial, vertex count, G, normalized curvature, cost and margin. The boundary profile `(11,11,11,11,11,15)` at ambient71 has margin `−5597/248832`, also matching independently. It is a scalar limitation of this envelope, not a realized cusp or polytope.

The unchanged all-rank B/C/D feature bounds are independently evaluated at every piecewise-linear breakpoint and their limits. On each branch the normalized cost has form `a+b/m`, so endpoints and the limiting value suffice. Their maxima are160/23,75/11 and163/23, respectively. This is an all-rank branch argument, not a finite numerical cutoff.

## Finite budget and local interfaces

The [new finite8 budget](fresh_finite8_d12_b10_coupon_budget_20261017.md) applies through ambient75, covering all three new dimensions. The review independently enumerates counts of D12,D13,other inherited large components,A17 andB10 first, and only then all admissible dimensions and the two small-rank endpoints. The resulting13,034 exact rows equal the producer's rows. Every surplus is nonpositive. This domain ordering differs from the producer's dimension-first loop splitting a total large-component count into its D12/D13 subcounts.

The component slacks and disjoint coupon increments also pass source/mathematical review: the D12 short-leaf upgrade costs at most `a*b`, the A2/B10 halves cost at most `a*v/2`, and the inherited A4/D13-neighbor charge remains `a*c/4`. Their rank constraint and affine endpoint reduction are valid. No geometric transfer follows merely from that numerical budget; those transfers are checked separately.

The review reconstructs the four canonical order13 affine diagrams and their distances without project graph helpers. It checks all elementary triangle pairs, changed short8 pairs, selected marked-image supports and immediate-factor interfaces against the primary. The resulting costs are A65,B75,C73,D74.

- The [cycle13 theorem](fresh_cycle13_ideal5_transfer_20261017.md) replaces all13 distance-six pairs, including their use in the original short8-first stage. Its fixed simple leaf and degree-two neighbor force the single-contact corner component to be finite even with arbitrary spectators. The remaining profiles force both isolated heights to1 before the unscaled energy10/13 is used. The strict spectrum excludes that value. The full source/proof was read; the six/16 arithmetic cases were not separately rerun here.
- The [B12 singleton theorem](fresh_order13_affine_b_singleton_transfer_20261017.md) seals the selected radical cores and supplies whole finite coupon components. The four saved conditional distance11 cases are checked separately: two have a finite B13 unit, and two receive the now-budgeted whole A1/D12-short-leaf unit. Its strict universal inner4 gaps1/5 and1/6 discharge the separate immediate-D14 provenance. Full source/proof review found no missing contact or indirect-tail assumption; the underlying801-pattern census remains inherited.
- The [universal D12 leaf theorem](fresh_affine_d12_leaf_pair_transfer_20261017.md) covers every use of its four removed pairs. At its zero boundary, the positive spectator quadratic form and disjoint supports force both contact vectors to vanish. With one-sided contact, choosing the corner that excludes the contacting root avoids any spectator classification. With two-sided contact, equal heights are established before the strict energy-one theorem is applied; its short-path or two-whole-A1 outcomes give the stated payments. Whole actual components and strict edge hypotheses were checked in the source/proof read; the96 literal configurations were not rerun here.
- The [D16 weighted theorem](fresh_d17_partial_weight_transfer_20261017.md) changes twelve order17 pairs to weight1/2. Their complete saved support consists of twelve sealed cores, each with an existing half at a distinct finite corner. The review checks the full136-entry weight vector, its exact sum120, unchanged triangle/short8 unit weights, and all twelve saved payment totals. The former zero-weight singleton retains its prior proof. A half-weight is never treated as a binary unit.

The canonical graph checks and exact saved interfaces are independent computations. The associated full source/mathematical reads are separate from the geometric censuses they bind.

## Global composition and scope

The finite cost8n, ideal factor-fiber multiplicity `V/(r+1)`, and fractional pair weights compose linearly. All triangle totals remain at least2 and quadrilateral totals at least1. The new local replacements cover the formerly paid order13 short, selected and immediate cases; unchanged cases retain the earlier complete face theorem.

Every reduced ideal profile at ambient72–74 has positive margin, and the finite slacks are18,28,37. The inherited global shape inequality has the positive constant `5(n−1)`. Summing the face payments therefore contradicts the vertex budgets in these three dimensions. The completed [dimension74 theorem](fresh_ideal5_piecewise_dimension74_20261017.md) excludes every dimension at least75. This proves the unrestricted ceiling71 under the stated finite-sided, finite-volume real-hyperbolic Coxeter hypotheses, without arithmeticity, simple-cusp or factor-count restrictions.

The negative boundary margin does not prove optimality of71. No new example, compact-arithmetic dimension result, external refereeing or proof-assistant formalization is asserted. The compact arithmetic partial models remain separate from this geometric proof.

## Verification status

All4,913 scalar rows,13,034 budget rows and composition interfaces passed
comparison with the primary. The final generated certificate binds the
frozen primary source, data and report. A fresh separate-process saved
replay, run directly by the root agent as session15598, terminated
exit0/PASS on the complete reconstruction. All33 direct source bindings,
Python syntax and local report links pass.

A separate full internal source/mathematical read passed the independent
partition enumeration, three curvature calculations, all-rank branch bounds,
alternative finite-budget enumeration, canonical graph reconstruction and
weighted local composition. That read did not duplicate the4,913 scalar
calculations or inherited geometric censuses. Source, data and report are
final/frozen. All reviews are internal; no external refereeing or
proof-assistant formalization is claimed.
