# Unrestricted finite-volume dimension at most66

**Every finite-sided, finite-volume real hyperbolic Coxeter polytope has
dimension at most66.** The new proof excludes dimensions67 and68 by
reducing unnecessary long-pair payments in whole finite rank17 components.
The completed [dimension68 theorem](fresh_ideal5_weighted_dimension68_20261017.md)
excludes all higher dimensions. There is no arithmetic, compactness,
simple-cusp, cusp-factor or facet-count restriction.

The [source](../scripts/fresh_ideal5_weighted_dimension66_20261017.py) and
[certificate](fresh_ideal5_weighted_dimension66_20261017.json) assemble the
local theorem, exact finite budget and both complete ideal comparisons.
The final verification scope is recorded below.

## Reducing finite payments without changing ideal weights

The [complete local theorem](fresh_finite_d17_long_pair_support_20261017.md)
changes only same-component unit payments on whole finite A17,B17,D17.
Their costs become114,105,127, down from135,136,136. The old A17 endpoint
deletion remains in force. All finite pairs at distance at most8, and
all distance9 pairs with a marked leaf of the whole component, remain.
Every cross-component coupon and every ideal pair weight is unchanged.

These deletions have complete face provenance, not merely an absence in
a selected list. The entire2,614 high-core overlap domain splits into
2,391 cases with whole component bound16,74 protected cases,9 adaptive
cases and140 bound17 cases. Every induced finite rank17 map is enumerated
for the latter:48 A maps,18 B maps and168 D maps. Their marked pairs are
retained. Native labels are recovered from the scale-invariant quantity
4*g_ij^2/(g_ii*g_jj), so a rational root gauge cannot change the type.
Low-core, disjoint-core and universal short-path branches are separately
covered by their complete inherited classifications and proofs.

The local checker also reconstructs174 complete old and later literal
configurations and applies all three new rules simultaneously. Their
totals remain at least1. The sole debit is the already-removed A17
endpoint in old case662; both finite D17 alternatives remain, totaling2.
Thus the D17 change preserves the very payments used to justify the
earlier A17 deletion.

Universal transfer proofs are inspected at their actual conclusions.
Their unbounded-spectator alternatives use finite short paths, distance9
at a whole-component leaf, explicitly sealed small components, or unchanged
cross-component coupons. The long-cycle B/D leaf claims follow from
inverse-diagonal sealing; a possible extension of an A endpoint instead
forces the positive spectator correction already in that proof. Later
rank17 literal units are explicitly retained. A former permissive
"whole finite component at most17" phrase is never silently imported
after deleting one of its units.

Consequently triangles retain at least2 and quadrilaterals at least1.
There is no new coupon cost and no change in actual cusp multiplicity.

## Exact finite maxima for the actual new rule

The full variant of the [rank17 finite-forest calculation](fresh_rank17_reduced_forest_frontier_20261017.md)
assumes A17 base114, D17 base127 and B17 base at most127. The local
theorem now supplies those hypotheses, with B17 actually105. The other
finite types and all coupons are unchanged. The conditional budget is
therefore applicable to this actual face-payment rule.

The type reduction retains A17 separately because it has coupon roles;
non-coupon B17 is dominated by D17. The two A11 cycle coupons remain
included through the unchanged A11-to-D11 domination theorem. There is
no double counting of a half that was upgraded to one.

Enumerating every1,608 large-component count vector through rank75 leaves
675 additive Bellman tables. Each table gives exact upper bounds with
attaining predecessors, and every global combination is checked. The
[separate standard-library verifier](fresh_rank17_reduced_forest_frontier_20261017_review.md)
imports no producer arithmetic. For the full variant it independently
checks the type/coupon domination, every count vector,1,626,075 Bellman
upper alternatives,50,625 attaining predecessors,22,504 global upper
alternatives and all76 attaining forests including rank0.

At ambient67 the exact finite maximum is504, attained by the single
component D67. At68 it is512, attained by D68. These equal8n−32 and are
strictly below the corresponding curvature coefficients. The old4D17
obstruction has disappeared; its new cost is508 at rank68.

## The ideal rule is unchanged

For affine rank r the normalized factor cost remains

|r|c(r)|
|---|---:|
|1 through10|r/2|
|11|65/12|
|12|75/13|
|13|43/7|
|14|32/5|
|15|107/16|
|16|120/17|
|17 and larger|163/23|

At a cusp with link vertex count V, the
[actual factor-fiber theorem](fresh_cusp_factor_payment_framework_20261017.md)
multiplies a factor's pair cost by V/(r+1), even at nonsimple cusps.
For d=n−1 and factor ranks r_i, put

    H(t)=product_i(1+t+...+t^r_i), V=H(1),
    G=sum_j [t^j]H(t)*min(j,d-j).

The ideal margin is(5G−dV)/(2V) minus the sum of c(r_i).
The [all-rank convex-order lemma](fresh_distance14_unrestricted_dimension156_20261017.md)
and the completed68 concentration argument remain valid without change:
keep arbitrary copies of10,11,13,14,15,16,17; at most one12; at most one
part between1 and9; and at most one part at least18. The duplicate12 move
raises cost by23/1092. Low-rank concentration stops at10, not11. Every
move increases a bounded integer squared-rank sum and decreases the margin.

Both complete terminal domains are freshly enumerated and evaluated:

|Ambient dimension|Profiles|Minimum ideal margin|Minimizing ranks|
|---:|---:|---|---|
|67|1,761|851/165888|(11,11,11,11,11,11)|
|68|1,906|411499/3234816|(11,11,11,11,11,12)|

All3,667 margins are positive. The primary uses the frozen exact68
polynomial implementation; a separate global reviewer reconstructs the
domains and arithmetic independently.

## Two strict contradictions and the higher-dimensional tail

Let D5 sum5 minus the number of sides over all two-faces, and let v_f
count finite vertices. Face coverage gives

    D5 <= M_n*v_f + sum_(ideal v) V_v*sum_i c(r_vi).

The inherited [shape inequality](fresh_cusp_link_shape_curvature_20261017.md)
gives

    D5 >= beta_n*v_f + sum_(ideal v) (5G_v-dV_v)/2 + 5(n−1).

The exact coefficients are:

|n|M_n|beta_n|Finite slack|Global additive constant|
|---:|---:|---:|---:|---:|
|67|504|1023/2|15/2|330|
|68|512|527|15|335|

Every finite and ideal difference is positive, and the global constant
is positive. The two inequalities contradict one another in each
dimension, including when all vertices are finite or all are ideal.
The completed ceiling68 excludes n>=69, so no dimension at least67 is
possible. This proves the ceiling66.

At ambient66 this same comparison still fails in two separate ways:
D66 has finite cost496 against curvature495, and the formal cusp profile
(13,13,13,13,13) has ideal margin−1565/10976. These are comparison
obstructions, not realized polytopes or an optimality claim. No ceiling65
is asserted.

Compact arithmetic remains at most21 and arbitrary compact Coxeter
polytopes at most29. The effective compact-facet ambient ceiling is66;
the older auxiliary proof remains79. No finite-volume construction or
new example dimension follows from this proof.

## Verification scope

The local deletion theorem has completed generation, fresh replay and
two full internal mathematical/source reviews, including the universal
finite17 interfaces. Its complete inherited all-label censuses and
supporting-root sealing proofs retain their distinct verification scope.
The finite-budget calculation has a separate exhaustive certificate
checker. Primary global generation and a fresh separate-process replay
both terminated exit0/PASS. The
[independent global reconstruction](fresh_ideal5_weighted_dimension66_20261017_review.md)
agrees with every one of the3,667 polynomial, vertex-count, curvature,
cost and margin rows, using its own optional-first enumeration, sliding
convolution, cumulative-composition and absolute-deviation identities.
It independently re-enumerates all234 finite rank17 selected induced maps,
reconstructs the canonical paid sets, recomputes174 saved simultaneous
literal totals and directly counts the attaining D67/D68 pair costs.
The reviewer report records final bound generation and replay.

Separate full internal source/mathematical reviews passed the conditional
budget's application to the actual local rule, the distinct A17 coupon
role, unchanged A11 domination, cusp multiplicities and concentration,
both strict global contradictions, and the higher-dimensional tail.
The negative boundary66 profiles remain correctly distinguished from
actual realizations. No identified proof obligation remains in this
composition. Source, data and report are final/frozen.

No external refereeing, proof-assistant formalization or publication
priority is claimed.
