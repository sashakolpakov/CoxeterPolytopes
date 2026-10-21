# Unrestricted finite-volume dimension at most69

**Every finite-sided, finite-volume real hyperbolic Coxeter polytope has
dimension at most69.** The argument excludes dimension70 using two new
local projection obstructions and an exact cusp-profile comparison. The
completed [dimension70 theorem](fresh_ideal5_weighted_dimension70_20261017.md)
supplies all higher dimensions. Arithmeticity, compactness, simple cusps,
cusp-factor restrictions and a fixed number of facets are not required.

The [source](../scripts/fresh_ideal5_weighted_dimension69_20261017.py) and
[certificate](fresh_ideal5_weighted_dimension69_20261017.json) bind the local
theorem and all1,339 reduced cusp profiles. Verification status appears below.

## The two changed payments

Retain the entire preceding finite8/ideal5 rule, including all order13
transfers, order16/order17 half-weights and finite coupons. In affine B13
of order14 remove the singleton pair(0,10); in affine B14 of order15 remove
(1,11). The [complete local theorem](fresh_b14_b15_singleton_completion_transfer_20261017.md)
reduces their costs87 to86 and97 to96, with no new finite payment.

Every affected selected-core occurrence forces a missing terminal root
with a doubled bond to an old marked node. The terminal cannot be another
old corner root: the other old components are orthogonal to that marked
component, and the two omitted black roots are absent from the selected
corner. It is therefore a genuinely new supporting white root of square2.

The old core H is nonsingular and Lorentzian. Its complete necessary
contact domain splits into two ends with entrywise-positive inverse blocks.
Touching both ends supplies an already-paid short path. Otherwise the
forced doubled contact alone gives a projection lower bound2(H^-1)_ii.
All64 affected maps have this bound at least5/2, strictly greater than2.
Additional contacts on the same end only increase the projection. Thus
the actual affine occurrence is impossible in an otherwise unpaid face.
No unspecified full-matrix contacts are assigned zero in this argument.

There are23 affected order14 contexts with46 maps and11 order15 contexts
with18 maps; five contexts occur in both lists. Every elementary triangle,
short8 and basic ideal payment remains at unit weight. The inherited
low-core, disjoint-core, protected, immediate and finite alternatives are
unchanged for the precise domain reasons in the local theorem. Its complete
selected-map support exhausts the remaining changes. Consequently all
triangles still receive at least2 and all quadrilaterals at least1.

The [finite budget](fresh_finite8_d12_b10_coupon_budget_20261017.md) stays8n
through dimension75. The independently reconstructed13,034 budget endpoints
are inherited through the completed70 theorem. No new finite budget is
needed for the present dimension conclusion.

## Ideal costs and complete concentration

The order14 costs are now A70,B86,C84,D82, with maximum86. The order15
costs are A75,B96,C95,D96, with maximum96. Every other cost is unchanged.

|Affine rank r|Normalized factor cost c(r)|
|---|---:|
|1 through11|r/2|
|12|75/13|
|13|43/7|
|14|32/5|
|15|107/16|
|16|120/17|
|17 and larger|163/23|

At a cusp with V vertices in its Euclidean link, a factor of order r+1
multiplies its pair-weight sum by V/(r+1). Hence the total cusp cost is
at most V times the sum of c(r_i), including at nonsimple cusps.

For total cusp rank d=n-1, set

    H(t)=product_i(1+t+...+t^r_i), V=H(1),
    G=sum_j [t^j]H(t)*min(j,d-j).

The normalized curvature is(5G-dV)/(2V), and the ideal margin subtracts
the sum of c(r_i). The inherited convex-order lemma decreases this
curvature under concentration of ranks. The complete preceding reduction
domain remains valid. Low ranks concentrate up to11 with unchanged linear
cost. Two ranks12 become11 and13; the cost gain is now19/182, still positive.
Two ranks at least18 concentrate until one is17, preserving their constant
cost. Every nontrivial move increases the bounded integer sum of squared
parts, proving termination.

It suffices to keep arbitrary copies of11,13,14,15,16,17, at most one12,
at most one low rank1 through10 and at most one rank at least18. No further
rank-elimination rule is assumed. At total rank69 this domain contains
1,339 profiles, all with strictly positive exact ideal margin. The minimum
is1819/51840, attained at(11,11,11,11,11,14).

## Zero finite slack and the strict global contradiction

Let D5 sum5 minus the number of sides over all two-faces. Triangle and
quadrilateral coverage gives the upper bound

    D5 <= 8n*v_f + sum_(ideal v) V_v*sum_i c(r_vi).

The inherited [shape inequality](fresh_cusp_link_shape_curvature_20261017.md)
gives

    D5 >= beta_n*v_f + sum_(ideal v) (5G_v-dV_v)/2 + 5(n-1).

At n=70, beta_n=70*64/8=560=8n. Thus the finite coefficient has **zero
slack**, not negative slack. Each ideal coefficient has strictly positive
slack, and the additional global constant is345. Subtracting the upper
bound from the lower bound is therefore strictly positive even if there
are no ideal vertices. This contradiction excludes dimension70. The
completed ceiling70 excludes n>=71, proving the ceiling69.

At ambient69 the formal profile(11,11,11,11,11,13) still has negative margin
-2029/13608. This is a limitation of the present comparison, not a realized
cusp or a claim of optimality. The finite coefficient also needs stronger
control for a further reduction. No bound68 is asserted.

Compact arithmetic remains at most21 and arbitrary compact Coxeter
polytopes at most29. The effective compact-facet ambient ceiling is69;
the older auxiliary proof remains79. No new finite-volume construction or
new example dimension follows from this proof.

## Verification status

Primary generation and a fresh separate-process saved replay both terminated
exit0/PASS on all1,339 profiles and both changed local interfaces. The
[independent standard-library reconstruction](fresh_ideal5_weighted_dimension69_20261017_review.md)
agrees with every profile polynomial, curvature, cost and margin, the
minimum and the negative boundary profile. It also reconstructs both
canonical B targets, their entire saved affected-map support,13 distinct
old core inverses and signatures, actual preferred-corner contact domains,
positive end blocks and all64 forced terminal inequalities. Its own report
records its final generation and fresh replay.

A separate full internal source/mathematical review passed the primary
local composition, concentration moves, finite-budget equality, global
additive strictness and higher-dimensional tail. The frozen complete
core/map censuses and unchanged finite-budget reconstruction remain
inherited, with their stated scopes. All reviews are internal; no external
refereeing or proof-assistant formalization is claimed. Source, certificate
and report are final/frozen.
