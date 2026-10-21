# Unrestricted finite-volume dimension at most70

**Every finite-sided, finite-volume real hyperbolic Coxeter polytope has
dimension at most70.** The proof uses the new weighted order16 local theorem
and an exact comparison in dimension71. The completed
[dimension71 theorem](fresh_ideal5_weighted_dimension71_20261017.md) supplies
all higher dimensions. There is no arithmeticity, compactness, simple-cusp,
cusp-factor or facet-count restriction beyond finite sidedness.

The [source](../scripts/fresh_ideal5_weighted_dimension70_20261017.py) and
[certificate](fresh_ideal5_weighted_dimension70_20261017.json) bind the complete
local transfer and compare all1,432 reduced cusp profiles. Verification and
review status are recorded below.

## The changed local rule

Retain the entire preceding finite8/ideal5 face rule, including the order13
transfers, finite coupons and twelve half-weights in order17. The only
change is in affine D15, a factor with16 nodes. Its ten still-paid pairs
at distance greater than9 now have weight1/2. The
[complete local proof](fresh_d16_partial_weight_transfer_20261017.md)
reduces its total weight from111 to106 without adding any finite coupon.

Every elementary triangle payment and every short8 unit remains unchanged.
The low-core, disjoint-core, protected-feature, immediate-factor and adaptive
finite alternatives retain their old proofs, for the explicit domain reasons
given in the local report. The complete induced-map census contains54
affected selected contexts and328 affected maps. The local proof controls
actual affine-factor extensions, including arbitrary contacts and indirect
tails, before fixing their completed Gram matrices. Of66 distinct literal
completions,20 require multiple core-contact roots,14 require projection
square greater than2, and22 have two negative directions. These cannot
occur. The ten remaining factors are already whole and their cores are
sealed against all further supporting white roots.

Each surviving configuration has a different actual finite corner carrying
an existing half-payment: whole A4 against a whole D12 fork center or short
leaf, or whole D5 against a whole D11 short leaf. These are whole actual
components, including in the presence of arbitrarily many spectators. The
finite half and the remaining ideal half total1. Thus all triangles still
receive at least2 and all quadrilaterals at least1.

The unchanged [finite budget](fresh_finite8_d12_b10_coupon_budget_20261017.md)
is8n through dimension75, including dimension71. Its13,034 exact endpoint
cases and mathematical composition were independently reconstructed in the
[preceding review](fresh_ideal5_weighted_dimension71_20261017_review.md).
No new finite expense is added here.

## The ideal-factor envelope

The order16 costs are now A80, B107, C103 and D106. Hence the maximum is107,
and the rank15 coefficient falls from111/16 to107/16. All other ranks retain
their preceding costs:

|Affine rank r|Normalized cost c(r)|
|---|---:|
|1 through11|r/2|
|12|75/13|
|13|87/14|
|14|97/15|
|15|107/16|
|16|120/17|
|17 and larger|163/23|

For an actual cusp with factor ranks r_i and V vertices in its Euclidean
link, each pair in a factor of order r_i+1 occurs V/(r_i+1) times. Fractional
weights respect this exact counting formula. Therefore the cusp cost is
at most V times the sum of c(r_i). No simple-cusp assumption is used.

## Complete reduction and exact comparison

For a cusp of total rank d=n-1, define

    H(t)=product_i(1+t+...+t^r_i), V=H(1),
    G=sum_j [t^j]H(t)*min(j,d-j),
    normalized curvature=(5G-dV)/(2V).

The ideal margin is this curvature minus the sum of c(r_i). The preceding
convex-order lemma decreases curvature under concentration of two ranks.
We retain exactly its complete reduction domain. Low ranks concentrate up
to11 without changing their linear cost, leaving at most one rank1 through10.
Two12 ranks become11 and13, increasing cost by16/91. Two ranks at least18
concentrate until one is17, preserving their common constant cost. These
moves do not involve the changed rank15 cost. Each move strictly increases
the bounded integer sum of squared parts, so the procedure terminates.

Every partition consequently reduces to repeated11,13,14,15,16,17, at most
one12, at most one rank1 through10, and at most one rank at least18. This
complete domain need not be minimal. In particular repeated15 parts are
retained and no additional elimination rule is assumed.

At ambient dimension71 the total rank is70. The complete reduced domain has
1,432 profiles. Every exact ideal margin is positive; the minimum is

    23641/338688,

attained at(11,11,11,11,13,13). The finite-vertex slack is

    beta_71 - 8*71 = (70*66)/8 - 568 = 19/2.

To finish, let D5 be the sum of5 minus the number of sides over all two-faces.
If T and Q count triangular and quadrilateral faces, then

    D5 <= 2T+Q <= total corner payments
       <= 8n*v_f + sum_(ideal v) V_v*sum_i c(r_vi).

The inherited [shape inequality](fresh_cusp_link_shape_curvature_20261017.md)
gives

    D5 >= beta_n*v_f + sum_(ideal v) (5G_v-dV_v)/2 + 5(n-1).

At n=71, every finite and ideal coefficient on the lower-bound side exceeds
the respective upper-bound coefficient, and the additive constant is350.
This is impossible. The completed ceiling71 excludes every n>=72; together
these arguments establish the ceiling70.

At ambient70, the formal profile(11,11,11,11,11,14) still has negative margin
-1637/51840. It contains no rank15 factor, so further improvement of the
order16 rule alone does not remove this obstruction. This is a limitation
of the present comparison, not a realized cusp, a polytope or an optimality
claim for70. No bound69 is asserted.

Compact arithmetic remains at most21, and arbitrary compact Coxeter
polytopes remain at most29. The effective compact-facet ambient ceiling
becomes70; its older auxiliary proof remains79. No new finite-volume example
or new construction dimension has been established.

## Verification status

Primary generation and a fresh separate-process saved replay both terminated
exit0/PASS on all1,432 comparisons and the saved local interfaces. The
[independent standard-library reconstruction](fresh_ideal5_weighted_dimension70_20261017_review.md)
agrees with every profile, its full polynomial and curvature, the minimum,
the negative boundary profile and the changed local interface. It also
independently checks all328 saved affected maps, the nine distinct full Grams
underlying the nontrivial completions, and all40 sealed corners and their
whole finite half-payments. Its own report records its generation and replay.

A separate full internal source/mathematical read passed the new local
composition, the unchanged complete concentration domain, the finite budget,
the higher-dimensional tail and the global contradiction. The local
geometric census is an explicit frozen input, not independently enumerated
by this global checker or the new reviewer. The unchanged13,034 finite-budget
endpoint reconstruction remains inherited from the frozen71 review.
All reviews are internal; external mathematical refereeing and
proof-assistant formalization are absent. Source, certificate and report
are final/frozen.
