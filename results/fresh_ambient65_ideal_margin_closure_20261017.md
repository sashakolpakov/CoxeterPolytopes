# All cusp margins at dimension65 are positive

The proved local improvements now resolve the entire ideal-vertex
comparison at ambient dimension65. Every cusp has normalized surplus
at least61/1936 and absolute surplus at least445/23. Consequently every
finite-sided, finite-volume real hyperbolic Coxeter polytope whose vertices
are all ideal has dimension at most64. The unrestricted ceiling remains65:
its finite-vertex budget still exceeds curvature by7.

The [source](../scripts/fresh_ambient65_ideal_margin_closure_20261017.py)
and [certificate](fresh_ambient65_ideal_margin_closure_20261017.json)
check all1,650 terminal profiles. This replaces the earlier
[conditional ideal targets](fresh_order13_b73_scalar_frontier_20261017.md)
by compatible proved local rules and an exact budget for the new coupon.

## Compatible local rules and complete profile domain

The order15
[B](fresh_order15_affine_b_residual_terminal_transfer_20261017.md),
[C](fresh_order15_affine_c_terminal_projection_transfer_20261017.md) and
[D](fresh_order15_affine_d_half_coupon_transfer_20261017.md) refinements
give costs A75/B94/C91/D94 without added finite cost. The
[order17 D transfer](fresh_order17_affine_d_three_quarter_transfer_20261017.md)
gives A85/B111/C108/D119 after adding the precisely typed E6/D11 quarter
coupon. Its [complete finite budget](fresh_e6_d11_quarter_coupon_forest_budget_20261017.md)
leaves every maximum through75 unchanged. These changes affect distinct
ideal orders, preserve their short and triangle interfaces, and introduce
only that one new nonnegative finite coupon. The earlier order13 maximum74
and all other payments remain in force.

For affine rank r, the resulting normalized envelopes are:

|Rank|Cost c(r)|
|---|---:|
|1 through10|r/2|
|11|16/3|
|12|74/13|
|13|6|
|14|94/15|
|15|107/16|
|16|7|
|17 and larger|163/23|

Use the inherited all-rank convex-order concentration theorem. Concentrate
two ranks below10 until one disappears or the other reaches10. Concentrate
two ranks at least18 until one reaches17. These moves decrease normalized
curvature and preserve total cost, because the respective cost ranges are
linear and constant. Each increases the bounded sum of squared ranks.

The complete terminal domain therefore has arbitrary multiplicities of
10 through17, at most one1 through9 and at most one rank at least18.
There is no extra duplicate-rank reduction. All1,650 terminal profiles of
total rank64 are enumerated.

For each profile define H(t)=product_i(1+t+...+t^r_i), V=H(1), and
G=sum_j[t^j]H(t)*min(j,64-j). The exact normalized margin is

    m = (5G-64V)/(2V) - sum_i c(r_i).

Its minimum is61/1936 at(9,10,10,10,10,15). Every margin is positive.
All polynomials, V and G agree with the frozen complete-domain predecessor,
and all new margins agree with its formerly conditional simultaneous target.

## Absolute surplus and a stronger global constant

The absolute surplus V*m has minimum445/23, attained by the single rank64
profile. This lower bound also applies before concentration. Indeed, for
0<a<=b and a positive allowed transfer delta,

    (a+1)(b+1)-(a-delta+1)(b+delta+1)
      = delta*(b-a+delta) > 0.

Thus V decreases under every concentration, with a vanished factor having
value1. Normalized margin also decreases, and its terminal value is
positive. Their product consequently decreases. The minimum over the
terminal domain bounds the absolute surplus for every actual cusp.

The [shape-curvature theorem](fresh_cusp_link_shape_curvature_20261017.md)
has a quantitative extremum argument that admits a sharper constant here,
and in fact for all the same finite-volume Coxeter polytopes. At every
vertex figure, finite or ideal, H is a product of polynomials
1+t+...+t^r, so every coefficient h_j for0<=j<=d is at least1. Its proof gives

    E-G = sum_(j>d/2) h_j*(2j-d)
        >= sum_(j>d/2)(2j-d)
        = floor((d+1)^2/4).

At each of the two distinct global extrema the cutting hyperplane misses
the entire vertex figure, giving E rather than G uncut edges. Their
combined extra contribution to D5 is therefore at least
5*floor(n^2/4), strengthening the earlier5(n-1). No new geometric
hypothesis is needed. At n=65 this constant is5280.

## The resulting theorem and remaining finite obstruction

Write v_f and v_i for the numbers of finite and ideal vertices. The
compatible triangle/quadrilateral payment theorem gives

    D5 <= 487*v_f + sum_(ideal v) V_v*sum_i c(r_vi).

The strengthened shape inequality gives

    D5 >= 480*v_f + sum_(ideal v) (5G_v-64V_v)/2 + 5280.

Combining them with the absolute ideal surplus proves the necessary
dimension65 inequality

    7*v_f >= 5280 + (445/23)*v_i.

In particular, v_f=0 is impossible. The completed
[unrestricted ceiling65](fresh_ideal5_weighted_dimension65_20261017.md)
already excludes dimensions66 and higher, so the all-ideal ceiling is64.
The same inequality requires at least755 finite vertices in any remaining
dimension65 example, or at least758 if it has an ideal vertex. These are
necessary conditions, not claims that such examples exist.

There is no upper bound on v_f in this argument. In particular the
inequality does not eliminate dimension65 in general. The formal finite
forest D17+3D16 still costs487 against480. Arithmeticity, simplicity of
cusps and a facet-count bound are not assumptions of the all-ideal result.
No new finite-volume example or arithmetic dimension bound follows.

## Verification status

The local input theorems and the complete finite budget have finished
generation, fresh replay and separate internal mathematical reviews.
Final generationa94cea and separate saved replay967ad6 terminated exit0/PASS
on all1,650 rows. The
[independent scalar reconstruction](fresh_ambient65_ideal_margin_closure_20261017_review.md)
uses a different enumeration and arithmetic implementation and agrees on
every saved scalar field and both unique extrema. Its separate replay passes.

A separate full source/report mathematical composition review reread the
shape-curvature and factor-fiber proofs, final finite budget and independent
scalar scope. It checked the5280 extremum constant, product monotonicity,
absolute445/23 surplus, simultaneous payment rules and the all-ideal64
conclusion, finding no unresolved step. It did not duplicate the scalar
arithmetic, geometric censuses or finite-budget replay. The geometric
census and all-rank convex-order theorem retain their explicitly stated
inherited review scopes. No external review or proof-assistant
formalization is claimed.

    python3 scripts/fresh_ambient65_ideal_margin_closure_20261017.py
    python3 scripts/fresh_ambient65_ideal_margin_closure_20261017.py --verify
