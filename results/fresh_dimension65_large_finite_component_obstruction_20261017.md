# Dimension65 requires many finite vertices with large components

Let P be a finite-sided, finite-volume real hyperbolic Coxeter polytope
of dimension65. Let v_L count its finite vertices whose finite Coxeter
diagram has an irreducible component of rank at least18, let v_S count
its other finite vertices, and let v_i count its ideal vertices. Then

    7*v_L >= v_S + 5280 + (445/23)*v_i.

Equivalently,

    161*v_L >= 23*v_S + 121440 + 445*v_i.

In particular P must have at least755 finite vertices with a component
of rank at least18, or at least758 if it has any ideal vertex. More
generally, the lower bound is

    v_L >= ceil((23*v_S + 121440 + 445*v_i)/161).

Consequently, **if every finite vertex has all irreducible components
of rank at most17, the dimension is at most64**. This includes the
all-ideal case, where the hypothesis on finite vertices is vacuous.
The unrestricted dimension ceiling remains65. These are necessary
conditions, with no assertion that dimension65 examples exist.

The [binder](../scripts/fresh_dimension65_large_finite_component_obstruction_20261017.py)
and [certificate](fresh_dimension65_large_finite_component_obstruction_20261017.json)
compose completed geometric and exact counting theorems. They do not
run a new geometric census, finite knapsack calculation or ideal
partition enumeration.

## The compatible actual payment rule

Start with the completed
[dimension65 ideal-margin rule](fresh_ambient65_ideal_margin_closure_20261017.md).
Its ideal weights give every cusp absolute surplus at least445/23.
Apply the frozen
[finite rank16/17 local transfer](fresh_finite_d16_d17_long_pair_half_transfer_20261017.md):
the same-component costs become A16=116, B16=104, D16=118 and D17=125.
The old A17/B17 costs114/105 remain. Add only a half-coupon between any
mark in a whole A1 or A2 component and the unique long-arm neighbor of
the fork center of a separate whole D14 component. The added totals
are1/2 per A1/D14 pair of components and1 per A2/D14 pair.

The local theorem covers the complete affected selected support,
arbitrary supporting whites in the16 radical cases, all216 literal
replacement configurations, and the remaining universal interfaces.
Its whole-component sealing is what permits the typed finite coupon.
All current ideal weights, short8 and protected-leaf9 interfaces, and
every other coupon remain. Thus triangles still receive at least2,
quadrilaterals at least1, and every weight is nonnegative.

The [exact finite budget](fresh_d16_d17_d14_coupon_forest_budget_20261017.md)
counts precisely this rule, including the earlier E6/D11 quarter coupon.
It gives the following bounds for each finite vertex in dimension65:

|Finite vertex class|Total payment upper bound|Curvature coefficient|
|---|---:|---:|
|Arbitrary, hence also class L|487|480|
|Every component has rank at most17, class S|479|480|

The second bound is attained by the formal finite forest D17+3D16;
the unrestricted bound is attained by D65. These are exact abstract
forest maxima, not assertions that either local forest occurs in a
finite-volume dimension65 polytope. The two classes partition the
finite vertices; a class L vertex is counted once even if it has
several components of rank at least18.

## Global incidence comparison

Put D5=sum_F(5-number_of_sides(F)), summed over the two-faces of P.
For an ideal vertex v, write its actual Euclidean link as a product
of simplices of ranks r_vi summing to64. Set

    H_v(t)=product_i(1+t+...+t^r_vi),
    V_v=H_v(1),
    G_v=sum_j [t^j]H_v(t)*min(j,64-j),
    kappa_v=(5*G_v-64*V_v)/2.

Let C_v=V_v*sum_i c(r_vi) be the ideal payment envelope. The actual
factor-fiber multiplicity gives this upper bound also at nonsimple
cusps. The complete1,650-profile theorem, its all-rank concentration
argument, and the independent scalar reconstruction prove

    kappa_v - C_v >= 445/23

for every actual cusp. The finite local change leaves this comparison
intact. In particular it adds no ideal charge.

Complete triangle/quadrilateral coverage and nonnegative payments give

    D5 <= 487*v_L + 479*v_S + sum_(ideal v) C_v.

The strengthened shape-curvature theorem in the same frozen package gives

    D5 >= 480*(v_L+v_S) + sum_(ideal v) kappa_v + 5280.

Here480=(65-1)(65-5)/8 is the odd-dimensional finite coefficient.
The constant5280 is5*floor(65^2/4). It comes from the two distinct
global extrema: every simplex/product-simplex H coefficient is at
least1, so each extremum has E-G at least1056. It is independent of
the types of those extrema and is not counted a second time in the
per-vertex ideal surplus.

Subtracting the upper bound from the lower bound gives

    0 >= -7*v_L + v_S + sum_(ideal v)(kappa_v-C_v) + 5280,

which proves the displayed theorem. The integer conclusions follow by
rounding upward:5280/7 lies strictly between754 and755, while
(5280+445/23)/7 lies strictly between757 and758. Positive v_S or
additional ideal vertices only increase the required count.

If v_L=0, the right-hand side of the necessary inequality is strictly
positive and the left-hand side zero. Dimension65 is therefore
impossible for the stated class. The completed
[unrestricted ceiling65](fresh_ideal5_weighted_dimension65_20261017.md)
already excludes dimensions66 and higher, yielding the class ceiling64.
No assertion about closure of the class under faces or truncation is
needed. Arithmeticity, compactness, simple cusps and a facet-count
bound are not hypotheses.

This argument supplies no upper bound on v_L. The large finite forest
D65 still exceeds curvature by7, so the necessary inequality does not
exclude dimension65 in general. It sharpens the earlier finite-vertex
count by requiring all755 counted vertices to have a component of
rank at least18; it does not improve the unrestricted ceiling.

## Verification and inherited scopes

The binder checks the exact frozen input hashes, matching geometric
and budget rules, both finite maxima and their typed witnesses, the
unchanged ideal minimum, the strengthened constant, all coefficient
subtractions and integer thresholds. It hash-binds the complete scalar
rows to their
[independent reconstruction](fresh_ambient65_ideal_margin_closure_20261017_review.md).
No1,650-row arithmetic or finite-budget verification is rerun here.

The finite local packet completed generation, fresh replay and a full
root source/report/helper/interface review. After that packet was
frozen, another agent independently completed a full source/report
and native-label/induced-map-helper read, finding no gap in the complete
selected domains, arbitrary-white radical sealing, whole A1/A2--D14
roles, simultaneous old17/new16 debits and later rank guards. That
later read did not repeat maps, core arithmetic,216-case totals or
the upstream all-label census; it is recorded here without altering
the frozen local report. The exact finite budget has its separate
source/report review and full saved-table upper-and-attainment replay.
The shape, concentration and geometric reductions retain their stated
inherited review scopes.

Composition generation138aee and separate fresh saved replay776a09
both terminated exit0/PASS. Python AST, all16 direct bindings, all193
inherited input bindings and all7 local report links pass.

Root completed a separate full source/report composition review. It
checked the match between the frozen geometric rule and finite budget,
unchanged ideal charges, both D5 inequalities, the exact integer
coefficients and755/758 thresholds, and the class ceiling64 using the
unrestricted65 higher tail. It also checked that the extremum correction
is counted once, that no upper bound on v_L is inferred, and that no
closure under taking faces is assumed. That review did not independently
repeat the finite-budget or1,650-profile arithmetic. Source, certificate
and report are final/frozen. No external refereeing or proof-assistant
verification is claimed.

    python3 scripts/fresh_dimension65_large_finite_component_obstruction_20261017.py
    python3 scripts/fresh_dimension65_large_finite_component_obstruction_20261017.py --verify
