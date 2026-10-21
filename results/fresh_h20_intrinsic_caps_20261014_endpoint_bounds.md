# Endpoint caps alone cannot complete the golden-scale cusp seed

For the intrinsic **affine A11 × affine A2 × affine E6** seed with Euclidean factor inradii `(1,φ,1)`, where `φ=(1+√5)/2`, every acute Coxeter cap supplying a finite or ideal opposite endpoint of an original cusp edge misses two explicit points. This includes arbitrary real cap normals, every positive radius, all finite endpoint labels, and the isolated dihedral families and parallel limits.

Consequently no finite-volume polytope can be formed by retaining this seed and adding **only such endpoint caps**. This is a restriction on the construction method, not on the intrinsic seed: an interior sphere cap centered at either missed point can lie outside the endpoint families. In particular the new radius-one caps at half-vertices are deliberately outside this theorem.

The [exact certificate](fresh_h20_intrinsic_caps_20261014_endpoint_bounds.json) uses the frozen [complete endpoint census](fresh_h20_intrinsic_caps_20261014_endpoint.md). It checks 16,272 transported finite-pool profiles, including semidefinite cases, and handles every large-label tail continuously. No numerical truth value or finite label cutoff enters the proof.

## The two points and a strict spacelike witness

Let C be the product of the three Euclidean simplexes, with incenter zero. Its outward unit normals are n_i and offsets are b_i, equal to1 in the A11/E6 factors and φ in the triangle. Let z₀ be the vertex omitting walls `[4,23,22]`. The points are

`x₀=0`, and `x₁=z₀/2`.

Both are strictly inside C. For x₁ the gaps `b_i−n_i·x₁` are half the offsets at the nineteen retained walls. At the three omitted walls they are respectively `13/2`, `2φ`, and `5/2`. Its square is `x₁²=76+√5/2`.

For every permitted endpoint sphere with center w and radius R, the exact bounds give

`|x₀−w|²−R² > 100`, and `|x₁−w|²−R² > 1/4`.

Thus the result is stronger than a feasible ideal point. In the standard Lorentz model with form `|x|²+z²−t²`, assign the formal upper-half-space height square `y²=−1/8` and put

`X=(x₁, (x₁²−1/8−1)/2, (x₁²−1/8+1)/2)`.

Then `X²=1/8>0`. Each vertical-wall pairing is the negative of its positive gap. A unit sphere-wall normal has pairing

`−(|x₁−w|²−R²−1/8)/(2R)<0`.

Hence X is strict at every allowed wall. Its time coordinate is positive, and its pairing with the upper-half-space axis point at height H is `−(H²+x₁²−1/8)/2<0`. For any finite collection of sphere caps, take H larger than all radii to obtain a common strict timelike point. The strict spacelike escape proves infinite volume for the resulting polytope.

## Bounding every radius

For one fixed endpoint profile at a product vertex z, solve `n_i·y=cos(π/m_i)` on its nineteen incident walls. The sphere center is `w=z−Ry`, and the complete endpoint census gives `p=y²≤1`.

No additional vertical facet is lost in this parameterization. The three positive affine relations force an acute added wall to have nonpositive cusp pairing. If that pairing were zero, every incident-wall pairing would vanish; their full horospherical rank would force the added normal to lie on the cusp null line, contradicting its positive square. Every new acute spacelike facet is therefore a sphere cap, oriented to retain the cusp above it.

For each of the three omitted walls, the signed distance-to-radius ratio has the exact form

`A_f/R−B_f`, with `A_f>0` and `B_f≥0`.

Every admissible ratio is zero or at least `1/2`: a nonright finite angle has cosine at least `1/2`, while parallel and divergent walls have ratio at least one. Define

`R₀=min_f A_f/(B_f+1/2)`.

If R is larger than R₀, at least one omitted ratio must be exactly zero. Therefore all possible radii are contained in `(0,R₀]` together with the at most three exceptional values `A_f/B_f` for positive B_f. The checker retains every exceptional value whose other two ratios are zero or at least `1/2`, then takes the largest remaining bound. Equality cases and right contacts are included. The relaxation intentionally retains many ratios that are not Coxeter angles; this can only make the exclusion harder.

For the central-E6 deletion profiles with zero B values in both the triangle and E6 factors, a further all-label gap improves the bound. The E6 ratio is `4/R>0`, so `R≤8`. The triangle ratio is then at least `3φ/8>1/2`; the next finite Coxeter cosine is `1/√2`. It follows that `R≤3φ√2`. This removes the only five relaxed profile occurrences that could otherwise cover x₁. The checker applies the improvement wherever its hypotheses hold, not only to those five occurrences.

The balls are nested as R increases: their centers move at speed `|y|≤1`, while their radius increases at speed one. It suffices to test the largest relaxed radius. With d the nineteen gaps of a test point at the selected vertex and B the square-two edge Gram, write

`F=2dᵀB⁻¹d`, `t=−2dᵀB⁻¹cos`, and `p=2cosᵀB⁻¹cos`.

The ball power is `F+2Rt+(p−1)R²`. At a bound `R=A/C`, the checker evaluates the exact numerator

`F C²+2ACt+(p−1)A²`.

All 16,272 literal finite-pool cases exceed `C²/4`. The power is a concave polynomial in R and its value at R=0 also exceeds `1/4`, so the same uniform lower bound holds throughout the relaxed interval. The incenter is fixed by the seed symmetries; its 260 representative tests exceed `100 C²`. Arithmetic is in the exact multiquadratic basis of the endpoint census, using rational root enclosures solely to prove strict signs. The smallest approximate diagnostic margin is `0.3277907`; the theorem uses only the exact rational lower bound `1/4`.

## The infinite dihedral family

Only the middle-E6 deletion has an isolated A1. Put `c=cos(π/m)` for its sole contact, and include the parallel limit c=1. It is safe to enlarge the parameter range to `1/2≤c≤1`.

Here the E6 omitted ratio is `6/R−c/2`, while the untouched triangle gives `R≤6φ`. The exceptional E6-right radius `12/c≥12` is impossible because `6φ<12`. Therefore

`R≤12/(1+c)`.

For x₁, the checker verifies the distance square F to each of the 108 relevant vertices exceeds100, and the isolated-wall gap is exactly `1/2`. At the limiting radius, `(1+c)²` times the ball power exceeds

`232c²+188c−44`.

Writing `c=(1+t)/2`, with `0≤t≤1`, turns this into `108+210t+58t²`. Division by `(1+c)²≤4` gives power greater than27. For the incenter, F exceeds300 and the isolated-wall gap is1. After subtracting `100(1+c)²`, the remaining lower polynomial is `320c²+376c+56>0`. Thus both point bounds cover all labels in the tail and its parallel limit without sampling them.

## Consequence for the constructive search

Closing all 252 original edges is insufficient for finite volume. Even allowing every endpoint contact profile and radius compatible with the retained seed leaves these points exposed. Caps through other faces or through the interior are necessary if the seed is kept.

A radius-one cap centered at x₁ cuts the saved escape immediately. It cannot itself close an original cusp edge: its center is farther than one from every original product vertex. The separate construction checks its actual seed angles and joint compatibility; this endpoint theorem makes no claim about that expanded family.

```sh
python scripts/fresh_h20_intrinsic_caps_20261014_endpoint_bounds.py --verify
```

The [producer](../scripts/fresh_h20_intrinsic_caps_20261014_endpoint_bounds.py) binds the frozen endpoint source and census. Generation and saved replay passed. The conclusion is restricted to the golden scale, the retained seed, and caps incident to an original cusp edge; it is not an obstruction to arbitrary intrinsic caps or wall replacement.

The [independent review](fresh_h20_intrinsic_caps_20261014_endpoint_bounds_review.md) also passed. Its separate SymPy arithmetic rechecks all 260 representative radius bounds, 90 right-angle exceptions, eight tightenings, five critical margins, and the symbolic dihedral-tail bounds. The full transported 16,272-case census is covered by source/proof review rather than a duplicate implementation.
