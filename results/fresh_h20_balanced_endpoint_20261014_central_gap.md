# The central cap remains isolated from every original-edge endpoint cap

**At both scales `(1,1,1)` and `(1,1,phi)`, the radius-two central cap together with any collection of original-edge endpoint caps still has a uniform strict escape.** Thus a finite-volume completion retaining this core must use a cap meeting another face or the interior of the alcove. This is a complete endpoint-family statement, not a sampled bank obstruction, and it does not exclude the intrinsic cusp seed itself.

The [complete endpoint census and radius bound](fresh_h20_balanced_endpoint_20261014.md) prove that every such cap sphere, of center `w` and radius `R`, satisfies `P=|w|²−R²≥68` and `R≤14`; at golden third-factor scale, `P>75`. The radius bound includes all exceptional right and pi/3 omitted-wall radii. The [corollary checker](../scripts/fresh_h20_balanced_endpoint_20261014_central_gap.py) binds those exact data and rechecks every final bound before applying the geometry below. Its [certificate](fresh_h20_balanced_endpoint_20261014_central_gap.json) is separate from the frozen census and radius source.

The normalized negative pairing with the radius-two central cap is `(P−4)/(4R)`. It is at least `8/7>1` at unit scale and above `71/56>1` at golden scale. Hence every endpoint cap is strictly dotted to the central cap. The argument below supplies an explicit common gap rather than relying solely on the disconnected finite-contact graph.

For a test point `x`, Cauchy gives

```text
|x−w|²−R² ≥ P+|x|²−2|x|√(P+R²).
```

For the points used here, `P≥Pmin>|x|²`, so the right side increases with `P` and decreases with `R`. It is therefore bounded below by replacing `P` with its uniform lower bound and `R` with14. This monotonicity is essential: the center norm `|w|` is not separately assumed bounded by its value at minimal `P`.

At unit scales take `x=3z0/20`, where `z0` is the vertex omitting walls `(0,7,14)`. Since `|z0|²=196`, its norm is `|x|=21/10`. The central ball power is `41/100`. Using `√264<65/4`, every endpoint ball has power strictly above

```text
68+441/100−(21/5)(65/4) = 104/25.
```

At golden scale take `x=z0/8`. Its squared norm is `(49+21phi)/16`, strictly between5 and `(23/10)²`. The central power is greater than1. Using `√271<33/2`, every endpoint power is strictly above

```text
75+5−2(23/10)(33/2) = 41/10.
```

Both points lie strictly inside the product alcove because their multipliers of `z0` lie strictly between zero and one; the checker records all22 positive seed gaps. In the null coordinates `D·U=−1`, set

```text
X=(x,(|x|²−1/4)/2,1).
```

Then `X²=1/4`. Its central pairing is negative, every endpoint pairing is negative because its ball power exceeds1/4, and every seed pairing is negative. The pairing with the common timelike center `(0,128,1)` is negative. These strict future spacelike witnesses prove infinite volume for any finite compatible endpoint selection, even though the certificate allows all endpoint caps simultaneously when deriving the escape.

Interior caps, caps that do not end an original cusp edge, and wall replacement remain available. The explicit endpoint caps cutting the old half-vertex point are consistent with this theorem: that point was not the uniform escape used here.

```sh
python scripts/fresh_h20_balanced_endpoint_20261014_central_gap.py --verify
```
