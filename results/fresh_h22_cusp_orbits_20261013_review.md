# Independent source review of the H22 cusp-orbit completion

**PASS.** I reviewed the [exact verifier](../scripts/fresh_h22_cusp_orbits_20261013_verify.py) and ran

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_cusp_orbits_20261013_verify.py --verify
```

The replay passed for 15 completed sampled orbits, 1,065,708 null directions and 2,252,588 actual facets. This is a source review and fresh replay of the existing arithmetic, not a second implementation. Compatibility of the earlier 121,172 walls is an explicitly hash-anchored, previously reviewed dependency.

Closure under the saved integral isometries, together with the earlier-predecessor test, proves connectivity from each representative and completeness under those generators. It also transfers the representative's norm, future orientation, height, old-wall pair types and old cusp-link structure. Closure of the earlier timelike cap orbit is inherited from its verified certificate, so checking every earlier orbit cap against each new representative transfers to every new/old orbit pair.

The integer null-product bound is valid with the following necessary detail: all rays are primitive. Representatives have gcd one; integral isometries of the nondegenerate form have determinant `±1`, hence preserve primitivity. Distinct archived vectors therefore cannot be positive proportional future null vectors. Their Lorentzian pairing is a strictly negative integer and is at most `−1`. For height-one caps this permits only right, parallel or dotted same-sign pairs; mixed-height and height-`1/4` distinct-ray pairs are dotted. The two signs over one ray are parallel. No unenumerated finite angle is being inferred from an approximate comparison.

The bounds against the previous signed timelike caps are sufficient: `U·delta ≤ −6` at height one gives pairing at most `−3+sqrt(5/2)<−1`; `U·delta ≤ −2` at height `1/4` gives at most `−1+sqrt(5/2)/4<−1/4`. Both are strictly dotted. The corresponding center bounds and center-cap height at most four establish dotted relations uniformly for every integer label `m≥3`.

Each original cusp link has only positive affine components of total Euclidean rank 20. Its new pair sums to `2 delta`, adds an orthogonal affine `A1`, and raises the quotient rank to 21, the full H22 horosphere dimension. Every other new ray and every previous nonincident cap is strictly negative there. Positive component relations force every recession vector to vanish; the common strict center supplies a nonempty link interior. Thus all listed links are bounded. Strict-center projection and nonpositive pairings prove that every retained inequality is an actual facet.

Finally, the rational saved point has positive norm, future orientation and strict inequalities against all walls. Squared integer inequalities handle both signs of the old radical-height caps; the new rational-height inequalities and `hp+4|z|<0` handle every new cusp cap and every center label uniformly. The chamber therefore remains infinite-volume. A separate direct signed count confirms that 28 null directions cut the previous escape, each through exactly one sign: 28 cutting signed facets as well.

The conclusion is limited to the 15 sampled cusp orbits and the displayed chamber family. It is not an exhaustive old-cusp census, a new finite-volume example, or an obstruction to arbitrary further walls.
