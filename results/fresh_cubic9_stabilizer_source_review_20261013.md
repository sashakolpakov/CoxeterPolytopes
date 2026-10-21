# Independent review of the proposed conductor-nine H3 stabilizer argument

The full saved certificate passes independent rerun. I reviewed the mathematical reduction and enumerator source, then ran [the final Fraction verifier](../scripts/fresh_cubic9_global_nonreflectivity_verify_20261013.py) with `--verify`. It repeated all 142 time/norm controls, checked all 136 dominant positive candidates at both endpoints, and verified all 38 oriented zero roots. The [saved certificate](fresh_cubic9_global_nonreflectivity_20261013.json) therefore proves nonreflectivity of the full H3 host, not merely infinite volume of a bounded prefix. The argument below explains the scope and the step from the finite census to the unbounded conclusion.

Let the form be `−αx² + Q(y)`, with `α>0` and `Q` positive definite at the distinguished real embedding. Orient a positive-height root `r=(x,y)` by `x>0`, and write its norm as `q>0` and its height as `h=αx²/q`. If a future timelike point `u=(u₀,z)` satisfies `(u,r)≥0`, Cauchy–Schwarz for Q gives

`α²x²u₀² ≤ (q+αx²) Q(z)`.

Because `Q(z)=αu₀²+(u,u)` and `(u,u)<0`, this implies

`h ≤ αu₀²/(-(u,u)) − 1`.

Thus only roots below that exact bound can prevent strict chamber membership. This applies separately to both p and gp.

For a point strictly inside the finite zero-height chamber, folding a root toward the dominant chamber increases its pairing with the point. Indeed, if a zero simple root a has `(r,a)>0`, then replacing r by `r−2(r,a)a/(a,a)` increases `(u,r)` because `(u,a)<0`. Finite Weyl reduction terminates in the dominant chamber, preserves the time coordinate and the norm, and is an integral lattice isometry. The same reflection increases pairings with both points when both are in that zero chamber. Consequently, any obstructing mirror yields an obstructing dominant candidate at the same height. A complete dominant census with strictly negative pairings therefore certifies strict membership in the full reflection chamber, including mirrors not accepted as facets by the prefix search.

An integral time-orientation-preserving isometry g permutes the full reflection arrangement. If p and gp are both strictly in the same chamber C, then `gC=C`. In dimension four, exact `trace(g)>4` proves g has infinite order, since every eigenvalue of a finite-order real matrix has absolute value one. The resulting infinite chamber stabilizer forces infinite index of the reflection subgroup in the arithmetic integral orthogonal group, so the full H3 host is nonreflective. This is a statement about the full norm alphabet; a norm-two-only stabilizer would not suffice.

I read [the trace-ellipsoid auditor](../scripts/fresh_cubic9_ellipsoid_audit_20261013.py) and its separate Cartesian enumerator. For a totally positive cap c, `Tr(x²/c)≤3` is necessary for the requested three embedding inequalities. The rational matrix `Q_ij=Tr(e_i e_j/c)` is positive definite; the code asserts its exact LDL reconstruction. Reverse-coordinate recursion enumerates all integral coefficient vectors in the ellipsoid. With offset `a/b`, its integer interval is exactly `ceil((-R−a)/b)` through `floor((R−a)/b)`, where R is the floor of the appropriate rational square root. The code's lower endpoint `−floor((R+a)/b)` is the same expression. Final exact embedding tests remove the extra ellipsoid points. I found no omitted-bound issue. Direct comparisons with the separate unpruned box alphabet for caps `1,2,β,2β,4` agreed (sizes `3,3,1,1,11`). Time coordinates still use the independently certified coefficient boxes.

The [forced-face argument](fresh_cubic9_high_label_faces_20261013.md) is a separate necessary-condition step. Removing three zero-height `I₂(18)` components from a hypothetical compact chamber gives a compact Coxeter face in the orthogonal `B+<1>` lower host. Every corresponding integral reflection extends by the identity to the direct-sum host. A full lower-host nonreflectivity theorem therefore obstructs compactness in the specific H9 `B⁴+<1>` host with the same α. The stronger finite-covolume exclusion in the final report is also valid: positive definiteness at another real place makes the form anisotropic over K, so its ambient arithmetic quotient is cocompact. Any finite-covolume integral reflection subgroup would have finite index in that ambient group and hence would be cocompact, contradicting the compact obstruction. This does not classify arbitrary H9 polytopes or exclude unrelated lattices.

The final replay also checks the proposed matrix and its inverse are integral, verifies its exact Gram congruence and trace `1+3T+T²>4`, verifies both endpoints are future timelike with height bounds strictly below 3000, and rejects a corrupted matrix and a nonregular endpoint. The replay completed in approximately 30 seconds.

```sh
PYTHONPATH=scripts python scripts/fresh_cubic9_global_nonreflectivity_verify_20261013.py --verify
```
