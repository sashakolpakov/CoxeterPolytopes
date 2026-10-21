# A local obstruction to a single-facet dimension lift

Replacing one facet by two sloping facets is a natural way to try to raise the dimension. The exact local Coxeter constraints exclude this particular construction for every facet of the stored Bugaenko H8, Vinberg H19, and Borcherds H21 chambers. This is an obstruction to the specified lift, not to arbitrary higher-dimensional polytopes.

## The unrestricted-label local argument

Normalize the selected old wall normal `r` to square two. Retain all other old normals at added coordinate zero, and replace `r` by

```text
r+ = s+ r + sqrt(2(1−s+²)) e,
r− = s− r − sqrt(2(1−s−²)) e,
```

where `e²=1`, `e` is orthogonal to the old space, and `0<s+,s−<1`. This covers unequal slopes as well as equal slopes, after positive rescaling of the new normals.

Suppose `u` is an old finite, nonright neighbor of `r`. The span of `r,u,e` is positive definite. Thus `r+,r−,u` has a positive definite Gram matrix, and its two pairs with `u` are nonright. A positive definite Coxeter triple cannot have all three pairs nonright: with normalized negative off-diagonal magnitudes `a,b,c≥1/2`, its determinant `1−a²−b²−c²−2abc` is nonpositive. Therefore the two replacement walls must be orthogonal. This gives `s+²+s−²=1`.

If the original angle is `π/m` and the replacement angles are `π/a,π/b`, then

`cos²(π/m) = cos²(π/a) + cos²(π/b) < 1`.

With integers `a,b≥3`, the only unordered possibilities are `(3,3)`, `(3,4)`, and `(3,5)`: if both are at least four the sum is at least one, and if one is three the other must be below six. These give original labels **4, 6, and 10**, respectively. This argument has no upper-label search cutoff. In particular, an old label three or five immediately forbids the lift.

For original label six, the two squared slopes are `1/3` and `2/3`. If `r` also has a retained parallel neighbor, the first new pair has squared normalized inner product `1/3`. This lies strictly between `cos²(π/3)=1/4` and `cos²(π/4)=1/2` and is not a Coxeter angle. Thus this case is also excluded.

## Exact source checks

The [certificate](fresh_dimension_wedge_20261013.json) checks the stored exact Gram matrices and supplies an obstructing neighbor for every facet:

| Source chamber | Facets excluded | Local reason |
|---|---:|---|
| Bugaenko H8 | 11 | Every facet has a finite label three or five |
| Vinberg H19 | 24 | A finite label three |
| Vinberg H19 | 12 | A finite label six together with a parallel neighbor |
| Borcherds H21 | 210 | Every facet has a finite label three |

```sh
python scripts/fresh_dimension_wedge_20261013.py --verify
```

The [checker](../scripts/fresh_dimension_wedge_20261013.py) uses exact rational or quadratic-field Gram entries. It does not assume that a floating-point angle close to an allowed label is equal to that angle.

An [independent source/argument review and stdlib replay](fresh_dimension_wedge_20261013_independent.md) also pass. The separate checker reconstructs the three source Gram matrices and all 257 neighbor obstructions with integer, Fraction, and separate quadratic-field operations.

## What must change in a future lift

The obstruction permits multiple old facets to move together, replacement normals whose old projections are not proportional to a single facet, or changing the starting chamber. The new coupled-cap and cusp-direction constructions use such additional freedom.

Simply deleting an obstructing old facet has its own limitation. If every new old-space projection lies in the positive cone of the retained old normals, then an omitted old facet normal gives a spacelike feasible point `(v,0)`, by the nonpositive old pairings. In any proposed full-dimensional chamber, let `H` be an actual strict interior point of the enlarged chamber. Then `(v,0)+εH` is strictly feasible and still spacelike for sufficiently small positive `ε`, including when some new wall has zero old-space projection. Such a deletion cannot produce finite volume without new projection directions.

For a noncompact source there is a separate necessary condition on this simple wedge: every old cusp must lie on the selected facet. Otherwise a cusp `δ` with `(δ,r)<0` admits a small nonzero added coordinate while satisfying the new walls; a small perturbation toward the center makes all other inequalities strict and leaves positive norm. This produces an infinite-volume escape. Neither this observation nor the local angle obstruction is a theorem against general dimension extension.
