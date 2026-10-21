# Conditional completion of the two E8-branch cusp links

The literal21-root diagram of the proposed E8 triangular branch contains two feasible affine F4 null rays. Each has incident link `affineF4+E8+A6`, of Euclidean rank18. **Three explicit new roots complete each link to `affineF4+affineE8+affineA7`, of rank19.** Both completions can be made simultaneously by sufficiently large positive heights.

This is an exact conditional existence certificate. It assumes exact promotion of the [numerical162-wall recipe](fresh_h20_e8_triangle_cluster_20261016.md), with the stated old dotted gates and a common strict future center. It does not itself promote that recipe, supply a certified168-wall realization, or establish finite volume. The [source](../scripts/fresh_h20_e8_cusp_completion_20261017.py) and [data](fresh_h20_e8_cusp_completion_20261017.json) pass generation and fresh saved replay and are final.

## Literal cusp incidence

The core is `T+[158,159,160,161]`, where T is the17-root support in the input recipe. It has exact determinant `-2(5+sqrt(2))` and inertia `(20,1,0)`. The E8 component is

`[158,12,13,14,15,61,87,151]`.

The two affine paths, with the displayed positive kernel normalization, are

| Cusp | Affine F4 path | Positive kernel | Remaining A6 |
|---|---|---|---|
|0|`[2,1,159,161,160]`|`(sqrt(2)/2,sqrt(2),3sqrt(2)/2,2,1)`|`[4,5,6,7,8,9]`|
|1|`[1,159,161,160,10]`|`(sqrt(2),2sqrt(2),3,2,1)`|`[3,4,5,6,7,8]`|

The resulting null vectors have mutual product−1. Their positivity and a strict chamber center make them nonzero future null vectors. Acuteness gives their nonpositive products against all retained walls. The source checks every core product exactly. Every old wall outside T has strictly negative product with each null vector because each path contains new partners that are strictly dotted to every such old wall. Thus the listed19 roots are the full incident set. Their Euclidean rank is `4+8+6=18`; the horospherical link is deficient by one dimension.

## Three new roots at one cusp

Write nu for either null vector. Let theta be the E8 highest root, with coefficients

`(6,5,4,3,2,4,2,3)`

in the E8 order above. It has square2 and only one nonzero simple-root product, equal to1 at root15. In the A6 component let omega_f and omega_l be the endpoint fundamental weights. Exact inversion gives

`omega_f²=omega_l²=6/7`, `omega_f.omega_l=1/7`, and `omega_f+omega_l=sum(A6 simple roots)`.

Choose a positive vector z orthogonal to all incident roots and to nu, with square8/7. Such a direction exists because the quotient `nu-perp/R nu` is positive definite of dimension19 and the old incident span has rank18. The script constructs z explicitly inside the literal21-dimensional Lorentz core, using the auxiliary timelike vector `nu0+nu1` and exact orthogonal projections. It does not assume that auxiliary vector is a strict chamber center.

Set

`bE=-theta`, `bA=-omega_f+z`, `bB=-omega_l-z`.

Their Gram is

```
2  0  0
0  2 -1
0 -1  2
```

The last two have only the required endpoint products−1 against A6. Adding arbitrary heights preserves all these products:

`rE=bE+hE nu`, `rA=bA+hA nu`, `rB=bB+hB nu`.

For strictly positive heights the affine relations are

`theta+rE=hE nu`,

`sum(A6)+rA+rB=(hA+hB)nu`.

Together with the original positive F4 relation, these give three affine components of ranks4,8,7 with the same future null ray. Their normals span `nu-perp`, of dimension20, and their horospherical chamber is the product of bounded Euclidean alcoves. There are22 incident roots and a full rank19 cusp link.

## Finite heights and simultaneous compatibility

For an old nonincident wall N, `nu.N<0`. A sufficient explicit bound is

`h > (b.N+sqrt(2 N²)+1)/(-nu.N)`.

It makes the new norm2 wall strictly dotted to N. If H is a strict future center for the original chamber, the further bound

`h > (b.H+1)/(-nu.H)`

preserves that center strictly. There are finitely many old walls, so a finite positive integer satisfying all these inequalities exists.

The two null vectors have negative mutual product. With a common height L, every pair consisting of one new wall from each cusp has product

`-a L²+b L+c`, with `a>0`.

The source saves all nine exact coefficient triples. The sufficient bound

`L >= 1+(abs(b)+abs(c+3))/a`

makes this product strictly less than−3. Taking the maximum of finitely many old-wall, center and cross-cusp bounds therefore gives simultaneous coherence and Coxeter compatibility. Every wall is actual: projecting the common strict timelike center onto each spacelike wall yields a timelike relative-interior point, since all other wall products are nonpositive.

## Why this does not close the unfinished H3 boundary

At cusp0 the T roots3 and10 have strictly negative null products; at cusp1 these roots are2 and9. The sufficient large-height construction makes every new wall dotted to those T roots. Consequently none of the six new walls meets the actual T-face. Its H3 boundary is unchanged. This certificate resolves the local cusp-link deficiency conditionally; the finite boundary still requires different endpoint contacts, and no global finite-volume conclusion follows.
