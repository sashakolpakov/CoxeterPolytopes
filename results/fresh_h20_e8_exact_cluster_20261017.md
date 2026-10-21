# Exact promotion of the162-wall E8 triangular branch

The four-wall [E8 recipe](fresh_h20_e8_triangle_cluster_20261016.md) now defines an **exact coherent Coxeter chamber in H20 with162 actual facets**. There are exactly four new finite vertices, bringing the inherited lower bound to188. A new E8-supported two-face is a compact triangle. Five previously proper cusps remain, and two additional ideal rays have deficient links `affineF4+E8+A6`.

The chamber has infinite volume. A strict spacelike escape is explicitly certified; the new finite face does not close the remaining open region. No complete ambient vertex census or new finite-volume example is claimed.

The [physical geometry export](../scripts/fresh_h20_e8_exact_cluster_20261017_geometry.py) replays the frozen [independent158-wall construction](fresh_h20_compact_quad_cap_20261016_review.md). Its [data](fresh_h20_e8_exact_cluster_20261017_geometry.json) are outward rational enclosures of the inherited exact algebraic vectors. The [new source](../scripts/fresh_h20_e8_exact_cluster_20261017.py) and [certificate](fresh_h20_e8_exact_cluster_20261017.json) then construct the four new roots without importing the numerical scout's arithmetic. Both pairs pass generation and fresh saved replay. This is a primary exact certificate; a separate review of these four new constructions is not part of this report.

## Exact definition and branch selection

Use the rank17 spherical support

`T=[1,2,3,4,5,6,7,8,9,10,12,13,14,15,61,87,151]`,

with type `A10+A4+A2+A1`. Each new normal has square2. Its prescribed full products are:

| Root | Nonright T products | Additional full products |
|---|---|---|
|158|−1 at12,61,151|95: `-75089/10000`;142: `-116551/10000`;145: `-36367/5000`|
|159|−1 at1|158:0;95: `-26659/2000`;142: `-209307/10000`|
|160|−1 at10|158:0;159:0;95: `-142043/10000`|
|161|none|158:0;159:`-sqrt(2)`;160:−1|

Unlisted T products are zero. First take the exact T projection using its rational inverse Gram. Three further linear equations determine an affine line in the four-dimensional Lorentz complement. Its square equation gives two algebraic branches. A projected old root31 supplies a nonzero kernel direction in every step. The residual computation uses the basis `[95,142,145,31]`, which avoids repeated interval cancellation in the larger physical coordinates.

These are exact definitions by inversion and positive square roots. Interval endpoints do not replace the vectors. Prescribed norm/contact equalities follow algebraically from projection and the quadratic equation; outward rational intervals prove nonzero pivots, positive radicals, and all strict inequalities. Previously imposed exact products are reused as identities in later frame Grams.

Some quadratic branches both satisfy the old dotted gates. The intended archived branch is selected by the following rational open intervals for the full product with145:

`159.145 in(-11.989,-11.988)`,

`160.145 in(-10.659,-10.658)`,

`161.145 in(-15.067,-15.066)`.

The intervals are selection inequalities, not extra equality constraints. Every chosen product enclosure lies strictly inside its interval. The source saves both branch tests and the complete algebraic recipe.

## All gates, coherence and finite vertices

All638 new/prior products are certified:564 are strictly dotted,67 are right, six have label3 and one has label4. Every new wall is strictly dotted to every old wall outside T. The inherited future timelike center remains strictly inside all162 halfspaces. Orthogonal projection of that center onto each wall gives a timelike point strict against every other wall, proving that every wall is an actual facet.

The full21-root core `T+[158,159,160,161]` splits into E8 and a13-cycle with one label4 bond. Its exact determinant is `-2(5+sqrt(2))`. Since a finite corner cannot contain a dotted pair, every new finite corner is supported inside this core. Deleting an E8 node leaves the indefinite cycle. Deleting one cycle node gives a path, and a standalone rational leaf-elimination census checks all13 possibilities. Exactly omission of1,159,160 or161 is positive definite. All four resulting points are constructed by projection, and their strict products against every other wall are certified. Their types are `E8+A12` or `E8+B12`, as listed individually in the input recipe.

More explicitly, omitting1 or160 gives `E8+B12`; omitting159 or161 gives `E8+A12`. The three corners that retain all of T form the new compact triangle on `T+[158]`. Its side normals159,160,161, after projection to that two-face, have Gram

```
 12/11   -1/11   -sqrt(2)
 -1/11   12/11   -1
 -sqrt(2) -1       2
```

The three two-by-two determinants are13/11,2/11,13/11, while the full determinant is negative. The three complete finite corners and strict old gates therefore certify the whole compact triangular face. Its intrinsic angles need not be Coxeter angles.

All old finite vertices persist, not just a sampled set. At any old finite vertex v, expand a new normal in the20 old spherical incident normals plus v. Its spherical coefficients are nonpositive because the spherical inverse is entrywise nonnegative. If the v coefficient were nonpositive, the new normal would pair nonnegatively with the common strict center, a contradiction. Hence its product with v is strictly negative. Applying this to all four new walls preserves the inherited lower bound184 and adds four distinct new vertices, giving188. The five previously proper cusp rays also have strictly negative products with every new root.

## Exact remaining escape and applicable cusp completion

The two positive affine F4 rays from the [conditional cusp certificate](fresh_h20_e8_cusp_completion_20261017.md) are now checked against all162 exact normals. Each has exactly19 incident roots and Euclidean link rank18. Thus the ambient hypotheses of that certificate are discharged for this branch. It supplies a simultaneous finite-height mechanism to complete both links, but this package does not list those six additional walls or call their existence an H3 boundary closure.

For the first null ray nu0, use the certificate's transverse vector z0 of square8/7 and the common strict center H. The exact vector

`E=nu0+(1/4096)z0+(1/4096)^3 H`

has certified square approximately `3.31738279589e-8`, negative product with H, and strictly negative product with all162 walls. The segment from H toward E reaches the future null boundary while every wall inequality is still strict, so the chamber has an open ideal boundary patch and infinite volume.

The earlier numerical H3 ray count and exterior-radius scores remain discovery data in their own frozen package. This certificate promotes the four exact normals, all actual facets, the complete new finite-corner census and the escape; it does not promote that numerical face enumeration to an exact cone census.
