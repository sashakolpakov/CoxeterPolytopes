# Independent D4/P2 geometry and literal-host audit

**PASS.** The separate [checker](../scripts/fresh_h9_d4_prism_coupling_20261014_review.py) and [certificate](fresh_h9_d4_prism_coupling_20261014_review.json) reconstruct the frozen [producer geometry](fresh_h9_d4_prism_coupling_20261014.json) using independent `Fraction` arithmetic in `Q(sqrt(5))`. They import no producer or other repository arithmetic. They also use a different complete ray-enumeration algorithm.

| Cone | Finite extreme rays | Ideal extreme rays | Spacelike extreme rays |
|---|---:|---:|---:|
| Compact H5 P2 source | 10 | 0 | 0 |
| Coupled 11-wall H9 seed | 14 | 0 | 16 |
| Coherent 15-wall extension | 49 | 0 | 89 |

Every exact projective ray, squared norm and complete incident-wall set agrees with the primary census. This is a complete **linear-cone extreme-ray census**, not a list of every point of the ideal boundary. The source is compact; both H9 chambers remain infinite-volume.

## Independent construction and complete clipping

The audit rebuilds the six-dimensional P2 form, its seven normals, the D4 Gram matrix, the outer weight of square one, and the coupled ten-dimensional form directly from their formulas. Exact LDL elimination gives signatures `(5,1)` and `(9,1)`. All eleven initial normals have square two, and the saved center is strictly interior and timelike.

For each cone, the checker first selects a full-rank set of defining inequalities. Their inverse gives every ray of the initial simplicial cone. It then adds each remaining halfspace. Old feasible rays survive. Every positive-negative ray pair is tested: it creates an extreme new ray exactly when the rank of their common active inequalities is `dimension−2`. This is the complete double-description update, including nonsimple intermediate cones. A common active face of dimension greater than two cannot create an extreme ray after only one further cut. The initial full-rank cone is pointed, and the strict center ensures that every intermediate cone is full-dimensional.

The extension's ray counts progress **10 → 30 → 50 → 70 → 106 → 138**. This avoids the primary's enumeration of all 5,005 nine-wall supports. At completion the checker again verifies each ray's feasibility and active rank nine. The sum of normals is timelike and has strictly negative pairing with every ray, giving a bounded projective slice and a consistent time orientation.

## Facets, polar acceptance and escape

All sixteen seed spacelike extreme rays are independently tested as polar normals. Exactly four pass, at saved polar indices `0,1,4,10`, and their six mutual relations are strictly dotted. The finite-angle table is complete here: a squared Coxeter cosine in `Q(sqrt(5))` has label among `2,3,4,5,6,10`; the other real cyclotomic degree-two possibilities require `sqrt(2)` or `sqrt(3)`.

All 105 extension pairs are checked: **81 right, nine label 3, one label 4, three label 5 and eleven dotted**. For each of the fifteen normals `n`, the exact projection

`h − (h,n)n/(n,n)`

is future timelike, lies on that wall, and satisfies every other inequality strictly. Thus all fifteen are actual facets, including the four additions.

The saved rational escape is independently decoded without floating-point arithmetic. Its square is positive, its center pairing is negative, and all fifteen wall pairings are strictly negative. The checker reproduces every archived pairing. A segment from the strict timelike center to this strict spacelike vector crosses the future light cone while retaining strict inequalities; hence there is an open ideal-boundary patch and the chamber has infinite volume. This also certifies infinite volume of the larger eleven-wall chamber.

## Literal arithmetic host: checked scope

The audit separately reconstructs the basis replacing one D4 root by the outer weight. Its Gram matrix is integral over `Z[phi]`, has determinant `(1−sqrt(5))/2`, a unit, and contains the saved norm-one basis vector. Its real signatures are `(9,1)` and positive definite at the conjugate embedding.

All fifteen reflection matrices are rebuilt in that basis and checked entry by entry against the certificate. Each matrix is integral over `Z[phi]`, preserves the Gram matrix and squares to the identity; therefore its inverse is integral too. This establishes the literal odd unimodular host and its reflection action.

The separate use of Bugaenko's classification/uniqueness and nonreflectivity theorem belongs to the main report. This audit does **not** independently reprove that external classification or infer a host isometry merely from displayed invariants. Its geometric infinite-volume conclusion needs no external classification, and it makes no claim excluding arbitrary nonintegral added caps.

Generation and a fresh saved replay pass:

```sh
python scripts/fresh_h9_d4_prism_coupling_20261014_review.py --verify
```

The review binds the primary source, primary JSON, the source's arithmetic dependency and its own checker. It does not bind the producer's prose report, so adding this review link there is safe.
