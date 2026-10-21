# A4-to-H4 sections through generic source edges

This pilot constructs genuine finite-volume H20 sections with an exact H4 component at a new finite corner. All 108 sampled proposals show a numerical forbidden angle at a neighboring finite vertex. Nine representative failures are certified exactly, one for each sampled source edge, by an actual obtuse ridge. This is a bounded pilot, not a continuous obstruction to this mechanism.

The files are the [numerical scout](../scripts/fresh_h20_h4_edge_sections_20261013_scout.py), its [frozen proposals](fresh_h20_h4_edge_sections_20261013_scout.json), the [standard-library exact checker](../scripts/fresh_h20_h4_edge_sections_20261013.py), and its [certificate](fresh_h20_h4_edge_sections_20261013.json). They preserve the previous cusp and all-right-edge files.

## Construction

Take an actual ambient edge with twenty positive-definite incident normals, whose Gram has a norm-two A4 component. At an endpoint `a` of that chain let `w_a` be the old dual weight, so `(w_a,r_i)=delta_ai` and `w_a²=4/5`. Set

```
q = (1+sqrt(5))/4,
gamma² = q - 4/5 = (5sqrt(5)-11)/20 > 0.
```

The orthogonal complement of the twenty incident normals is a Lorentzian plane. Normalize the edge's two future endpoint rays by `-(h,v)=1`, using the source's strict timelike center h, and choose an interior point P on their segment. If b is the second endpoint, put

```
z0 = b - (b,P)P/P²,
z = +/- sqrt(gamma²/z0²) z0,
e = w_a + z.
```

Here `P²<0`, `z0²>0`, `(e,P)=0`, and `e²=q`. The section is transverse to the ambient edge. All twenty incident facets survive at the resulting finite section vertex. Their projected Gram is

`H = G - delta_a delta_a^T/q`.

Only the endpoint diagonal of A4 changes. Its adjacent squared cosine becomes

`1/[2(2-1/q)] = (3+sqrt(5))/8 = cos²(pi/5)`,

so A4 changes exactly to H4; the other finite components are unchanged. Positive definiteness is checked with exact pivots. The transversality and independent projected normals also imply that the section meets the ambient strict interior, rather than merely touching an old boundary face.

Each resulting section has finite volume inherited from the frozen finite-volume H21 source. In the Klein model, its only ideal points can be inherited old ideal vertices, and the link there is a bounded Euclidean slice of the source cusp link. The obstruction found below concerns actual angles, not volume or the validity of the initial H4 corner.

## Bounded sample

After literal support deduplication, the saved source walk contains 6,492 edge supports, including 1,660 A4 components and 822 distinct unordered endpoint-orbit pairs supporting an A4 component. These are counts in the archived walk, not a new classification under a larger automorphism group.

The pilot chooses nine distinct endpoint-orbit pairs: one ideal/ideal, four ideal/finite, and four finite/finite. It uses one A4 component on each edge, both endpoints of that chain, the three projective positions `1/3,1/2,2/3`, and both radical signs: 108 literal proposals. At each proposal it follows the outgoing section edges to neighboring finite vertices and checks all 210 source inequalities before inspecting projected angles.

The numerical outcome is 107 actual-vertex obtuse-angle proposals and one squared-cosine value strictly between `1/4` and `1/2`. All are forbidden for Coxeter ridges. These numerical results remain numerical; the exact checker certifies the following nine midpoint, positive-radical cases only.

| Source endpoint-orbit pair | A4 endpoint | Released wall | New wall | Exact actual obtuse ridge |
|---|---:|---:|---:|---|
| 9, 9 | 6 | 17 | 23 | 6, 23 |
| 9, 235 | 8 | 13 | 14 | 8, 14 |
| 9, 239 | 6 | 17 | 23 | 6, 23 |
| 9, 245 | 6 | 17 | 23 | 6, 23 |
| 9, 208 | 8 | 13 | 14 | 8, 14 |
| 15, 38 | 40 | 2 | 26 | 26, 40 |
| 16, 298 | 14 | 5 | 4 | 4, 14 |
| 16, 359 | 8 | 5 | 6 | 6, 8 |
| 16, 196 | 24 | 14 | 32 | 24, 32 |

## Actual-ridge certificate

The exact arithmetic is `K=Q(sqrt(5))` with one positive square root `sqrt(gamma²/z0²)` for each recipe. Its field norm is negative, so the radicand is not a square in K; the quadratic-field divisions and coefficient comparisons are valid.

Let I be the inverse old edge Gram. Sherman–Morrison gives the section Gram inverse, and the section-edge dual vector for releasing old wall j simplifies to

`t_j = w_j - I_ja z/gamma²`.

The checker verifies this identity and all its pairings. For the proposed first intersected wall l it sets

`lambda = (P,r_l)/(t_j,r_l) > 0`,  `X=P-lambda*t_j`.

It checks `(X,e)=0`, `X²<0`, future orientation, and **every one of the 210 source inequalities**. Exactly the nineteen retained old edge walls and l vanish; every other inequality is strict. These signs prove this is the first intersection along the section edge, independent of the numerical minimum used to propose l.

The retained nineteen projected normals form a positive-definite principal submatrix of H. The new normal's Schur complement is

`(t_j,r_l)^2/t_j² + (P,r_l)^2/P²`

`= ((t_j,r_l)^2/t_j²) * X²/P² > 0`.

Thus all twenty incident section normals are independent and positive definite. X is a genuine simple finite vertex, so every pair of its incident facets gives an actual local ridge. For the displayed pair the checker proves positive projected pairing, positive diagonal entries, and positive two-by-two determinant. The angle is therefore genuinely obtuse. No rejection is based on a pair of redundant inequalities or disjoint supporting hyperplanes.

## Replay and scope

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_h4_edge_sections_20261013_scout.py
python scripts/fresh_h20_h4_edge_sections_20261013.py --verify
```

The numerical scout is reproducible but is not used as an exact predicate. The exact verifier reconstructs each of the nine normals, initial H4 corners, and bad vertices from the raw integer source and rational position, then rechecks the full wall sign list. Frozen source and input hashes are bound in the certificate.

The final saved-result replay passed for all nine cases. A separate [source and mathematical review](fresh_h20_h4_edge_sections_20261013_review.md), with a fresh invocation of the same verifier, also passed; its [hash-bound record](fresh_h20_h4_edge_sections_20261013_review.json) explicitly distinguishes this review from an independent field-arithmetic implementation.

The construction works directly with new finite vertices on ambient open edges. The nine exact sections fail globally. This does not exclude different positions, other edge orbits, a different local spherical update, or other section normals. No new finite-volume Coxeter polytope or new dimension is established.
