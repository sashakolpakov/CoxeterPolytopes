# A contact-producing repair of the new E7 chamber cusp

**Exact positive construction, still infinite volume.** Five norm-two caps extend the [90-wall chamber](fresh_h20_e7_finite_edgewalk_20261014.md) to 95 actual facets. They complete its designated deficient cusp to a bounded rank-19 Euclidean link, preserve all 84 verified finite vertices and both inherited bounded cusp links, and close 12 more outgoing edges at distinct finite endpoints. Six of these new endpoints also connect by actual edges to the repaired proper ideal vertex. The previous strict escape survives.

The [primary checker](../scripts/fresh_h20_e7_edge_cusp_repair_20261014_exact.py) and [exact certificate](fresh_h20_e7_edge_cusp_repair_20261014_exact.json) contain the complete pair, link and endpoint checks. The [independent review](fresh_h20_e7_edge_cusp_repair_20261014_review.md) reconstructs all five caps from the preceding geometry using separate sparse squarefree arithmetic. The two constructions' serialized normals agree exactly.

## Five exact quotient roots

Keep the 90-wall numbering and put `ν=N9+N88`. The previous frontier certificate supplies a vector `L` perpendicular to all 20 cusp-incident normals, with `L²=98/27` and `L·ν=0`. Define

```text
N90 = −N12−N13−N18−N31−N65 + ν,
N91 = −N15 + 2ν,
N92 = −√2(N17+N79)−N83 + ν,
N93 =  (3√3/7)L + (1+4√3/7)ν,
N94 = −(3√3/7)L + (10√3/7)ν.
```

The first three transverse parts are the negative highest root of A5, the negative A1 root, and the negative highest **short** root of the scaled B3 component. For the latter, the scaled old roots `√2N17, √2N79, N83` have squares `4,4,2`; their sum has square two and the required nonnegative products with the simple roots. The final two transverse parts have square two. Since every transverse part is perpendicular to ν, all five norms are two and their mutual products are independent of the chosen shifts.

These shifts produce actual finite contacts outside the old cusp link, rather than merely distant dotted walls. In particular N93 meets old wall 8 at π/3, and N94 meets old walls 0 and 10 at π/6. The checker supplies eight explicit representative ridge points, future timelike and strict on every other wall.

All 460 new unordered pairs pass: 107 right, ten π/3, one π/4, two π/6, five parallel and 335 dotted. The common strict timelike center remains. The acute center-projection argument proves all 95 facets actual; the five new facet points are also explicitly checked. All 84 previously verified finite points are strict at every new cap.

## Complete bounded link at the designated cusp

The full incident set at ν has 25 walls. Its connected components are affine D8, affine A5, affine C3, and three affine A1 factors, with Euclidean ranks

```text
8 + 5 + 3 + 1 + 1 + 1 = 19.
```

The proof does not infer boundedness from rank alone. For every component, the checker verifies a positive-definite principal minor, a strictly positive Gram kernel, and a weighted physical-vector sum equal to a **positive** multiple of ν. The component normal spaces span all 19 horospherical directions. A recession vector must therefore be orthogonal to every component and hence vanish. The common strict center gives a nonempty link interior, so the link is bounded and ν is a proper ideal vertex. Both previously completed cusp links remain unchanged because the five new caps are strict at their null rays.

## Twelve finite endpoints and six finite-to-ideal edges

The bounded continuation test examines the 541 previously unresolved germs at the 40 newest vertices of the 90-wall chamber. On each ray `V−tD`, where D is the dual of the released wall, it computes the first positive intersection with the five new caps. It accepts a new endpoint only if it is future timelike and has precisely the required 20 zero wall products, with all other 75 inequalities strict.

There are 12 distinct accepted endpoints: six on N92 and six on N93. Their norms and complete wall signs prove actual finite vertices; the nonzero released-dual product gives independence of the 20 incident normals. Strictness along the intervening segment proves twelve actual finite edges from the previous vertices. All previously certified 145 finite edges remain, since both old endpoints are strict at the five new caps.

Each of the six N93 endpoints shares a positive-definite 19-wall support with ν. Both endpoints are feasible, and every other wall is strict at at least one end, so the open segment is strict outside that support. These are six actual edges from finite vertices to the repaired proper ideal vertex. Thus six old corners now lie on certified finite–finite–ideal paths.

This is a bounded continuation census, **not a complete finite-vertex census for the 95-wall chamber**. It proves at least 96 finite vertices. It does not classify all remaining germs or all new cusp directions.

## Global status and verification

The inherited strict escape remains feasible at all 95 walls, with positive square

```text
259854917001622621 / 62500000000000000.
```

Consequently the chamber still has infinite volume. The result repairs one full cusp and closes specified edges; it does not establish global coverage.

The [numerical event scout](../scripts/fresh_h20_e7_edge_cusp_repair_20261014_scout.py) considered fixed quotient roots and a bounded set of nonincident contact labels. It proposed the shifts above; every acceptance in the primary certificate is exact and independent of the scout. There is no all-label exclusion claim. Exact generation and saved replay passed; the independent review is recorded in the linked certificate. Frozen predecessor files were preserved.
