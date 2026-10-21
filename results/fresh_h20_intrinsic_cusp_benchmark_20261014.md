# A different intrinsic cusp seed for H20

**Affine A6 × A6 × A7 is the most balanced next three-factor seed by worst vertex distance.** At unit inradius in every factor, all 392 product vertices have squared distance 196 from the incenter. This uniquely minimizes the largest squared vertex distance among all 229 rank-nineteen ADE products with at most three factors. It replaces the old A11 direction instead of rescaling it. This is a geometric choice for a next construction, not a finite-volume candidate.

The [standalone Fraction checker](../scripts/fresh_h20_intrinsic_cusp_benchmark_20261014.py) builds finite Cartan matrices for A1–A19, D4–D19, and E6–E8. It closes their finite reflection root sets, recovers the unique highest root and affine positive marks, and reconstructs every vertex of each unit-inradius alcove. Positive-definite finite Gram matrices and the positive affine relation prove that these are the full bounded simplexes. The [certificate](fresh_h20_intrinsic_cusp_benchmark_20261014.json) contains all 38 factor geometries and all 607 rank-nineteen products with at most four factors. No lattice integrality or arithmetic reflectivity condition enters this calculation.

| Product | Vertices | Vertex squared distance range | Minimum squared separation |
|---|---:|---:|---:|
| A11 × A2 × E6 | 252 | 302–446 | 12 |
| A6 × A6 × A7 | 392 | 196 | 84 |
| A4 × A5 × A5 × A5 | 1,080 | 125 | 40 |
| A9 × A10 | 110 | 385 | 180 |
| A19 | 20 | 1,330 | 760 |
| E7 × E6 × E6 | 392 | 45–711 | 48 |
| E8 × E7 × D4 | 360 | 65–1,667 | 36 |
| E8 × E8 × A3 | 324 | 90–2,490 | 24 |
| A1 repeated nineteen times | 524,288 | 19 | 4 |

All rows have unit inradii; the first row is therefore not the previous golden scaling. A4 × A5³ uniquely minimizes the worst distance among the 607 products with at most four factors, at the cost of a larger vertex set. The A1 product is an additional comparison, not part of that bounded enumeration.

The exceptional-factor alternatives behave differently. E6 has vertex squares 12–156, E7 has 21–399, and E8 has 40–1,240. E7 × E6² puts a few vertices close to the incenter, making finite contacts with a central cap more plausible, but leaves much farther outer vertices. Its smallest squared vertex distance, 45, is the smallest among the three-factor-or-fewer products. It is a useful second route if the aim is to couple the center to a nearby corner; it is not the best choice for uniform global covering. E8²+A3 does not improve the worst gap.

For a central hemisphere of radius two, unit-inradius vertical walls all meet it at pi/3. For a full equal-radius layer at all alcove vertices, nonpositive pairings require `R²≤d_min²/2`. The smaller value `d_min²/4` makes the nearest pair parallel; neither number by itself guarantees all other cap relations are Coxeter. These two radii are recorded separately. The earlier scale draft's confusion between them is not used here.

For A6² × A7, the necessary acute bound is `R²≤42`, with parallel-neighbor value 21. These are much larger than the unit-scale old seed's bounds 6 and 3. The prospective parallel-neighbor layer has simple compatible relations: all distinct centers have squared separation at least `84=4·21`, so pairs are parallel or dotted. Its omitted-facet distances are 7, 7, and 8, all larger than √21, so those relations are dotted. The radius-two central cap is also dotted to these vertex caps since `196>(2+√21)²`. The checker verifies these scalar inequalities. They identify a coherent starting mechanism for the next round; no full new chamber certificate or volume closure is claimed here.

The concrete next mechanism should place caps at interior face centers and couple them to this new layer, testing the full product for uncovered directions. The balanced vertex set makes that global test meaningful; adding one cap at the old remote escape would not. There is no finite-volume claim, no exclusion of other cusp types, and no assertion that the displayed geometric optimum is an optimum over all facet scales or all cap profiles.

```sh
python scripts/fresh_h20_intrinsic_cusp_benchmark_20261014.py --verify
```
