# The remaining H3 boundary after the exact171 replacement

The [source](../scripts/fresh_h20_replacement_boundary_20261017.py) numerically reconstructs the T-supported H3 face of the [frozen exact171 branch](fresh_h20_replacement_hexagon_20261017.md). Its [ledger](fresh_h20_replacement_boundary_20261017.json) separates genuine hyperbolic edges from intersections outside the Klein ball. This is a numerical boundary census, not another exact ambient certificate.

The final face has15 actual facets,17 finite vertices, and31 hyperbolic edges. Of these edges,21 have two finite endpoints,9 have one finite endpoint, and1 has no finite endpoint. No polyhedral vertex is numerically near the unit sphere. There are11 distinct exposed ideal limits of the open edges; these are not claimed to be proper cusps.

The ten open edges are

```
95/157, 95/171, 145/152, 145/157, 152/172,
159/160, 159/161, 160/161, 161/171, 171/172.
```

Only161/171 has no finite endpoint. The other nine each have one finite endpoint. In the graph whose nodes are facet labels and whose edges are these open hyperbolic edges, there is an outer six-cycle

```
95 — 157 — 145 — 152 — 172 — 171 — 95
```

joined by171/161 to the triangle159/160/161. Thus the current open-facet incidence graph is connected. This graph is not itself a proof that the ideal boundary is topologically connected. The separately recorded exterior polyhedral skeleton also has one connected component, but remains a numerical diagnostic.

Six facets are already compact:142 is a heptagon,150 a triangle,156 a quadrilateral,168 a triangle,169 a hexagon, and170 a triangle. In particular, no open edge is incident to169. Its complete hexagon is established by the exact predecessor certificate, while the other facet descriptions here are numerical. The remaining nine facets have open ends.

Every pair of positive-norm projected normals is tested. A genuine edge requires a positive-definite projected two-root Gram and a nonempty feasible line segment inside the unit ball. The entire affine line is intersected with all positive-side halfspaces, then clipped by the sphere; its resulting endpoints are checked against **all** retained halfspaces, including projected timelike ones. Coincident finite and ideal endpoints are merged numerically. The source also exhausts all three-plane intersections and checks boundedness of the positive-side Euclidean polytope with six coordinate linear programs.

The positive-side polytope has9 exterior vertices, with maximum squared radius about2.75097097. Those exterior vertices are not hyperbolic vertices. Timelike projected inequalities can change such outside-ball intersections without changing the hyperbolic face, so this ledger does not compare raw exterior counts from incompatible constraint conventions.

The deletion comparison is explicit. Removing158 alone loses its three finite corners; adding168 restores three corners on the opposite side of the159/160/161 region. The original eight open edges remain. Deleting31 changes no actual T-facet, hyperbolic edge, or vertex: its projection is future timelike. Adding169 and its three partners closes the169 hexagon and raises the finite-vertex count from10 to17, but leaves the ten open edges listed above. The original two open-facet incidence components become linked by the new161/171 edge.

For the next construction, the ledger favors a coordinated cap cluster addressing the outer cycle or the159/160/161 region, with every new facet's own boundary included in the test. Closing another isolated endpoint is insufficient evidence of whole-face progress. The paired open ends of157, adjacent to95 and145, provide a small incidence target; any proposed replacement must also account for the new facets it introduces.

The old residual frame is inherited numerically, and the five new coordinates are midpoints of the frozen exact definitions. Large inherited vectors exhibit floating-point cancellation in direct norm evaluation; the source records the absolute and coordinate-scaled residuals rather than treating those evaluations as exact norm proofs. Native norm identities and all Coxeter gates remain obligations supplied by the frozen exact171 package. Near-degeneracies use explicit numerical tolerances, and no all-parameter exclusion is inferred.

Status: **FINAL/FROZEN** after generation and fresh saved replay passed. Source syntax, direct source hashes, and report links were checked. The output is a numerical complete ledger in the stated fixed chart, with no new whole-H20 volume or vertex-count assertion.
