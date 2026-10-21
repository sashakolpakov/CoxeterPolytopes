# Independent review of the coupled 107-wall construction

The seven added walls form a coherent exact Coxeter extension of the frozen 100-wall chamber. One new wall cuts its old strict spacelike escape and terminates a previously open cusp edge. Six other walls close the S75 face as a compact right-angled hexagon. The enlarged ambient chamber still has infinite volume: a new literal rational strict spacelike escape survives all 107 walls.

The separate [checker](../scripts/fresh_h20_coupled_roles_20261014_exact_review.py) and [certificate](fresh_h20_coupled_roles_20261014_exact_review.json) import no arithmetic from the new producer. They reuse the earlier independent sparse-radical tower and reconstruct the physical vector recipes, rational projections and normalized Lorentz frame. New normalization radicals are handled by separately implemented outward rational square-root intervals. The [primary source](../scripts/fresh_h20_coupled_roles_20261014_exact.py) and [data](fresh_h20_coupled_roles_20261014_exact.json) are bound by hashes.

The replay checks all **721 new pairs: 133 right and 588 strictly dotted**, all seven actual new finite vertices, the common strict timelike center, and strict preservation of all three old bounded cusps. It also checks the renewed escape against every old and new wall. No old 100-wall pair census or complete old finite-vertex census is repeated.

## Escape-cutting cusp-edge wall

The saved germ has a positive-definite rank-19 support. Root 92 is multiplied by the positive factor `√2` when forming its rational Gram; this changes neither its hyperplane nor the projection subspace. The checker independently inverts that Gram and reconstructs the projection `W` of the rational vector `H0`, with `H0² = −2`.

For the old future null cusp `ν`, put `h = W·ν < 0`, `Z = −2hW + W²ν`, `R = Z − ν` and `P = Z + ν`. Direct arithmetic gives

`Z² = ν² = 0`, `Z·ν = −2h²`, `R² = 4h²`, `P² = −4h²`, and `R·P = 0`.

The 19 support pairings vanish. Every other old-wall pairing of `P` is strictly negative, and every other old-wall pairing of `R` is strictly dotted. Thus `P` is an actual finite endpoint on the same rank-19 edge as the retained bounded cusp. The old escape pairs positively with `R`. This wall has exactly 19 finite-angle neighbors, so its certified endpoint is its only finite vertex.

## Compact S75 face

The checker reconstructs the rational rank-18 support Gram and the physical orthogonal frame `C,U,V`, whose squares are `−d,s,t` with all three parameters positive. In the unit frame `C/√d,U/√s,V/√t`, the six normals have coefficients

`(√2,2,0)`, `(√2,1,√3)`, `(√2,−1,√3)`, `(√2,−2,0)`, `(√2,−1,−√3)`, `(√2,1,−√3)`.

Their squares are two; neighboring pairings are zero, next-neighbor pairings are minus four, and opposite pairings are minus six. Each adjacent intersection has future unit-time coordinate one and square `−1/3`. Its 18 support inequalities and two adjacent hexagon inequalities vanish; all other 87 inequalities are strict. The six resulting Klein vertices lie strictly inside the unit disk. Their bounded Euclidean hexagon is therefore precisely a compact hyperbolic right-angled face.

Each hexagon wall has exactly 20 finite-angle neighbors: the 18 support walls and its two adjacent hexagon walls. Those latter two are mutually dotted, so its only possible finite vertices are the two certified adjacent intersections. Together with `P`, this proves **exactly seven finite vertices involving a new wall**. The frozen lower bound of 134 old finite vertices is preserved, giving at least 141 overall; this is not a global vertex census.

All listed walls are actual facets. The independently verified common strict timelike center and pairwise acuteness permit its orthogonal projection onto any listed wall, producing a timelike point strict at every other wall. The same acute-extension argument preserves every old finite vertex; the three old proper cusp rays were also checked directly against all seven additions.

## Remaining global gap

The renewed rational witness has exact square

`38017800175279 / 500000000000 > 0`.

It lies in the future affine chart of the strict center and satisfies all 107 inequalities strictly. The segment from that center to the spacelike point crosses the light cone at a strict feasible point, producing an open ideal-boundary region. This proves ambient infinite volume.

The seven additions also leave the retained S61 plane unchanged. Its inverse support Gram is entrywise nonnegative and has diagonal `2/3` at wall 61. Acuteness therefore bounds each projected square by `q − (2/3)(r·N61)²`; strict dottedness gives a negative bound. The checker verifies the exact bound for `R` and the common-center orientation for all seven normals. Likewise `R` alone is redundant on S75: its inverse diagonal at wall 75 is five, and strict dottedness gives projected square less than `−9R²`. The compact S75 hexagon supplies the second geometric role; it does not supply global closure.

Generation and the saved independent replay pass. The primary old geometry is inherited through explicit manifests; the newly added pair, vertex, cusp and escape arithmetic is independently reconstructed. No finite-volume example or arithmeticity conclusion is asserted.
