# Independent root review: mixed dimensions eleven and twelve

18 October 2026. **PASS for the corrected weighted-density, vertex-ratio,
three-face and canonical four-face results.** No dimension exclusion follows.

The root read the complete corrected
[written proof](right_angled_mixed_weighted_triangles_research_20261018.md)
and its143-line
[exact checker](../scripts/check_right_angled_mixed_weighted_triangles_research_20261018.py).
A fresh read-only replay passed in chunk87cb24. The root independently
checked the substantive geometry and algebra described below. This review
covers that packet, not later proposed five-face equality extensions.

The truncation conventions are correct in the final certificate: b0 is
the finite vertex count plus2^(n-1) times the cusp count. It does not count
each old cusp twice. The exact h-polynomial coefficient elimination
therefore gives the corrected finite/cusp terms-462 in dimension12 and
-252 in dimension11. The earlier draft terms-463/-253 were not used in
this review. The common final ideal weight is2123/1280.

The root checked the face-incidence coefficients and the strict sign of
the residual using the external eight-dimensional facet bound82 and the
elementary minimum11 facets in dimension10. The bound82 is explicitly
stated in [Alexandrov2023, Section6](https://arxiv.org/html/2303.09533v1#S6).
The general truncation/Dehn--Sommerville framework and standard
finite/ideal links are published inputs. The new lower-dimensional
identities and eliminations are internal, rather than externally assumed.

All polygon inequalities, four minimal nonexceptional polygon types,
strict ratio eliminations and finite-corner triangle densities have the
correct directions. The final dimension12 ratio is

    (1474/21)*cusps < finite vertices < (1298/9)*cusps.

In dimension11 it is(268/7)*cusps < finite vertices < (236/3)*cusps.
The intervals are nonempty. They are necessary conditions, not new
existence claims or nonexistence theorems.

For the local geometry, the opposite wall of a one-finite/two-ideal
triangle has exactly the normal q-ni-nj. Three such walls are orthogonal.
The displayed opposite finite vertex satisfies every other wall strictly:
the three pair inequalities give a residual at most-a/2. All selected
walls have the stated orthogonality. The resulting five points lie in
the actual three-face and exhaust the known six-wall convex body.
The equality case of the universal six-facet minimum in dimension3 is
also sound: the area sum gives at least three cusps, and corner counting
gives at most three, forcing exactly the geometric triangular bipyramid.

The root checked the global three-face Euler identities, the non-bipyramid
gap and both resulting B3 bounds. Weighted triangle symmetrization in a
K4-free graph is valid: no term contains two nonadjacent supported
vertices, so transferring their total weight to one vertex preserves or
improves the objective. The resulting clique bound is n^3/27. The
derived strict bounds therefore really force a four-clique at a finite
corner in both ambient dimensions.

For the four-face, its final ideal point is the null sum of two
asymptotic outward normals. Acuteness checks every remaining inequality
at that point. The complete active-wall enumeration of the ten-wall
model gives five finite and five ideal vertices. All ten already lie
in the ambient polytope, so their convex hull proves the whole actual
four-face; no extra wall can cut a hidden portion of it.

The facet-neighbor/cusp-opposite injection recovers the universal
ten-facet minimum in dimension4. Equality forces every three-facet to
be the preceding bipyramid, and hence gives the same canonical four-face.
The gap to eleven facets, the incidence C>=11*a4-(n-3)*a3 and the
edge-deletion clique count are valid. Each canonical four-face has five
finite corners, so the factor5 is essential and correctly retained.

The standalone checker regenerates the DS rows and rational consequences,
and enumerates all20 and210 active-wall sets in dimensions3 and4. The
root reviewed its Gaussian elimination and feasibility tests. It does
not encode the universal geometric arguments, independently prove the
published facet bound, or establish a formal Lean theorem. Its status
flag about geometric review records issuance before this separate review.

The corrected packet therefore gives new internal necessary conditions
and a forced canonical four-face. The published mixed ceiling12 and
entirely ideal ceiling6 remain unchanged. No Lean or external human
review is claimed.
