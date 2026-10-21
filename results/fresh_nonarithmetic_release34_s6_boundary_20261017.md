# Exact intrinsic boundary after releasing wall 34

The changed 65-wall chamber reduces the number of geodesic open ends of the actual S6 face from **19 to 17**. Its new facets 63 and 64 are compact quadrilaterals, and 65 is a compact pentagon. The whole H3 face remains noncompact: one exterior region is unchanged, and removal of wall 34 opens three formerly compact edges elsewhere.

The [source](../scripts/fresh_nonarithmetic_release34_s6_boundary_20261017.py) and [certificate](fresh_nonarithmetic_release34_s6_boundary_20261017.json) give a complete intrinsic cone and edge census. They use the exact [changed-wall construction](fresh_nonarithmetic_s6_release34_exact_20261017.json) and the frozen [preceding 63-wall chamber](fresh_nonarithmetic_escape62_20261017_cap.json). The numerical boundary diagnostic is neither an input nor an acceptance test. No full ambient 65-wall cone is constructed here.

| Intrinsic quantity | Preceding 63 walls | Changed 65 walls |
| --- | ---: | ---: |
| Actual projected side walls | 16 | 18 |
| Future timelike, redundant projections | 41 | 41 |
| Finite projective vertices | 19 | 23 |
| Exterior projective rays | 9 | 9 |
| Ideal rays / nonsimple rays / singular triples | 0 / 0 / 0 | 0 / 0 / 0 |
| Compact geodesic edges | 22 | 28 |
| Edges with one open end | 13 | 13 |
| Edges with two open ends | 3 | 2 |
| Projective edges outside hyperbolic space | 4 | 5 |
| Geodesic open ends | 19 | 17 |

Here S6={0,1,2,3,5,6}, with spherical diagram H4+A2. The actual common-wall face is three-dimensional. The ordered roots [8,30,44,52] are right to S6 and have an invertible Gram matrix. They therefore form a basis of its Lorentzian four-dimensional perpendicular space. Label 52 denotes the original normal before the change and its explicit replacement afterward. The source checks the right incidences, positive definiteness of S6, and every inverse pivot.

For any ambient normal n, its products with these four roots are its projected covector. Multiplication by the inverse four-root Gram gives the projection coordinates and its square. The projected common center is future timelike and strict against every projected wall. Each positive-square projection defines an actual side: Schur projection of an acute system is acute, and projecting the strict center once more onto that side gives a strict finite point of the facet. Every other nonwhite projection has strictly negative square and negative product with the future center, so is future timelike. Its inequality is strict on the entire future H3 and is redundant there. There are no zero or null projected normals in either case.

The auxiliary projective cone retains exactly the positive-square side inequalities. In particular, its exterior rays need not agree with the slice of a complete ambient cone: a timelike projected inequality can cut exterior projective directions while being redundant on every hyperbolic point. Omitting such inequalities preserves precisely the intrinsic hyperbolic face and all its geodesic edges.

Let s be the sum of the positive-side normals. Its product with every nonzero vector of their feasible cone is strictly negative. Indeed, all side products are nonpositive, and equality of their sum would make every side product zero; the four basis sides then force the vector to vanish. Thus the affine section s·x=−1 is bounded. The computation also certifies that s itself is future timelike. This provides a complete projective chart containing the intrinsic hyperbolic face.

The source solves every triple of side equations together with that chart equation: 560 systems before the change and 816 afterward. All four-by-four systems have certified nonzero pivots. A triple is either rejected by a strictly positive side product, or accepted with every nonincident side product strictly negative. The accepted points have strictly negative or strictly positive Lorentz square; none is ideal. Thus every accepted vertex is simple and no degeneracy is discarded. Every vertex of a bounded three-dimensional polytope has three independent active side equations, so this exhausts its vertices and hence its cone rays. The ray counts are 28 and 32. The certificates include every rejection and every accepted point, its norm, and all nonincident side signs.

Two common side equations identify each projective edge. The complete vertex list gives exactly two endpoints for each such support; every ray has three incident edges, and both Euler identities are checked. Write aa, ab, bb for the two endpoint squares and mutual product. A finite endpoint and the negative two-by-two determinant identify an actual geodesic edge. With two exterior endpoints, the positive cone meets H3 exactly when ab<0 and aa·bb−ab²<0 in these strictly separated cases. For every such edge, the checker gives the positive combination with coefficients bb−ab and aa−ab, whose square is

\[
 (aa\,bb-ab^2)(aa+bb-2ab)<0.
\]

Its product with the future chart is negative. Every remaining exterior edge has either ab>0 or positive determinant and lies outside H3. No equality or unresolved sign is interpreted as a decision.

Five geodesic edges disappear with wall 34: compact [7,34], [8,34], [12,34], and one-open [34,47], [34,52]. Ten new geodesic edges are compact: [8,65], [30,63], [30,64], [44,63], [44,65], [52,64], [52,65], [63,64], [63,65], [64,65]. The six changes on retained labels are:

| Side pair | Before | After |
| --- | --- | --- |
| [7,12] | Compact | One open end |
| [7,47] | Compact | One open end |
| [8,12] | Compact | One open end |
| [8,44] | One open end | Compact |
| [30,44] | One open end | Compact |
| [30,52] | Two open ends | One open end |

The last row compares the old wall 52 to its replacement; it does not claim that the supporting line is unchanged. The compact intrinsic side faces are now 10,11,44,55,63,64,65. Side 12 was compact before removal of 34 and is open afterward.

The five-ray exterior component with supports [7,33,49], [30,33,49], [30,49,59], [49,58,60], [49,59,60] is unchanged. Its links [30,49] and [33,49] still cross H3 with two open ends; [49,59] and [49,60] remain outside.

The former four-ray region involving [8,44,52] and [30,44,52] is replaced by four exterior rays with supports [7,12,47], [8,12,52], [12,47,52], [30,47,52]. Their internal links [12,47], [12,52], [47,52] are all outside H3. This records the relocated open region rather than treating the three new compact polygons as closure of the whole face.

This is a new primary exact boundary computation using shared algebraic geometry and outward rational interval arithmetic, not an independent arithmetic reconstruction of the new roots. The parent construction's isolated algebraic parameter and exact contact identities are inherited and source-bound. Generation and fresh saved replay both passed. The retained strict ambient escape still proves infinite volume, and the intrinsic census itself leaves 17 open ends. No finite-volume example is asserted.
