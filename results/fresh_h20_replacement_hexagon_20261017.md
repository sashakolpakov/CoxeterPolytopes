# A 171-wall alternative with a compact hexagon

This construction explicitly deletes old walls31 and158 from the frozen168-wall chamber, retains the other166 roots, and adds five exact normals with labels168 through172. The resulting labels are not consecutive: there are171 actual walls. It is an alternative branch, not an extension retaining all168 roots.

The [source](../scripts/fresh_h20_replacement_hexagon_20261017.py) reconstructs the literal algebraic normals and uses outward rational intervals for strict inequalities. Its [certificate](fresh_h20_replacement_hexagon_20261017.json) passed generation and a fresh saved replay. The numerical discovery chain is separately preserved in the [old-edge profiles](fresh_h20_replaced_old_edge_profiles_20261017.md), [H4-edge continuation](fresh_h20_h4_edge_continuation_20261017.md), [two-wall replacement](fresh_h20_two_wall_replacement_cluster_20261017.md), and [hexagon closure](fresh_h20_h4_hexagon_closure_20261017.md). Those numerical calculations do not supply the exact sign proofs below.

The fixed spherical support is

\[
T=[1,2,3,4,5,6,7,8,9,10,12,13,14,15,61,87,151],
\]

of type A10+A4+A2+A1. New root169 joins61 at pi/3, so S=T+[169] has type A10+A4+A3+A1 and rank18. Its actual face is a hyperbolic plane. This nonright-T contact is the qualitative change from earlier right-support cap families.

The five definitions are sequential and exact:

| New label | Contacts with T | Three additional prescribed contacts |
|---|---|---|
|168|All right|159 right,160 pi/3,161 right; exact opposite-cap core formula|
|169|61 pi/3, all other T right|95 right,142 pi/5,152 pi/5|
|170|All right|142,152,169 all right|
|171|All right|95 pi/4,169 pi/3,161 pi/3|
|172|All right|152,169,171 all right|

Root168 is the branch p+sqrt((2-p²)/u²)u from the [exact opposite-cap theorem](fresh_h20_opposite_prism_gate_20261017.md). There p²=(-5+sqrt2)/23, u²=1/30, and p.u=0; the formula therefore has norm2. Its non-Coxeter product with the deleted158 is not treated as an admissible gate.

For each subsequent root, its T-projection is obtained from the exact rational inverse Gram of T. In T-perp, three prescribed products determine a particular vector p. Projecting a fixed old residual vector off that three-root span gives a one-dimensional kernel vector k. The two possible completions are

\[
r=r_T+p\pm\sqrt{\frac{2-r_T^2-p^2}{k^2}}\,k.
\]

The frame determinant, kernel norm, and positive radicand are certified. Rational open windows on an unused old product select the discovered branch. They impose no extra equality. Norm2 and the prescribed finite products are exact consequences of the projection construction; interval overlap is only a consistency check, not a proof that an approximate angle is exact. Every remaining retained contact must pass a strict dotted inequality.

The exact certificate contains840 new/prior gates:740 dotted,93 right,4 at pi/3,1 at pi/4, and2 at pi/5. Old/old gates are inherited from the frozen168 chamber after deletion. A single strict future timelike center is constructed as a rational linear combination of the residual old roots95,142,145,31, followed by a small rational displacement into the T halfspaces. The removed31 remains a legitimate algebraic coordinate vector, not a retained chamber wall.

For any retained normal n and strict center H, the point H-(H.n/n²)n lies on its wall, has negative norm, and pairs strictly negatively with every other normal by acuteness. Thus one common center and the acute Gram establish all171 actual facets; there is no inference from pairwise compatibility alone.

The complete side list of the S-face is95,142,152,170,171,172. Exactly those six old/new normals have positive squared norm after projection off S. Every other retained normal outside S has negative projected norm; its halfspace contains the entire relevant future hyperbolic plane because it contains the projected strict center. The six sides occur cyclically as

\[
95,142,170,152,172,171.
\]

Adjacent projected two-root Grams are positive definite; nonadjacent ones are strictly indefinite with negative product. The six20-root corners are positive definite, future, and feasible with every other retained wall strict. Their segments consequently form the bounded convex six-sided intersection. This certifies the whole hexagon, not just six isolated local corner diagrams.

Corner products are evaluated using the exact projection identities

\[
v=H-NG^{-1}h,\qquad v^2=v.H=H^2-h^TG^{-1}h,
\]

\[
v.n=H.n-(G^{-1}h)^T(N^Tn).
\]

This avoids cancellation in large physical coordinates at shallow finite corners. The stored physical-coordinate intervals are checked for consistency with those identities.

The construction also retains a strict spacelike escape E with E.H<0 and E.n<0 for every one of the171 walls. The segment from H to E reaches a future ideal point strictly inside every halfspace; a neighborhood of that ideal point remains in the chamber. Thus the chamber has infinite volume. A compact hexagon is genuine local closure, but not global finite volume. No complete ambient finite-vertex or cusp census is asserted, and no persistence claim is made for old corners or cusps incident to the two deleted walls.

Status: **FINAL/FROZEN**. Exact generation and fresh saved replay both passed with `PASS_EXACT_171_WALL_REPLACEMENT_WITH_COMPACT_HEXAGON`; syntax, source manifests, and report links were checked. A separate agent's narrow mathematical source read found no issue in the common-center, actual-facet, exact-angle, and complete-hexagon arguments. That read did not independently replay the840 contacts or the interval arithmetic and is not described as an independent arithmetic certificate.
