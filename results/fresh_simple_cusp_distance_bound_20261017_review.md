# Independent geometric and source review of the distance15 bounds

**PASS: no mathematical gap identified in the stated restricted conclusions.** The argument gives dimension at most125 when every ideal vertex is simple, and ambient dimension at most127 when a compact facet exists. The latter permits nonsimple ideal vertices elsewhere. Neither condition is asserted for all finite-volume Coxeter chambers, so the unrestricted bound163 is unchanged.

This is a separate mathematical/source review by exact_hub_counts. It does not duplicate the producer's 384 distance computations, 245 classical graph checks, saved replay, or the inherited quasi-Lannér classification census.

Reviewed [source](../scripts/fresh_simple_cusp_distance_bound_20261017.py) SHA256:
`0af68dae0bc9625bb7fabf1629117bbfb480bcde002de03f73700d1e64f75f62`.
Reviewed [certificate](fresh_simple_cusp_distance_bound_20261017.json) SHA256:
`6eceec3e7999b98511f996c1315d22fae767beb1af0e10c9a219eab6c2d842ee`.
The [proof report](fresh_simple_cusp_distance_bound_20261017.md) supplies the geometric argument; its later validation-status wording is not an arithmetic input to this note.

## Ideal quadrilaterals and the finite-black radius condition

The crucial extension from compact to ideal polygons is sound. The relative interior of an actual edge contains a hyperbolic point even if its endpoints are ideal. Its supporting normals therefore have a positive-definite Gram. Thus, for a quadrilateral with white support W, W plus any one black side is positive definite.

For two opposite sides, positive definiteness of their combined block would give an actual finite common point by acute projection. A singular PSD block gives a positive affine-kernel combination. Its pairing with a strict timelike center is strictly negative, so this null vector is nonzero and future directed. It is orthogonal to the PSD block and nonpositive against all other walls by acuteness. Hence it gives a feasible ideal point in both closed opposite sides. Closed opposite edges of a genuine quadrilateral are disjoint, including when endpoints are ideal. Both alternatives are impossible, so each opposite block is indefinite.

A minimal indefinite subdiagram in that block contains both black nodes. Deleting either one leaves a principal subdiagram of W plus a single side, hence is **positive definite**. The catalogue test must therefore use pairs of finite deletions, as the source does, rather than arbitrary pairs of affine-or-finite deletions.

The frozen complete order6–10 catalogue covers those cases. Any connected diagram of order at most5 has distance at most3 to a two-node set: a shortest path to the set avoids the other marked node. The low-order argument also covers arbitrary small labels, and a dotted order-two core has radius zero. There is no unexamined rank or label tail in the radius7 conclusion.

Two opposite minimal indefinite cores cannot be disjoint and orthogonal in an index-one ambient space. An overlap gives a path of length at most14 between the two black sets; an inter-core edge gives at most7+1+7=15. A shortest such path has no black internal vertex. Its endpoints are consecutive sides of the quadrilateral, so it is contained in their actual finite or simple-ideal corner diagram.

For triangles, every black deletion lies in a PSD corner. A minimal indefinite core contains all three marks and has at most10 nodes. Its minimal connecting subtree yields at least two mark-to-mark paths avoiding the third mark, each of length at most8. This is valid for ideal corners as well.

## All-rank graph costs

The path, cycle, single-fork and double-fork formulas correctly separate path–path, leaf–path and leaf–leaf pairs. In particular, the affine-D double fork with path order m>=2 has

`p_k=max(m−k,0)+4[k<=m]+2[k=2]+4[k=m+1]`,

with total order N=m+4. At k=2 the count is exactly N; the other cases are at most N as well. The exceptional m=1 affine-D4 star is handled separately. The remaining exceptional finite/affine graphs have order at most9 and at most4N pairs in total. Finite dihedral diagrams, affine A1, and the small path-shaped F/G/H diagrams are covered.

These explicit formulas and the finite classification, rather than the checked rank interval, establish the uniform distance15 cost at most15N. Spherical components add. A simple ideal vertex has a simplex Euclidean link, hence one connected affine diagram of order n; a finite vertex has a spherical forest of order n. No order-n assertion is made for a nonsimple cusp.

## Strict thresholds and scope

For an all-simple chamber, the projective closure is a bounded simple polytope, so the inherited strict h-vector inequality applies even when the hyperbolic chamber has ideal vertices. At n=126 the scalar lower bound equals15n, but the actual lower inequality is strict. This excludes126, giving at most125. The two displayed parity polynomials are nonnegative for every t>=0, not merely for the sampled t values in the checker.

For a compact facet of dimension d=n−1, its vertices are finite and its two-faces are actual ambient two-faces. Their angle weights form a subset of ambient rank-n pair weights, so the upper budget is15(d+1). The two parity differences become strictly positive starting at d=127 and128, giving d<=126 and ambient n<=127. The proof does not assume that the induced facet is itself a Coxeter polytope.

The comparison with the stronger published simple **all-ideal** bound is correctly separated from the mixed finite/ideal hypothesis here. No novelty, current-best literature status, arithmeticity, or external refereeing claim is made by this review.
